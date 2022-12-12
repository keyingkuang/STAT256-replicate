### Set Up ##############################################################

library(data.table)
library(Rlab)
library(LaplacesDemon)
library(sl3)
library(tmle3)
library(tmle3shift)
library(haldensify)
library(Rsolnp)
library(xgboost)
library(resample)

# global parameter initialization
true_psi <- 22.95

# learners/specs used in the TMLE and IPTW analysis

haldensify_lrnr <- Lrnr_haldensify$new(
  n_bins = 5, grid_type = "equal_mass",
  lambda_seq = exp(seq(-1, -9, length = 100))
)
hse_lrnr <- Lrnr_density_semiparametric$new(mean_learner = Lrnr_glm$new())
mvd_lrnr <- Lrnr_density_semiparametric$new(mean_learner = Lrnr_glm$new(),
                                            var_learner = Lrnr_mean$new())
g_learner <- Lrnr_sl$new(
  learners = list(haldensify_lrnr, hse_lrnr, mvd_lrnr),
  metalearner = Lrnr_solnp_density$new()
)

mean_lrnr <- Lrnr_mean$new()
glm_lrnr <- Lrnr_glm$new()
xgb_lrnr <- Lrnr_xgboost$new()
Q_learner <- Lrnr_sl$new(
  learners = list(mean_lrnr, glm_lrnr, xgb_lrnr),
  metalearner = Lrnr_nnls$new()
)

learner_list <- list(Y = Q_learner, A = g_learner)

### Simulation  #########################################################
#'@param n_obs the number of observations simulated each time, 50, 100, 150, 200
#'@param delta the shifted function predefined
#'@param case one of the four cases in the paper (1)g&Q correct (2)g correct (3)Q correct (4)both incorrect

simulation1 <- function(n_obs = 50, delta = 2, case = 1, n_boot = 50){
  result_mat <- matrix(NA, n_boot, 3)
  coverage_mat <- matrix(NA, n_boot, 3)
  colnames(result_mat) <- c("TMLE", "IPTW", "A-IPTW")
  colnames(coverage_mat) <- c("TMLE", "IPTW", "A-IPTW")
  for(i in 1:n_boot){
    W1 <- as.numeric(runif(n_obs, 0 , 1))
    W2 <- as.numeric(rbern(n_obs, 0.7))
    
    A <- as.numeric(rpois(n_obs, lambda = exp(3+0.3*log(W1) - 0.2*exp(W1)*W2)))
    
    Y_mean <- 1+0.5*A - 0.2*A*W2 + 2*A*tan(W1^2) - 2*W1*W2 + A*W1*W2
    Y_cov <- diag(rep(1, n_obs))
    Y <- rmvn(1, Y_mean, Y_cov)
    Y <- as.vector(Y)
    
    data <- data.table(W1, W2, A, Y)
    node_list <- list(W = c("W1", "W2"), A = "A", Y = "Y")
    
    # TMLE
    tmle_spec <- tmle_shift(shift_val = delta,
                            shift_fxn = shift_additive,
                            shift_fxn_inv = shift_additive_inv)
    tmle <- tmle3(tmle_spec, data, node_list, learner_list)
    
    tmle_est <- tmle$estimates[[1]]$psi
    tmle_covered <- true_psi >= tmle$summary$lower && true_psi <= tmle$summary$upper
    result_mat[i,1] <- tmle_est
    coverage_mat[i,1] <- tmle_covered
    
    # IPTW
    
    # building the learning task
    covars <- c("W1","W2")
    IPTW_task <- sl3_Task$new(
      data = data,
      covariates = covars,
      outcome = "A"
    )
    
    # inherit the same super learner (g_learner) as TMLE to compare
    
    IPTW_learner <- g_learner
    g_fit <- hse_lrnr$train(IPTW_task)
    g_preds <- g_fit$predict()
    
    g0 <- g_preds
    data_shifted <- copy(data)
    data_shifted$A = data_shifted$A - 2
    shifted_task <- make_sl3_Task(data_shifted, covariates = covars, outcome = "A")
    
    g_star <- g_fit$predict(shifted_task)
    
    iptw_est <- mean(g_star/g0*data$Y)
    
    # iptw is asymptotically linear estimator with IC D_IPTW
    # derive the variance 
    iptw_sd <- sd(g_star/g0*Y - iptw_est)
    iptw_covered <- true_psi >= iptw_est-1.96*iptw_sd/sqrt(n_obs) && true_psi <= iptw_est+1.96*iptw_sd/sqrt(n_obs)
    result_mat[i,2] <- iptw_est
    coverage_mat[i,2] <- iptw_covered
    
    ### A-IPTW #####################################
    data_upshift <- copy(data)
    data_upshift$A <- data_upshift$A + 2
    task_Q_upshift <- make_sl3_Task(
      data_upshift,
      covariates = c(covars, "A"),
      outcome = "Y"
    )
    
    task_Q <- make_sl3_Task(
      data,
      covariates = c(covars, "A"),
      outcome = "Y"
    )
    fit_Q <- Q_learner$train(task_Q)
    pred_Q <- fit_Q$predict()
    pred_Q_upshift <- fit_Q$predict(task_Q_upshift)
    
    aiptw_est <- mean(g_star / g0 * (data$Y - pred_Q) + pred_Q_upshift)
    
    # estimate the variance using influence curve (estimator as a function of empirical mean)
    aiptw_sd <- sd(g_star/g0* (data$Y - pred_Q) + pred_Q_upshift - aiptw_est)
    aiptw_covered <- true_psi >= aiptw_est-1.96*aiptw_sd/sqrt(n_obs) && true_psi <= aiptw_est+1.96*aiptw_sd/sqrt(n_obs)
    result_mat[i,3] <- aiptw_est
    coverage_mat[i,3] <- aiptw_covered

  }
  return(list(result_mat, coverage_mat))
}

result50 = simulation1(n_boot = 50)
result100 = simulation1(n_boot = 100)
result200 = simulation1(n_boot = 200)
result500 = simulation1(n_boot = 500)


write.csv(result50[[1]], file="~/STAT256/output/est50.csv")
write.csv(result50[[2]], file="~/STAT256/output/cov50.csv")
write.csv(result100[[1]], file="~/STAT256/output/est100.csv")
write.csv(result100[[2]], file="~/STAT256/output/cov100.csv")
write.csv(result200[[1]], file="~/STAT256/output/est200.csv")
write.csv(result200[[2]], file="~/STAT256/output/cov200.csv")
write.csv(result500[[1]], file="~/STAT256/output/est500.csv")
write.csv(result500[[2]], file="~/STAT256/output/cov500.csv")


cov_mat = matrix(NA, 4, 3)
colnames(cov_mat) = c("TMLE", "IPTW", "A-IPTW")
rownames(cov_mat) = c("50", "100", "200", "500")
cov_mat[1,] = colMeans(result50[[2]])
cov_mat[2,] = colMeans(result100[[2]])
cov_mat[3,] = colMeans(result200[[2]])
cov_mat[4,] = colMeans(result500[[2]])

est_mat = matrix(NA, 4, 3)
colnames(est_mat) = c("TMLE", "IPTW", "A-IPTW")
rownames(est_mat) = c("50", "100", "200", "500")
est_mat[1,] = colMeans(result50[[1]])
est_mat[2,] = colMeans(result100[[1]])
est_mat[3,] = colMeans(result200[[1]])
est_mat[4,] = colMeans(result500[[1]])

var_mat = matrix(NA, 4, 3)
colnames(var_mat) = c("TMLE", "IPTW", "A-IPTW")
rownames(var_mat) = c("50", "100", "200", "500")
var_mat[1,] = colVals(result50[[1]])
var_mat[2,] = colVars(result100[[1]])
var_mat[3,] = colVars(result200[[1]])
var_mat[4,] = colVars(result500[[1]])

write.csv(est_mat, file="~/STAT256/output/est.csv")
write.csv(cov_mat, file="~/STAT256/output/cov.csv")
write.csv(var_mat, file="~/STAT256/output/var.csv")

