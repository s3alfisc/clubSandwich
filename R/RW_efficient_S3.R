
Wild_bootstrap <- function(obj, constraints,  
                                 cluster, bootstraps = 2000,
                                 auxilliary_dist = r_Rademacher, 
                                 residual_adjustment = "CR2", 
                                 test_adjustment = "CR0", 
                                 random_seed = 2252,
                                 ...) {
  
  # fit expanded model by updating obj according to constraints
  
  ### Data needed to generate bootstrap samples
  
  cluster <- droplevels(as.factor(cluster))
  
  T_matrx <- model_matrix(obj)
  T_list <-  matrix_list(T_matrx, cluster, "row")
  
  U_matrx <- model.matrix(data = obj$data,object =  constraints) # question: add an argument called "data"?
  U_list <- lapply(split(U_matrx, cluster),matrix, ncol = ncol(U_matrx))

  UseWeight <- T
  if(UseWeight == T){
    W_list <- weightMatrix(obj, cluster)
    W_vector <- unlist(lapply(X = W_list,FUN = diag))
    W_matrx <- diag(x = W_vector, nrow = length(W_vector))
  }else{
    # W_vector <- rep(x = 1,nrow(X_matrx)) # No weights. All 1s.
    # W_matrx <- diag(x = W_vector, nrow = length(W_vector))
    # W_list <- lapply(split(W_vector, cluster), FUN = function(x){diag(x, nrow = length(x))})
  }
  

  J <- nlevels(cluster)# Number of groups
  
  
  ### Calculate alpha_tilde
  # alpha_tilde is the estimation for alpha under the null hypothesis
  
  alpha_tilde <- na.omit(coef_CS(obj))
  
  
  ###Calculate e_tilde
  #e_tilde is the residual under the null hypothesis. (Residual of Y on T)
  
  e_tilde_matrx <- as.matrix(residuals_CS(obj))
  e_tilde_list <- matrix_list(e_tilde_matrx, cluster, "row")
  
  ###Calculate U_dd_matrx. 
  
  sumTWT <- Reduce("+", Map(function(X,W) t(X) %*% W %*% X,
                            T_list, W_list))
  M_T_matrx <- solve(sumTWT)
  
  H_T_matrx <- T_matrx %*% M_T_matrx %*% t(T_matrx) %*% W_matrx # H = X * Q * X' * W
  ImH_T_matrx <- diag(rep(x = 1,nrow(T_matrx))) - H_T_matrx
  
  #U_dd_matrx is the residuals from the WLS regression of U on T
  
  U_dd_matrx <- ImH_T_matrx %*% U_matrx
  U_dd_list <- lapply(split(U_dd_matrx, cluster),matrix, ncol = ncol(U_dd_matrx))
  
  ####FULL model
  
  
  ###Calculate beta_hat. 
  #Beta_hat is the WLS estimation for beta under the MAIN/FULL model.
  
  sumU_ddWU_dd <- Reduce("+", Map(function(U_dd,W) t(U_dd) %*% W %*% U_dd,
                                  U_dd_list, W_list))
  
  M_Udd_matrx <- solve(sumU_ddWU_dd)
  
  sumU_ddWe_tilde <- Reduce("+", Map(function(U_dd_j,W_j,e_tilde_j) t(U_dd_j) %*% W_j %*% e_tilde_j,
                                     U_dd_list ,W_list, e_tilde_list))
  
  beta_hat <- M_Udd_matrx %*%sumU_ddWe_tilde
  
  ###Calculate e_hat.
  #e_hat is the residual of the MAIN/FULL model.
  e_hat <- e_tilde_matrx - U_dd_matrx %*%beta_hat
  
  ###Full model
  
  obj_fullmodel <- list(study_orig_id = obj$study_orig_id,
                        b.r = T,
                        data.full = list(e.r = e_hat,
                                         avg.var.eff.size = obj$data.full$avg.var.eff.size,
                                         r.weights = obj$data.full$r.weights,
                                         userweights = obj$data.full$userweights),
                        Xreg = cbind(U_matrx,T_matrx),
                        user_weighting = obj$user_weighting,
                        N = obj$N
  )
  
  class(obj_fullmodel) <- "robu"
  
  vcov_full <- vcovCR(obj = obj_fullmodel, cluster = cluster, type = test_adjustment)
  ###Calculate Adjustment matrix for estimation of VCOV matrx.

  A_value <- attr(vcov_full, "adjustments")
  
  A_list <- switch(test_adjustment,
                   CR0 = lapply(split(rep(1, length(cluster)), cluster),diag),
                   CR1 =  lapply(split(rep(A_value, length(cluster)), cluster),diag),
                   CR1S = lapply(split(rep(A_value, length(cluster)), cluster),diag),
                   CR2 = A_value,
                   CR3 = A_value,
                   CR4 = A_value)
  
  
  # calculate Wald test for updated model
  vcov_full_beta <- vcov_full[1: ncol(U_matrx),1: ncol(U_matrx)]

  Qstat <- as.numeric(t(beta_hat) %*% chol2inv(chol(vcov_full_beta)) %*% beta_hat)
  
  
  #Bootstrap preparation
  ###Calculate adjustment matrix B_j. This keeps constant across all of the bootstrap sets.
  
  vcov_null <- vcovCR(obj = obj, cluster = cluster, type = residual_adjustment)
  
  B_value <- attr(vcov_null, "adjustments")
  
  B_list <- switch(test_adjustment,
                   CR0 = lapply(split(rep(1, length(cluster)), cluster),diag),
                   CR1 =  lapply(split(rep(B_value, length(cluster)), cluster),diag),
                   CR1S = lapply(split(rep(B_value, length(cluster)), cluster),diag),
                   CR2 = B_value,
                   CR3 = B_value,
                   CR4 = B_value)
  
  
  ###Calculate corrected bootstrap residual f_j = B_j * e_tilde_j
  
  f_list <- Map(f = function(B_j, e_tilde_j) {B_j %*% e_tilde_j},
                B_list , e_tilde_list)

  ####E_j =  t(U_dd_j) %*% W_j
  E_list <-  Map(f = function(U_dd_j, W_j){t(U_dd_j) %*% W_j},
                 U_dd_list, W_list)
  
  ####G_j =  E_j %*% A_j
  G_list <- Map(f = function(E_j, A_j){ E_j %*% A_j},
                E_list , A_list)
  

  ##Here is the function that generate bootstrap Q_statistics.
  Bootfun_Get_Qstat <- function(){
    
    ##Get bootstrap residuals
    ###Step one: randomly generate auxilliary random variables 
    eta_vec <- auxilliary_dist(n = J)
    eta_list <- as.list(eta_vec)
    
    ####Step two: Get bootstrap residual  e_b_j = f_j * eta_j
    e_b_vec <- unlist(Map(f = function(residual, randomVec){residual*randomVec},
                          f_list, eta_list))
    
    ##Get coefficient estimations for each bootstrap sample set.
    e_b_list <- lapply(split(e_b_vec, cluster),matrix, ncol = 1)
    sumEe_b <- Reduce("+", Map(function(E_j,e_b_j) E_j %*% e_b_j,
                               E_list, e_b_list))
    beta_hat_b <- M_Udd_matrx %*% sumEe_b
    
    ##Calculate bootstrap residuals based on the full model.
    e_hat_b_vec <- ImH_T_matrx %*% e_b_vec - U_dd_matrx %*% beta_hat_b
    
    ##Calculate bootstrap vcov matrix estimations.
    e_hat_b_list <- lapply(split(e_hat_b_vec, cluster), matrix, ncol = 1)
    sumGeeG <- Reduce("+", Map(function(G_j,e_hat_b_j) tcrossprod(G_j %*% e_hat_b_j),
                               G_list, e_hat_b_list))
    
    vcov_b <- M_Udd_matrx %*% sumGeeG %*% M_Udd_matrx
    
    ##Calculate Q statistic for each bootstrap sample
    Qstat_b<- t(beta_hat_b)%*%chol2inv(chol(vcov_b))%*%beta_hat_b
    return(Qstat_b)
  }
  
  if (!is.null(random_seed)) set.seed(random_seed)
  # Run the bootfun multiple times.
  Qstat_b_list <- replicate(n = bootstraps, 
                            expr = Bootfun_Get_Qstat(),
                            simplify = F)
  Qstat_b_vec <- unlist(Qstat_b_list)
  
  
  ####Result Analysis
  #95% quantile.
  Significance_limit <- sort(Qstat_b_vec)[floor(bootstraps*0.95)]
  
  Boot_p_val <- sum(Qstat_b_vec > Qstat)/bootstraps
  
  #Determine if its significant.
  Significance <- ifelse(Qstat>Significance_limit,"TURE","FALSE")
  
  res_df <- data.frame(Qstat = Qstat,
                     Boot_p_val = Boot_p_val,
                     Significance_limit = Significance_limit,
                     Significance = Significance)

  return(list(Boot_Qf_vec = Qstat_b_vec,
              res_df= res_df))

}

r_Mammen <- function(n) {
  pts <- c(-(sqrt(5) - 1), sqrt(5) + 1) / 2
  prob <- (sqrt(5) + 1) / (2 * sqrt(5))
  sample(x = pts, size = n, replace = TRUE, prob = c(prob, 1 - prob))
}

r_Rademacher <- function(n) {
  sample(x = c(-1L, 1L), size = n, replace = TRUE)
}

r_sixpoint <- function(n) {
  six_points <- c(-sqrt(3/2), -1, -sqrt(1 / 2), sqrt(1 / 2), 1, sqrt(3 / 2))
  sample(x = six_points, size = n, replace = TRUE)
}