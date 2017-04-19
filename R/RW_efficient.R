#----------------------------------------------
# auxilliary distribution functions
#----------------------------------------------

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

#----------------------------------------------
# wild bootstrap function
#----------------------------------------------

#' Cluster-wild bootstrap F test
#' 
#' \code{vcovCR} returns a sandwich estimate of the variance-covariance matrix 
#' of a set of regression coefficient estimates.
#' 
#' @param obj Fitted null model from which to calculate residuals for use in the
#'   cluster-wild bootstrap algorithm.
#' @param constraints formula specifying the set of additional parameter 
#'   restrictions to be tested.
#' @param bootstraps Integer specifying number of bootstrap replicates to 
#'   generate.
#' @param auxilliary_dist Function for generating auxilliary random variables 
#'   for perturbing the residuals in the cluster wild bootstrap algorithm. The 
#'   default is \code{r_Rademacher}, which samples from a two-point distribution 
#'   with equal mass at -1 and 1. A user-specified function may be passed 
#'   instead.
#' @param clusters Integer vector indicates the group id of each observation.
#' @param target Optional matrix or vector describing the working 
#'   variance-covariance model used to calculate the \code{CR2} and \code{CR4} 
#'   adjustment matrices. If a vector, the target matrix is assumed to be 
#'   diagonal. If not specified, \code{vcovCR} will attempt to infer a value.
#' @residual_adjustment Character string specifying which small-sample 
#'   adjustment to use when estimating the residuals under the null model. 
#'   Available options are the same as the \code{type} argument of 
#'   \code{\link{vcovCR}}. Default is \code{"CR2"}.
#' @test_adjustment. Character string specifying which small-sample adjustment 
#' to use when calculating the Wald test statistic. Available options are the
#' same as the \code{type} argument of \code{\link{vcovCR}}. Default is
#' \code{"CR0"}.
#' @param ... Additional arguments passed to \code{\link{vcovCR}} when 
#'   calculating the adjustment matrices used to estimate the residuals.
#' @inheritParams vcovCR
#'   
#' @description
#' @details
#' @references
#' @return
#' @examples 
#' 
#' 
#' @export


Wild_bootstrap <- function(obj, constraints,  
                           cluster, bootstraps = 2000,
                           auxilliary_dist = r_Rademacher, 
                           residual_adjustment = "CR2", 
                           test_adjustment = "CR0", 
                           random_seed = 2252,
                           ...) {
  
  # fit expanded model by updating obj according to constraints
  
  ### Data needed to generate bootstrap samples
  Study_id <- obj$X.full[,1]
  cluster <- Study_id
  
  T_matrx <- obj$Xreg
  T_list <- lapply(split(T_matrx, Study_id),matrix, ncol = ncol(T_matrx))

  U_matrx <- model.matrix(data = obj$data,object =  constraints)
  U_list <- lapply(split(U_matrx, Study_id),matrix, ncol = ncol(U_matrx))
  
  Y_matrx <- obj$data.full$effect.size
  Y_list <- lapply(split(Y_matrx, Study_id),matrix, ncol = 1)
  
  UseWeight <- T
  if(UseWeight == T){
    W_vector<- obj$data.full$r.weights
  }else{
    W_vector <- rep(x = 1,nrow(X_matrx)) # No weights. All 1s.
    
  }

  W_matrx <- diag(x = W_vector, nrow = length(W_vector))
  W_list <- lapply(split(W_vector, Study_id), FUN = function(x){diag(x, nrow = length(x))})
  
  m <- length(Y_list)# Number of groups
  
  
  ### Calculate alpha_tilde
  # alpha_tilde is the estimation for alpha under the null hypothesis
  #Actually:
  #alpha_tilde <- obj$b.r
  sumTWT <- Reduce("+", Map(function(X,W) t(X) %*% W %*% X,
                             T_list, W_list))
  M_T_matrx <- solve(sumTWT)
  sumTWy <- Reduce("+", Map(function(T_j,W_j,y_j) t(T_j) %*% W_j %*% y_j,
                            T_list ,W_list, Y_list))

  alpha_tilde <-  M_T_matrx %*% sumTWy
  
  
  ###Calculate e_tilde
  #e_tilde is the residual under the null hypothesis. (Residual of Y on T)
  
  H_T_matrx <- T_matrx %*% M_T_matrx %*% t(T_matrx) %*% W_matrx # H = X * Q * X' * W
  ImH_T_matrx <- diag(rep(x = 1,nrow(T_matrx))) - H_T_matrx
  
  e_tilde_matrx <- ImH_T_matrx %*% Y_matrx
  #Or:
  #e_tilde_matrx <- Y_matrx - T_matrx %*% alpha_tilde
  #Or:
  #e_tilde_matrx <- obj$data.full$e.r
  e_tilde_list <- lapply(split(e_tilde_matrx, Study_id), matrix, ncol = 1)
  
  ###Calculate U_dd_matrx. 
  #U_dd_matrx is the residuals from the WLS regression of U on T
  
  U_dd_matrx <- ImH_T_matrx %*% U_matrx
  U_dd_list <- lapply(split(U_dd_matrx, Study_id),matrix, ncol = ncol(U_dd_matrx))
  
  
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
  
  if (test_adjustment %in% c("CR0","CR1","CR1S")) {
    if(test_adjustment == "CR0"){
      #All identity matrices
      A_list <-  lapply(split(rep(1, length(Study_id)), Study_id),diag)
    }else{
      ##Identity matrices * A_value
      A_value <- attr(vcov_full, "adjustments")
      A_list <-  lapply(split(rep(A_value, length(Study_id)), Study_id),diag)
    }
  }else{
    
    A_list <- attr(vcov_full, "adjustments")
  }
  
  


  # calculate Wald test for updated model
  vcov_full_beta <- vcov_full[1: ncol(U_matrx),1: ncol(U_matrx)]
  
  
  Qstat <- as.numeric(t(beta_hat) %*% chol2inv(chol(vcov_full_beta)) %*% beta_hat)
  
  
  
  # prepare for bootstrapping
  ###Calculate adjustment matrix B_j. This keeps constant across all of the bootstrap sets.
  
  vcov_null <- vcovCR(obj = obj, cluster = cluster, type = residual_adjustment)
  
  
  if (residual_adjustment %in% c("CR0","CR1","CR1S")) {
    if(residual_adjustment == "CR0"){
      #All identity matrices
      B_list <-  lapply(split(rep(1, length(Study_id)), Study_id),diag)
    }else{
      ##Identity matrices * A_value
      B_value <- attr(vcov_null, "adjustments")
      B_list <-  lapply(split(rep(B_value, length(Study_id)), Study_id),diag)
    }
  }else{
    
    B_list <- attr(vcov_null, "adjustments")
  }
  
  
  ###Calculate corrected bootstrap residual f_j = B_j * e_tilde_j
  
  f_list <- Map(f = function(B_j, e_tilde_j) {B_j %*% e_tilde_j},
                B_list , e_tilde_list)
  
  ###Generate bootstrap samples
  #res_list_null <- e_tilde_list
  

  Get_e_b_vec <- function(){
    #Step one: randomly generate auxilliary random variables 
    eta_vec <- auxilliary_dist(n = m)
    eta_list <- as.list(eta_vec)
    
    #Step two: Get bootstrap residual  e_b_j = f_j * eta_j
    e_b_vec <- unlist(Map(f = function(residual, randomVec){residual*randomVec},
                         f_list, eta_list))

    return(e_b_vec)
  }
  
  set.seed(random_seed)
  #Get Bootstrap data list: a list of e_b_vec, bootstrap residuals.
  e_b_vec_list <- replicate(n = bootstraps, 
                           expr = Get_e_b_vec(),
                           simplify = F)
  
  
  ###Calculate bootstrap Beta estimations
  
  ####E_j =  t(U_dd_j) %*% W_j
  E_list <-  Map(f = function(U_dd_j, W_j){t(U_dd_j) %*% W_j},
                 U_dd_list, W_list)
    
  ####G_j =  E_j %*% A_j
  G_list <- Map(f = function(E_j, A_j){ E_j %*% A_j},
                E_list , A_list)
  
  Get_beta_hat_b <- function(e_b_vec){
    
    e_b_list <- lapply(split(e_b_vec, Study_id),matrix, ncol = 1)
    sumEe_b <- Reduce("+", Map(function(E_j,e_b_j) E_j %*% e_b_j,
                               E_list, e_b_list))
    beta_hat_b <- M_Udd_matrx %*% sumEe_b
    
    return(beta_hat_b)
  }
  beta_hat_b_list <- Map(f = Get_beta_hat_b, e_b_vec_list)
  
  
  ###Calculate bootstrap residuals based on the full model.
  
  e_hat_b_vec_list <- Map(f = function(e_b_vec, beta_hat_b){ImH_T_matrx %*% e_b_vec - U_dd_matrx %*% beta_hat_b},
                          e_b_vec_list, beta_hat_b_list)
  
  ###Calculate bootstrap vcov matrix estimations.
  
  Get_vcov_b <- function(e_hat_b_vec){
    e_hat_b_list <- lapply(split(e_hat_b_vec, Study_id), matrix, ncol = 1)
    sumGeeG <- Reduce("+", Map(function(G_j,e_hat_b_j) tcrossprod(G_j %*% e_hat_b_j),
                                    G_list, e_hat_b_list))
    
    vcov_b <- M_Udd_matrx %*% sumGeeG %*% M_Udd_matrx
  }
  
  vcov_b_list <- Map(f = Get_vcov_b, e_hat_b_vec_list)
  
  Qstat_b_list <- Map(f = function(beta_hat_b,vcov_b){t(beta_hat_b)%*%chol2inv(chol(vcov_b))%*%beta_hat_b},
                      beta_hat_b_list, vcov_b_list)
  
  Qstat_b_vec <- unlist(Qstat_b_list)

  ####Result Analysis
  
  
  #95% quantile.
  Significance_limit <- sort(Qstat_b_vec)[floor(bootstraps*0.95)]
  
  Boot_p_val <- sum(Qstat_b_vec > Qstat)/bootstraps
  
  #Determine if its significant.
  Significance <- ifelse(Qstat>Significance_limit,"TURE","FALSE")
  
  res_vec <- round(c(Qstat = Qstat,
                     Boot_p_val = Boot_p_val,
                     Significance_limit = Significance_limit),digits = 3)
  # 
  # return(list(e_b_vec_list = e_b_vec_list[[2]],
  #             predict_mean = T_matrx%*%alpha_tilde,
  #             beta_hat_b_list = beta_hat_b_list[[2]],
  #             vcov_b_list = vcov_b_list[c(1,2,3)],
  #             vcov_full_beta = vcov_full_beta,
  #             e_hat_b_vec_list = e_hat_b_vec_list[[1]],
  #             f_list = f_list,
  #             e_tilde_list = e_tilde_list))
  # 
  return(list(Boot_Qf_vec = Qstat_b_vec,
              res_vec= c(res_vec,Significance=Significance)))

  #res <- data.frame(Wald_statistic = Q, p_val = p_val)
  #return(res)
}
