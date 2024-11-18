#' #' This is the main function. EM algorithm function and change the beta prior to find better results.
#'
#' @param r_t_pure Read counts of Input sample
#' @param r_t_mix Read counts of IP sample
#' @param beta the hyper parameter of Beta prior
#' @param alpha_a the initial alpha_a
#' @param alpha_b the initial alpha_b
#' @param tau_a_vector the vector of the initial tau_a
#' @param l the effective length of site t
#' @param prop_vector the vector of the prop of the beta hyper parameter
#' @param max_iteration max iteration
#' @param stop_stepsize the threshold for stop
#'
#'
#' @return a list with name c("tau_a_bayes","tau_b_bayes","alpha_a_bayes","alpha_b_bayes",
#' "tau_a_no_bayes","tau_b_no_bayes","alpha_a_no_bayes",
#' "alpha_b_no_bayes","num_step",
#' "difference_value_tau","log_likelihood_bayes",
#' "log_likelihood_no_bayes","tau_a_bayes_vector",
#' "tau_a_no_bayes_vector","q_matrix_mix_a","q_matrix_mix_b"
#' ,"beta_a","beta_b","lm_res","tau_a_bayes_vector_aeeip","log_likelihood_bayes_vector",
#' tau_a_no_bayes_vector_aeeip","log_likelihood_no_bayes_vector","lm_tau_vector")
#' @export
#' @examples
#' set.seed(1)
#' data <-simulateData(T=1000, hyper_a_1=20,hyper_b_1 =0.2, hyper_a_2=50,hyper_b_2 =0.05,
#' pure_rate=0.99,mix_rate=0.1,mode=2,replicates=3)
#'
#' res <- AEEIP_test(r_t_pure=data[[3]][[1]],r_t_mix=data[[4]][[1]],beta=2,alpha_a=NA,
#' alpha_b=NA,tau_a_vector=seq(0,0.99,by=0.01),l=1,prop_vector=0.5,
#' max_iteration =1e6,stop_stepsize=1e-7)
#' res <- AEEIP_test(r_t_pure=data[[3]][[1]],r_t_mix=data[[4]][[1]],beta=100,alpha_a=NA,
#' alpha_b=NA,tau_a_vector=1-1e-5,l=1,prop_vector=seq(0.01,0.99,by=0.01),
#' max_iteration =1e6,stop_stepsize=1e-7)
AEEIP_test <- function(r_t_pure,r_t_mix,beta,alpha_a=NA,alpha_b=NA,tau_a_vector,l,
                       prop_vector=seq(0.01,0.99,by=0.01),
                       max_iteration =1e6,stop_stepsize=1e-7){
  options(warn=-1)


  if(length(tau_a_vector)>1&length(prop_vector)>1){
    stop("tau_a_vector and prop_vector can not be vector simultaneously")
  }

  if(length(tau_a_vector)*length(prop_vector)==1){
    stop("tau_a_vector and prop_vector are not vector, you can use EM_function")
  }
  if(is.na(alpha_a)){
    #temp <- rexp(length(r_t_pure))
    #alpha_a <- temp/sum(temp)
    alpha_a <- rep(1/length(r_t_pure),length(r_t_pure))
  }

  if(is.na(alpha_b)){
    #temp <- rexp(length(r_t_pure))
    #alpha_a <- temp/sum(temp)
    alpha_b <- rep(1/length(r_t_pure),length(r_t_pure))
  }



  if(length(tau_a_vector)>1){
    tau_a_bayes_vector_aeeip  <-log_likelihood_bayes_vector<-tau_a_no_bayes_vector_aeeip <- c()
    log_likelihood_no_bayes_vector<- lm_tau_vector<-lm_tau_vector2<- c()
    for( i in 1:length(tau_a_vector)){
      res <- EM_function(r_t_pure,r_t_mix,beta,alpha_a=alpha_a,alpha_b=alpha_b,
                         tau_a=tau_a_vector[i],l,prop=prop_vector,
                         max_iteration =max_iteration,stop_stepsize=stop_stepsize)
      tau_a_bayes_vector_aeeip [i] <- res$tau_a_bayes
      log_likelihood_bayes_vector [i]<- res$log_likelihood_bayes
      tau_a_no_bayes_vector_aeeip [i] <- res$tau_a_no_bayes
      log_likelihood_no_bayes_vector [i]<- res$log_likelihood_no_bayes
      lm_tau_vector [i] <- res[["lm_res"]][["coefficients"]]
    }




    res <- EM_function(r_t_pure,r_t_mix,beta,alpha_a=alpha_a,alpha_b=alpha_b,tau_a_vector[
      which.max(tau_a_no_bayes_vector_aeeip)],l,prop=prop_vector,
      max_iteration =max_iteration,stop_stepsize=stop_stepsize)

    res$tau_a_bayes_vector_aeeip <- tau_a_bayes_vector_aeeip
    res$log_likelihood_bayes_vector <-log_likelihood_bayes_vector
    res$tau_a_no_bayes_vector_aeeip <- tau_a_no_bayes_vector_aeeip
    res$log_likelihood_no_bayes_vector <- log_likelihood_no_bayes_vector
    res$lm_tau_vector <- lm_tau_vector
  }

  if(length(prop_vector)>1){
    tau_a_bayes_vector_aeeip  <-log_likelihood_bayes_vector<-tau_a_no_bayes_vector_aeeip <- c()
    log_likelihood_no_bayes_vector<- lm_tau_vector<- lm_tau_vector2 <- c()
    for( i in 1:length(prop_vector)){
      res <- EM_function(r_t_pure,r_t_mix,beta,alpha_a=alpha_a,alpha_b=alpha_b,tau_a=tau_a_vector,l,prop=prop_vector[i],
                         max_iteration =max_iteration,stop_stepsize=stop_stepsize)
      tau_a_bayes_vector_aeeip [i] <- res$tau_a_bayes
      log_likelihood_bayes_vector [i]<- res$log_likelihood_bayes
      tau_a_no_bayes_vector_aeeip [i] <- res$tau_a_no_bayes
      log_likelihood_no_bayes_vector [i]<- res$log_likelihood_no_bayes
      lm_tau_vector [i] <- res[["lm_res"]][["coefficients"]]
    }




    res <- EM_function(r_t_pure,r_t_mix,beta,alpha_a=alpha_a,alpha_b=alpha_b,tau_a=tau_a_vector,l,prop=prop_vector[
      which.max(tau_a_no_bayes_vector_aeeip)],
      max_iteration =max_iteration,stop_stepsize=stop_stepsize)

    res$tau_a_bayes_vector_aeeip <- tau_a_bayes_vector_aeeip
    res$log_likelihood_bayes_vector <-log_likelihood_bayes_vector
    res$tau_a_no_bayes_vector_aeeip <- tau_a_no_bayes_vector_aeeip
    res$log_likelihood_no_bayes_vector <- log_likelihood_no_bayes_vector
    res$lm_tau_vector <- lm_tau_vector
  }





  res
}
