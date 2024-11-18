#' EM algorithm function
#'
#' @param r_t_pure Read counts of Input sample
#' @param r_t_mix Read counts of IP sample
#' @param beta the hyper parameter of Beta prior
#' @param alpha_a the initial alpha_a
#' @param alpha_b the initial alpha_b
#' @param tau_a the initial tau_a
#' @param l the effective length of site t
#' @param prop prop of beta
#' @param max_iteration max iteration
#' @param stop_stepsize the threshold for stop
#'
#' @return a list with name c("tau_a_bayes","tau_b_bayes","alpha_a_bayes","alpha_b_bayes",
#' "tau_a_no_bayes","tau_b_no_bayes","alpha_a_no_bayes",
#' "alpha_b_no_bayes","num_step",
#' "difference_value_tau","log_likelihood_bayes",
#' "log_likelihood_no_bayes","tau_a_bayes_vector",
#' "tau_a_no_bayes_vector","q_matrix_mix_a","q_matrix_mix_b"
#' ,"beta_a","beta_b","lm_res")
#' @export
#'
#' @importFrom stats lm
#' @examples
#'
#' set.seed(1)
#' data <-simulateData(T=1000, hyper_a_1=20,hyper_b_1 =0.2, hyper_a_2=50,hyper_b_2 =0.05,
#' pure_rate=0.99,mix_rate=0.1,mode=2,replicates=3)
#'
#' res <- EM_function(r_t_pure=data[[3]][[1]],r_t_mix=data[[4]][[1]],2,
#' tau_a=0.9,l=1,prop=0.5)
#'
#'
EM_function <- function(r_t_pure,r_t_mix,beta,alpha_a=NA,alpha_b=NA,tau_a,l,prop=0.5,
                        max_iteration =1e6,stop_stepsize=1e-7){

  options(warn=-1)
  tau_b <- 1-tau_a
  if (length(r_t_mix)!=length(r_t_pure)){
    stop("the length of mix read and pure read are not same")
  }

  if(anyNA(r_t_pure)|anyNA(r_t_mix)){
    stop("r_t_pure or r_t_mix contain NA")
  }

  if(length(alpha_a)==1){
    #temp <- rexp(length(r_t_pure))
    #alpha_a <- temp/sum(temp)
    alpha_a <- rep(1/length(r_t_pure),length(r_t_pure))
  }

  if(length(alpha_b)==1){
    #temp <- rexp(length(r_t_pure))
    #alpha_b <- temp/sum(temp)
    alpha_b <- rep(1/length(r_t_pure),length(r_t_pure))
  }



  N_pure <- sum(r_t_pure)
  N_mix <- sum(r_t_mix)
  alpha_mix_sample <- r_t_mix/N_mix
  tau_a_vector <- c()
  tau_a_no_bayes_vector <- c()
  beta_a <- prop*beta
  beta_b <-  (1-prop)*beta


  if(beta_a*beta_b<1){
    stop("beta_a or beta_b should >=1")
  }
  for (num_step in 1:max_iteration ){


    ## E- step

    q_matrix_mix_a <-  r_t_mix*(alpha_a*tau_a/l)/(alpha_a*tau_a/l+alpha_b*tau_b/l)
    q_matrix_mix_b <-  r_t_mix*(alpha_b*tau_b/l)/(alpha_a*tau_a/l+alpha_b*tau_b/l)



    q_matrix_mix_a [is.na(q_matrix_mix_a)] <- 0
    q_matrix_mix_b [is.na(q_matrix_mix_b)] <- 0

    ## M- step

    A <- sum(q_matrix_mix_a,na.rm = TRUE)+beta_a-1
    B <- sum(q_matrix_mix_b,na.rm = TRUE)+beta_b-1
    beta_sum <- beta_a+beta_b
    tau_a_new <- A/(A+B)
    tau_b_new <- 1-tau_a_new
    tau_a_vector[num_step] <- tau_a_new
    tau_a_no_bayes_vector[num_step] <- sum(q_matrix_mix_a,na.rm = TRUE)/N_mix


    alpha_a_new <- (r_t_pure+q_matrix_mix_a)/(N_pure + (tau_a_new*N_mix))
    alpha_a_new[is.na(alpha_a_new)] <- 0
    #alpha_a_new <- alpha_a_new /sum(alpha_a_new)
    #alpha_a_new <- r_t_pure/sum(r_t_pure)
    alpha_b_new <- q_matrix_mix_b/(tau_b_new*N_mix)
    alpha_b_new[is.na(alpha_b_new)] <- 0
    #alpha_b_new <- alpha_b_new /sum(alpha_b_new)
    # stop when the difference is low
    difference_value_tau <- abs(tau_a_new-tau_a)



    if (abs(tau_a_new-tau_a)<=stop_stepsize) {break}


    tau_a <- tau_a_new
    tau_b <- tau_b_new
    alpha_a <- alpha_a_new
    alpha_b <- alpha_b_new


  }




  #to calculate the likelihood
  log_alpha_a<- log(alpha_a_new/l)
  log_alpha_a[is.infinite(log_alpha_a)] <- 1e-10
  log_alpha_mix<- log(alpha_a_new*tau_a_new/l+
                        alpha_b_new*tau_b_new/l)
  log_alpha_mix[is.infinite(log_alpha_mix)] <- 1e-10

  constant <- r_t_pure*log_alpha_a+r_t_mix*log_alpha_mix
  constant <- sum(constant,na.rm = TRUE)
  log_likelihood <- constant


  alpha_a_no_bayes <-  (r_t_pure+q_matrix_mix_a)/(N_pure + (sum(q_matrix_mix_a)))
  alpha_b_no_bayes <-  (q_matrix_mix_b)/((sum(q_matrix_mix_b)))
  tau_a_no_bayes <- sum(q_matrix_mix_a)/N_mix
  tau_b_no_bayes <- 1-tau_a_no_bayes



  #to calculate the likelihood
  log_alpha_a_no_bayes<- log(alpha_a_no_bayes/l)
  log_alpha_a_no_bayes[is.infinite(log_alpha_a_no_bayes)] <- 0
  log_alpha_mix_no_bayes<- log(alpha_a_no_bayes*tau_a_no_bayes/l+
                                 alpha_b_no_bayes*tau_b_no_bayes/l)
  log_alpha_mix_no_bayes[is.infinite(log_alpha_mix_no_bayes)] <- 0

  constant_2 <- r_t_pure*log_alpha_a_no_bayes+r_t_mix*log_alpha_mix_no_bayes
  constant_2 <- sum(constant_2,na.rm = TRUE)
  log_likelihood_no_bayes  <- constant_2




  alpha_m <- r_t_mix/sum(r_t_mix)

  lm_y <- (alpha_m -alpha_b_no_bayes)
  lm_x <- (alpha_a_no_bayes-alpha_b_no_bayes)
  lm_y <- as.numeric(lm_y)
  lm_x <- as.numeric(lm_x)
  data = data.frame(cbind(lm_y,lm_x))

  lm_res <- lm(lm_y~lm_x-1,data)


  list_with_beta <- list(
    tau_a_new,tau_b_new,
    alpha_a_new,alpha_b_new,
    sum(q_matrix_mix_a)/N_mix,1-sum(q_matrix_mix_a)/N_mix,
    alpha_a_no_bayes,alpha_b_no_bayes,num_step,
    difference_value_tau,log_likelihood,
    log_likelihood_no_bayes,tau_a_vector,
    tau_a_no_bayes_vector,q_matrix_mix_a,q_matrix_mix_b
    ,beta_a,beta_b,lm_res)

  names(list_with_beta) <- c(
    "tau_a_bayes","tau_b_bayes",
    "alpha_a_bayes","alpha_b_bayes",
    "tau_a_no_bayes","tau_b_no_bayes",
    "alpha_a_no_bayes",
    "alpha_b_no_bayes","num_step",
    "difference_value_tau","log_likelihood_bayes",
    "log_likelihood_no_bayes","tau_a_bayes_vector",
    "tau_a_no_bayes_vector","q_matrix_mix_a","q_matrix_mix_b"
    ,"beta_a","beta_b","lm_res")

  list_with_beta
}

