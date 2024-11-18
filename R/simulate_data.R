#' Simulate Data function base on uniform or negative binomial distribution
#' @param T number of site
#' @param hyper_a_1 hyper parameter, minimum for uniform, mean for negative binomial for input
#' @param hyper_b_1 hyper parameter, maximum for uniform, size for negative binomial for input
#' @param hyper_a_2 hyper parameter, minimum for uniform, mean for negative binomial for IP
#' @param hyper_b_2 hyper parameter, maximum for uniform, size for negative binomial for IP
#' @param pure_rate the rate of pure sample (Input sample) from RNA fragment A
#' @param mix_rate the rate of mix sample (IP sample) from RNA fragment A
#' @param mode 1 means uniform distribution, 2 means negative binomial
#' @param replicates number of replicates
#'
#' @return list(input_real_list,IP_real_list,input_list,IP_mix_list)
#' input_real_list : real input
#' IP_real_list : real IP
#' input_list : mix input
#' IP_mix_list : mix IP
#'
#' @export
#' @importFrom stats rnbinom runif
#'
#' @examples
#'
#'
#' set.seed(1)
#' data <-simulateData(T=6000, hyper_a_1=0,hyper_b_1 =100, hyper_a_2=0,hyper_b_2 =50,
#' pure_rate=0.99,mix_rate=0.1,mode=1,replicates=3)
#'
#' set.seed(1)
#' data <-simulateData(T=6000, hyper_a_1=20,hyper_b_1 =0.2, hyper_a_2=50,hyper_b_2 =0.05,
#' pure_rate=0.99,mix_rate=0.1,mode=2,replicates=3)

simulateData <- function( T=2000, hyper_a_1=0,hyper_b_1 =1000,hyper_a_2=0,hyper_b_2 =1000,
                          pure_rate=0.999,mix_rate=0.24,mode=1,replicates=3){


  # simulate the RNA fragments for modification and non_modification
  if (mode==1){
    # follow discrete uniform distribution
    r_t_unmodification <-  as.integer(runif(T,min=hyper_a_1+1,max=hyper_b_1+1))

    r_t_modification <-  as.integer(runif(T,min=hyper_a_2+1,max=hyper_b_2+1))
  }else{
    # follow negative binomial distribution
    r_t_unmodification <-  rnbinom(T,mu=hyper_a_1, size=hyper_b_1)

    r_t_modification <-  rnbinom(T, mu=hyper_a_2, size=hyper_b_2)

  }



  # count the total number
  N_unmodification <- sum(r_t_unmodification)

  N_modification <- sum(r_t_modification)


  # calculate the alpha for both
  alpha_unmodification <- r_t_unmodification/N_unmodification

  alpha_modification <- r_t_modification/N_modification






  input_real_list <- list()
  for (i in 1: replicates){
    prop <- runif(1,0.9,1.1)

    input <- table(c(sample(1:T,round(prop*N_unmodification*1),
                            prob=alpha_unmodification,replace = TRUE),
                     sample(1: T,round(prop*N_modification*(0)),
                            prob=alpha_modification,replace = TRUE)))
    input_name <- as.numeric(names(input))
    input_value <- as.numeric(input)
    input <- rep(0,T)
    input[input_name] <- input_value
    input_real_list[[i]] <- input
  }

  IP_real_list <- list()
  for (i in 1: replicates){
    prop <- runif(1,0.9,1.1)

    IP <- table(c(sample(1:T,round(prop*N_modification*(0)),
                         prob=alpha_unmodification,replace = TRUE),
                  sample(1: T,round(prop*N_modification*(1)),
                         prob=alpha_modification,replace = TRUE)))
    IP_name <- as.numeric(names(IP))
    IP_value <- as.numeric(IP)
    IP <- rep(0,T)
    IP[IP_name] <- IP_value
    IP_real_list[[i]] <- IP

  }

  input_list <- list()
  for (i in 1: replicates){
    prop <- runif(1,0.9,1.1)

    input <- table(c(sample(1:T,round(prop*N_unmodification*pure_rate),
                            prob=alpha_unmodification,replace = TRUE),
                     sample(1: T,round(prop*N_modification*(1-pure_rate)),
                            prob=alpha_modification,replace = TRUE)))
    input_name <- as.numeric(names(input))
    input_value <- as.numeric(input)
    input <- rep(0,T)
    input[input_name] <- input_value
    input_list[[i]] <- input
  }

  IP_mix_list <- list()
  for (i in 1: replicates){
    prop <- runif(1,0.9,1.1)

    IP <- table(c(sample(1:T,round(prop*N_modification*mix_rate),
                         prob=alpha_unmodification,replace = TRUE),
                  sample(1: T,round(prop*N_modification*(1-mix_rate)),
                         prob=alpha_modification,replace = TRUE)))
    IP_name <- as.numeric(names(IP))
    IP_value <- as.numeric(IP)
    IP <- rep(0,T)
    IP[IP_name] <- IP_value
    IP_mix_list[[i]] <- IP

  }

  list(input_real_list,IP_real_list,input_list,IP_mix_list)
}
