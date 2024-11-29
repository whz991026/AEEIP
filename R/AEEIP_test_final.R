#' This is the main function to extract the final result from AEEIP_test.
#'
#' @param AEEIP_test_list list of AEEIP_test result of different replicates
#' @param r_t_pure Read counts of Input sample same as AEEIP_test_list
#' @param r_t_mix Read counts of IP sample same as AEEIP_test_list
#' @param beta the hyper parameter of Beta prior same as AEEIP_test_list
#' @param alpha_a the initial alpha_a same as AEEIP_test_list
#' @param alpha_b the initial alpha_b same as AEEIP_test_list
#' @param tau_a the initial tau_a same as AEEIP_test_list
#' @param l the effective length of site t same as AEEIP_test_list
#' @param prop_vector the vector of the prop of the beta hyper parameter same as AEEIP_test_list
#' @param max_iteration max iteration same as AEEIP_test_list
#' @param stop_stepsize the threshold for stop same as AEEIP_test_list
#'
#' @return final list, and two dataframe for visualization
#' @export
#'
#'
AEEIP_test_final <- function(AEEIP_test_list,r_t_pure,r_t_mix,beta,alpha_a=NA,alpha_b=NA,tau_a,l,
                             prop_vector=seq(0.01,0.99,by=0.01),
                             max_iteration =1e6,stop_stepsize=1e-7){

  # plot
  data.frame_1 <- as.data.frame(cbind(AEEIP_test_list$tau_a_no_bayes_vector_aeeip,
                                    AEEIP_test_list$tau_a_bayes_vector_aeeip))
  colnames(data.frame_1) <- c("tau_a_no_bayes","tau_a_bayes")

  # from the plot we can see that is near to 0.1


  data.frame_2 <- as.data.frame(cbind(AEEIP_test_list$tau_a_no_bayes_vector_aeeip,
         AEEIP_test_list$tau_a_no_bayes_vector_aeeip-AEEIP_test_list$lm_tau_vector))
  colnames(data.frame_2) <- c("tau_a_no_bayes","diff")

  #We selected the prop where the learned LM coefficient and the tau coefficient
  #changed. The corresponding index is determined as follows:

  # we select the smallest tau the derivative >0.01
  index <- min(intersect(which(diff(AEEIP_test_list$tau_a_no_bayes_vector_aeeip-
                            AEEIP_test_list$lm_tau_vector)/
                       diff(AEEIP_test_list$tau_a_no_bayes_vector_aeeip)>0.01),
                         which(abs(AEEIP_test_list$tau_a_no_bayes_vector_aeeip-AEEIP_test_list$lm_tau_vector)<=0.01)))+
    1 # +1 because the derivative is length-1

  if(length(index)==0){
    stop("no suit value for index")
    }


  res_final <- EM_function(r_t_pure=r_t_pure,r_t_mix=r_t_mix,beta=beta,alpha_a=alpha_a,
                           alpha_b=alpha_b,tau_a=tau_a,l=l,prop=prop_vector[index],
                           max_iteration =max_iteration,stop_stepsize=stop_stepsize)
  if(index<=2){
    print("please reduce the size parameter or the prop_vector of the AEEIP_test_list")
    print(paste0("the antibody bias is lower than ",res_final[["tau_a_bayes"]]))
  }else{
    print(paste0("the antibody bias is around ",res_final[["tau_a_bayes"]]))
  }
  return(list(res_final,data.frame_1,data.frame_2))
}
