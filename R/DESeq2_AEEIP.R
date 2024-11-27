#' DESeq2 with correct bias
#'
#' @param AEEIP_test_final_list list of AEEIP_test_final result of different replicates
#' @param input input read counts of different replicates
#' @param IP IP read counts of different replicates
#' @param times sample times
#' @param seed sample seed
#' @import DESeq2
#' @importFrom stats median
#' @importFrom stats rmultinom
#'
#' @return DESeq2 results
#' @export
#'
#'
DESeq2_AEEIP <- function(AEEIP_test_final_list, input, IP,times=100,seed=1){

  #check dim
  if (dim(input)[2]!=dim(IP)[2]){
    stop("the dimension of input and Ip should be same")
  }

  if (dim(input)[2]!=length(AEEIP_test_final_list)){
    stop("the dimension of input and Ip should be same as the length of the list")
  }

  counts <- cbind(as.data.frame(input),as.data.frame(IP))

  data_size <- as.matrix(counts)


  # first log
  log_data <- log(data_size)
  log_data [is.infinite(log_data)] <- NA
  log_mean <- rowMeans(log_data)
  log_s <- log_data-log_mean

  # then exp
  s_size <- exp(apply(log_s,2,function(x)median(x,na.rm=TRUE)))
  names(s_size) <- NULL


  rep_num <- dim(input)[2]

  s_size [(rep_num+1):(2*rep_num)] <- s_size [(rep_num+1):(2*rep_num)]*
    (unlist(lapply(AEEIP_test_final_list,function(x) x[["tau_b_bayes"]])))


  IP_new <- IP
  set.seed(seed)
  for (i in 1:rep_num){
    sample_prob=input[,i]/sum(input[,i])
    IP_new[,i] <- apply(IP_new[,i]-rmultinom(times,floor(sum(IP[,i])*AEEIP_test_final_list[[i]][["tau_a_bayes"]]),sample_prob),1,function(x) as.integer(median(x)))
  }
  IP_new[IP_new<0]<-0

  counts <- cbind(as.data.frame(input),as.data.frame(IP_new))

  design = data.frame(
    Trt = c(rep("input", dim(as.data.frame(input))[2]),
            rep("IP", dim(as.data.frame(IP_new))[2])))
  model = ~  Trt
  #
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = design,
                                design = model)



  sizeFactors(dds) = s_size

  dds <- DESeq(dds, test = "Wald")
  res_DESeq2_mix_corrected = DESeq2::results(dds)
  return(res_DESeq2_mix_corrected)




}
