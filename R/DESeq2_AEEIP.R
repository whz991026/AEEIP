#' DESeq2 with correct bias
#'
#' @param AEEIP_test_list list of AEEIP_test result of different replicates
#' @param input input read counts of different replicates
#' @param IP IP read counts of different replicates
#' @import DESeq2
#' @importFrom stats median
#'
#' @return DESeq2 results
#' @export
#'
#'
DESeq2_AEEIP <- function(AEEIP_test_list, input, IP){

  #check dim
  if (dim(input)[2]!=dim(IP)[2]){
    stop("the dimension of input and Ip should be same")
  }

  if (dim(input)[2]!=length(AEEIP_test_list)){
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
    (unlist(lapply(AEEIP_test_list,function(x) x[["tau_b_bayes"]])))



  for (i in 1:rep_num){
    IP[,i] <- round(IP[,i]-sum(IP[,i])*AEEIP_test_list[[i]][["alpha_a_bayes"]])
  }



  IP[IP<0] <-0

  counts <- cbind(as.data.frame(input),as.data.frame(IP))

  design = data.frame(
    Trt = c(rep("input", dim(as.data.frame(input))[2]),
            rep("IP", dim(as.data.frame(IP))[2])))
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
