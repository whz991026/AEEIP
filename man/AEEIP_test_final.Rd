% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AEEIP_test_final.R
\name{AEEIP_test_final}
\alias{AEEIP_test_final}
\title{This is the main function to extract the final result from AEEIP_test.}
\usage{
AEEIP_test_final(
  AEEIP_test_list,
  r_t_pure,
  r_t_mix,
  beta,
  alpha_a = NA,
  alpha_b = NA,
  tau_a,
  l,
  prop_vector = seq(0.01, 0.99, by = 0.01),
  max_iteration = 1e+06,
  stop_stepsize = 1e-07
)
}
\arguments{
\item{AEEIP_test_list}{list of AEEIP_test result of different replicates}

\item{r_t_pure}{Read counts of Input sample same as AEEIP_test_list}

\item{r_t_mix}{Read counts of IP sample same as AEEIP_test_list}

\item{beta}{the hyper parameter of Beta prior same as AEEIP_test_list}

\item{alpha_a}{the initial alpha_a same as AEEIP_test_list}

\item{alpha_b}{the initial alpha_b same as AEEIP_test_list}

\item{tau_a}{the initial tau_a same as AEEIP_test_list}

\item{l}{the effective length of site t same as AEEIP_test_list}

\item{prop_vector}{the vector of the prop of the beta hyper parameter same as AEEIP_test_list}

\item{max_iteration}{max iteration same as AEEIP_test_list}

\item{stop_stepsize}{the threshold for stop same as AEEIP_test_list}
}
\value{
final list, and two dataframe for visualization
}
\description{
This is the main function to extract the final result from AEEIP_test.
}
