
###
### Sample size calculations for detecting a specific alternative 
### 
###



#' Sample size calculations for detecting a specific alternative 
#'
#' This function extracts the binary treatment assignment vector `Z`.
#' For the details of the method implemented by this function, see the
#' references.
#' 
#' 
#' @param assign The assignment vector for each observation of the data.
#' @param J A vector of cluster labels for each observation of the dataset.
#' 
#' @return A list which contains the following item:
#' \item{Z}{ A  list of the binary treatment assignment variable for each observation. }
#' 
#' 
#' @author Kosuke Imai, Department of Statistics, Harvard University
#' \email{imai@harvard.edu}, \url{https://imai.fas.harvard.edu/};
#' Zhichao Jiang, School of Public Health and Health Sciences, University of Massachusetts Amherst
#' \email{zhichaojiang@umass.edu};
#' Karissa Huang, Department of Statistics, Harvard College
#' \email{krhuang@college.harvard.edu}
#' @references Zhichao Jiang, Kosuke Imai (2020).
#' \dQuote{Statistical Inference and Power Analysis for Direct and Spillover Effects in Two-Stage Randomized Experiments}, \emph{Technical Report}.
#' @keywords two-stage randomized experiments
#' @name ZSRE
#' 
#' 
#' 
#' @export ZSRE
#' 
#' 
#' 

ZSRE <- function(assign, J){
  label <- unique(J)
  
  ## total number of clusters
  n <- length(label)
  Z <- vector("list", n)
  for(i in 1:n){
    Z[[i]] <- assign[J == label[i]]
  }
  return(Z)
}


