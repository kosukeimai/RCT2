
###
### Sample size calculations for detecting a specific alternative 
### 
###



#' Sample size calculations for detecting a specific alternative 
#'
#' This function extracts the necessary outcome vector `Y`  for further analysis.
#' For the details of the method implemented by this function, see the
#' references.
#' 
#' 
#' @param J A vector of cluster labels for each observation of the dataset.
#' @param var A vector of the variable of interest, from which the binary treatment vector and outcome vector are to be extracted.
#' 
#' @return A list which contains the following item:
#' \item{Y}{ A list of the outcome variable for each observation. }
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
#' @name YSRE
#' 
#' 
#' 
#' @export YSRE
#' 
#' 
#' 

YSRE <- function(J, var){
  label <- unique(J)
  
  ## total number of clusters
  n <- length(label)
  Y <- vector("list", n)
  for(i in 1:n){
    Y[[i]] <- var[J == label[i]]
  }
  return(Y)
}


