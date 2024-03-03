#' Mixed Effect Similarity Matrix Regression Model (SMRmix) 
#' 
#' Integrates heterogeneous microbiome datasets and tests the association between 
#' microbial composition and an outcome of interest.
#'
#'
#' @param otu.tab Numeric matrix of template microbiome count dataset. Rows are
#' samples, columns are taxa.
#' @param n.break.ties Number of replicates to break ties when ranking relative
#' abundances. Defaults to \code{100}.
#' @param mode A character indicating the modeling approach for relative abundances.
#' If \code{'parametric'}, a parametric model involving fitting a generalized
#' gamma distribution is used. If \code{'nonparametric'}, the nonparametric
#' approach involving quantile matching is applied. Note that a parametric
#' model is required if library sizes or characteristics of taxa will be modified.
#' Defaults to \code{'nonparametric'}.
#' @return Returns a list that has components:
#' @param formula A two-sided linear formula describing both the fixed-effects 
#' and random-effects part of the model.
#' @param data An optional data frame containing the variables named in \code{formula}.
#' @param Kernels A list that either consists of multiple lists of kernel matrices, or 
#' a single list of kernel matrices. Denote Kernels\eqn{=\{K_1,\cdots,K_P\}}, where 
#' \eqn{P} is the number of kernel metrics under consideration. Each \eqn{K_p} (\eqn{p\in \{1,\cdots,P\}}) 
#' is itself a list of \eqn{K} kernel matrices, where \eqn{K} is the number of 
#' heterogeneous studies. Each study \eqn{k} (\eqn{k\in \{1,\cdots,K\}}) contributes a kernel
#' matrix of dimension \eqn{n_k} by \eqn{n_k} to each \eqn{K_p}.
#' @param times Number of perturbation replicates. Defaults to \code{2000}.
#' @param var.type A character indicating the model assumption for the variance of
#' study-specific random effect for continuous outcomes. The default option (\code{'different'}) 
#' assumes the variance \eqn{\delta_k^2} vary across different studies. The option 
#' \code{var.type='same'} assumes a constant variance across all studies.
#' @param outcome A character indicating whether the outcome of interest is \code{continuous}
#' or \code{binary}. Default to \code{continuous}.
#'
#' @return Returns a list that has components:
#' \item{p_values}{Individual p-values for each candidate kernel metric.}
#' \item{omnibus_p}{Omnibus p-value considering multiple kernel metrics.}
#'
#' @import harmonicmeanp
#' @import stats
#' 
#' @export
SMRmix = function(formula, data = NULL, Kernels, times = 2000, var.type = "different", 
                  outcome = "continuous"){
  
  H0.lmer <- formula
  temp <- as.character(H0.lmer)
  H0.r <- as.formula( paste("~",sub("\\).*", "", sub(".*\\(", "", temp[3])))  )
  temp[3] <- gsub('\\+\\s*\\(.*?\\)', '',  temp[3], perl=TRUE)
  H0.lm <- as.formula(paste(temp[2], temp[1], temp[3]))
  
  X1 = NULL
  tmp = unlist(strsplit( gsub("\\s", "",temp[3]), "\\+")[[1]])[-1]
  if (length(tmp) !=0 ) {
    for (i in 1:length(tmp)) {
      X1 = cbind(X1, eval(as.name(tmp[i])))
    }
  }
  
  # how many kernels(whether omnibus)
  if (typeof(Kernels[[1]]) == "list") { 
    n.kernel = length(Kernels) 
    n.study = length(Kernels[[1]])
    n.sam = unlist(lapply(Kernels[[1]], function(x) dim(x)[1]))
  }
  if (typeof(Kernels[[1]]) == "double") {
    n.kernel = 1
    n.study = length(Kernels)
    n.sam = unlist(lapply(Kernels, function(x) dim(x)[1]))
    Kernels = list(Kernels)
  }
  
  
  if (outcome == "continuous") {  
    p.ind = SMRmix.cont(H0.lm, H0.r, data, Kernels, n.sam, n.study, y, X1,
                        n.kernel, times, var.type)  
  }
  if (outcome == "binary") { 
    p.ind = SMRmix.bi(H0.lm, H0.r, data, Kernels, n.sam, n.study, y, X1,
                      n.kernel, times) 
  }
  
  if (length(p.ind) == 1) {
    names(p.ind) = names(Kernels)
    return(p.ind) 
  }
  
  if (length(p.ind) > 1) {
    p.ind.pos = ifelse(p.ind == 0, 10^-8, p.ind)
    p.hmp = p.hmp(p = p.ind.hmp, w = NULL, L = length(p.ind.hmp), 
                  w.sum.tolerance = 1e-6, multilevel = F)[1]
    
    p.omnibus = c(p.hmp)
    names(p.omnibus) = c("HMP")
    names(p.ind) = names(Kernels)
    
    return(list(p_values = p.ind, omnibus_p = p.omnibus)) 
  }
  
}


