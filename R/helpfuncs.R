#' Help functions for SMRmix
#' 
#' @import expm
#' @importFrom nlme lme
#' @import MASS
#' @importFrom lme4 glmer
#' @import stats
#' @noRd
# SMRmix for continuous outcomes
SMRmix.cont = function(H0.lm, H0.r, data, Ks, n.sam, n.study, y, X1,
                       n.kernel, times, var.type) {
  
  
  score = matrix(0, nrow = times, ncol = n.kernel)  #perturbed statistics
  true = rep(0, n.kernel)                           #test statistic
  
  # whether equal individual-level variance
  if (var.type == "different") {
    lme = lme( H0.lm ,random = H0.r,  data, weights = varIdent(form = H0.r))
    delta2 = (1/unique(attributes(lme$modelStruct$varStruct)$weights))^2 * as.numeric(nlme::VarCorr(lme)[2,1])
    sigma2 = as.numeric(nlme::VarCorr(lme)[1,1])
  }
  
  if (var.type == "same") {
    lme = lme(H0.lm, random = H0.r, data)
    delta2 = as.numeric(nlme::VarCorr(lme)[2,1])
    delta2 = rep(delta2, n.study)
    sigma2 = as.numeric(nlme::VarCorr(lme)[1,1])
  }
  
  fity = predict(lme, re.form =~ 0)
  
  for (k in 1:n.study) {
    
    index = sum(n.sam[1:k])
    Kk = lapply(Ks, FUN = function(K) {K[[k]]})
    Sk =  lapply(Kk, FUN = function(K) {K[lower.tri(K)]})
    fityk = fity[(index - n.sam[k]+1):(index)]
    yk = y[(index - n.sam[k]+1):(index)]
    
    invV = CalV(n.sam[k], sigma2 = sigma2, delta2 = delta2[k])   # inverse var-cov
    weightk = lapply(Sk, FUN = function(Sk) {t(Sk)%*%invV} ) 
    mvcov = matrix(sigma2, n.sam[k], n.sam[k])
    diag(mvcov) = sigma2+delta2[k]
    X1k = cbind(1, X1[(index - n.sam[k]+1):(index),])
    Psqrt = sqrtm(mvcov)  - X1k%*% solve(t(X1k) %*% solve(mvcov) %*% X1k) %*% t(X1k) %*% solve(sqrtm(mvcov))
    
    for (m in 1:times) {
      perturb = mvrnorm(n = 1,mu = rep(0, n.sam[k]), Sigma =  diag(n.sam[k]))
      perZk = tcrossprod(Psqrt %*% perturb, Psqrt %*%perturb)
      perZk = perZk[lower.tri(perZk)]
      score[m, ] = score[m, ] + unlist( lapply(weightk, function(x) sum(x*(perZk-delta2[k]))) )
    }
    z = outer(yk-fityk,yk-fityk)
    
    true = true + unlist(lapply(weightk, function(x) sum(x*(z[lower.tri(z)]-delta2[k]))))
  }
  p.ind = rep(NA, n.kernel)
  
  for (i in 1:n.kernel) {
    p.ind[i] = mean(abs(true[i]) < abs(score[,i]))
  }
  return(p.ind)
}


# SMRmix for binary outcomes
SMRmix.bi = function(H0.lm, H0.r, data, Ks, n.sam, n.study, y, X1,
                     n.kernel, times) {
  
  score = matrix(0, nrow = times, ncol = n.kernel)  #perturbed statistics
  true = rep(0, n.kernel)                           #test statistic
  
  lme = glmer(formula.H0, family = binomial)
  sigma2 = as.numeric(lme4::VarCorr(lme))
  fity = predict(lme, re.form =~ 0, type = "response")
  
  for (k in 1:n.study) {
    
    index = sum(n.sam[1:k])
    Kk = lapply(Ks, FUN = function(K) { K[[k]] })
    Sk =  lapply(Kk, FUN = function(Kk) { Kk[lower.tri(Kk)] })
    
    fityk = fity[(index - n.sam[k]+1):(index)]
    invV = CalV_binary(n.sam[k], sigma2 = sigma2, mu0 = fityk)
    yk = y[(index - n.sam[k]+1):(index)]
    
    weightk = lapply(Sk, FUN = function(Sk) {t(Sk)%*%invV} ) 
    
    W = matrix(1, n.sam[k], n.sam[k])*sigma2 + diag(1/(fityk * (1-fityk)))
    remove.col = which(apply(X1[(index - n.sam[k]+1):(index),], 2, function(x) length(unique(x))) == 1)
    if (length(remove.col) > 0) {
      X1k = X1[(index - n.sam[k]+1):(index), - remove.col, drop =FALSE]
    } else {
      X1k = cbind(rep(1, n.sam[k]), X1[(index - n.sam[k]+1):(index),])
    }
    
    Psqrt = sqrtm(W)- X1k%*% chol2inv(chol(t(X1k) %*% solve(W) %*% X1k)) %*% t(X1k) %*% solve(sqrtm(W))
    
    mvcov = outer( fityk*(1-fityk), fityk*(1-fityk))*sigma2
    diag(mvcov) = fityk*(1-fityk)
    
    eigen = eigen(mvcov)
    eigen$values[eigen$values<0] = 0
    mvcov = eigen$vectors %*% diag(eigen$values) %*% t(eigen$vectors)
    
    D = diag(fityk*(1-fityk))
    for (m in 1:times) {
      perturb = mvrnorm(n = 1,mu = rep(0, n.sam[k]), Sigma =  diag(n.sam[k]))
      perZk = tcrossprod(D%*% Psqrt %*% perturb, D%*% Psqrt %*%perturb)
      perZk = perZk[lower.tri(perZk)]
      score[m, ]  = score[m, ] + unlist(lapply(weightk, function(x) sum(x*(perZk-mvcov[lower.tri(mvcov)]))))
      
    }
    z = outer(yk-fityk,yk-fityk)
    
    true = true + unlist( lapply(weightk, function(x) sum(x*(z[lower.tri(z)]-mvcov[lower.tri(mvcov)]))) )
    
  }
  p.ind = rep(NA, n.kernel)
  
  for (i in 1:n.kernel) {
    p.ind[i] = mean(abs(true[i]) < abs(score[,i]))
  }
  return(p.ind)
}


# calculate variance for continuous outcomes
CalV = function(n, sigma2, delta2) {
  V = matrix(2*sigma2^2, ncol = n*(n-1)/2, nrow = n*(n-1)/2)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      m = (2*n-i)*(i-1)/2+j-i
      if (j != n) {
        V[m, (m+1):(m+n-j)] = 2*sigma2^2+sigma2*delta2
        V[m, ((2*n-j)*(j-1)/2+1):((2*n-j)*(j-1)/2+n-j)] = 2*sigma2^2+sigma2*delta2
      }
      
      if (j-i != 1 && i != n-1) {
        for (k in 1:(j-i-1)) {
          V[m, ((2*n-k-i)*(k+i-1)/2+j-k-i)] = 2*sigma2^2+sigma2*delta2
        }
      }
      
    }
  }
  
  diag(V) = 2*sigma2^2+2*sigma2*delta2+delta2^2
  V[lower.tri(V)] = t(V)[lower.tri(V)]
  chol2inv(chol(V))
}


# calculate variance for binary outcomes
CalV_binary = function(n, sigma2, mu0) {
  gprime = mu0*(1-mu0)
  tmp = outer(gprime, gprime) 
  tmp = tmp[lower.tri(tmp)]
  V = outer(tmp,tmp)*2*sigma2^2
  tmp2 = outer((1-2*mu0)*gprime,(1-2*mu0)*gprime) + outer(gprime,gprime^2) + outer(gprime^2,gprime)
  tmp2 = tmp2[lower.tri(tmp2)]
  diag(V) = sigma2*tmp2+2*sigma2^2*tmp^2+tmp
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      m = (2*n-i)*(i-1)/2+j-i
      if (j != n) {
        V[m, (m+1):(m+n-j)] = 2*sigma2^2*gprime[i]^2*gprime[j]*gprime[(j+1):n]+ sigma2*gprime[i]*gprime[j]*gprime[(j+1):n]
        V[m, ((2*n-j)*(j-1)/2+1):((2*n-j)*(j-1)/2+n-j)] = 2*sigma2^2*gprime[j]^2*gprime[i]*gprime[(j+1):n]+ sigma2*gprime[j]*gprime[i]*gprime[(j+1):n]
      }
      if (j-i != 1 && i != n-1) {
        for (k in 1:(j-i-1)) {
          V[m, ((2*n-k-i)*(k+i-1)/2+j-k-i)] = 2*sigma2^2*gprime[j]^2*gprime[i]*gprime[i+k]+ sigma2*gprime[j]*gprime[i]*gprime[i+k]
        }
      }
      
    }
  }
  V[lower.tri(V)] = t(V)[lower.tri(V)]
  eigen = eigen(V)
  eigen$values[eigen$values<0] = -eigen$values[eigen$values<0]
  V = eigen$vectors %*% diag(eigen$values) %*% t(eigen$vectors)
  chol2inv(chol(V))
}


