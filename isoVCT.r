######################
###Data Preparation
######################


###Generation of Label matrix (X matrix in the article)
##labId = The array indices of the case group
##label = The label of samples
##p = The number of isoforms of one gene
labGen <- function(labId, label, p){
  label <- matrix(0, nrow = length(label)*p, ncol = p)
  for(j in 1: length(labId)){
    label[c(c((labId[j]-1)*p+1):c((labId[j]-1)*p+p)), ] <- diag(p)
  }
  return(label)
}






#######################
#####Calculation
#######################


###Calculation of U statistics
##dis = c('Poi', 'NB'), The distribution assumption of data
##label = The label of samples
##p = The number of isoforms of one gene
##fit = Result of fitting the NB model
##y = The vector of isoforms expression levels
##phi = The vector of diagonal matrix
UEst <- function(dis, label, p, fit, y){
  labId <- which(label == 1)
  labelMat <- labGen(labId, label,p)
  if(dis == 'Poi'){
    xMat <- matrix(y-predict(fit, type = 'response'), 
                   ncol = p*length(label)) %*% labelMat
  }else{
    xMat <- matrix(as.vector(y-predict(fit, type = 'response'))*phi, 
                   ncol=p*length(label)) %*%labelMat
  }
  return(xMat%*%t(xMat))
}


###Permutation the empirical of U statistics
UPer <- function(dis, label, p, fit, y){
  labId <- which(sample(label) == 1)
  labelMat <- labGen(labId, label, p)
  if(dis == 'Poi'){
    xMat <- matrix(y-predict(fit, type = 'response'), 
                   ncol = p*length(label)) %*% labelMat
  }else{
    phi <- rep(1/(1+fit@beta/fit@theta), p)
    xMat <- matrix(as.vector(y-predict(fit, type = 'response'))*phi, 
                   ncol=p*length(label)) %*%labelMat
  }
  return(xMat%*%t(xMat))
}






#######################
###Data Analysis
#######################


###variance component testing for isoforms
##data = isoform expression data
##nper = The number of permutation
isoVCT <- function(data, dis, nper){
  require(lme4)
  p <- dim(data)[2] - 1
  if (p <= 1){
    stop('Number of isoform is one; use glm().')
  }
  if (dis != 'Poi' & dis != 'NB'){
    stop('dis is "Poi" or "NB".')
  }
  group <- rep(c(1: p), nrow(data))
  y <- matrix(t(as.matrix(data[, -1])), ncol = 1)
  if(dis == 'Poi'){
    ##Fit the mixed model under the H0 
    fitPoi <- glmer(y ~ (1|group), family = 'poisson')
    ##Calculation of P-value of statistics U
    U <- UEst(dis = dis, label = data[, 1], p = p,  
              fit = fitPoi, y = y)
    UDis <- replicate(nper, UPer(dis = dis, label = data[, 1], p = p,  
                                 fit = fitPoi, y = y))
    pval <- sum(as.numeric(U) <= UDis)/nper
  }else{
    ##Fit the mixed model under the H0
    fitNB <- try(glmer.nb(y~(1|group)), silent=T)  
    if(inherits(fitNB, 'try-error')){
      stop('The data is unfit to the NB.')
    }else{
      ##Calculation of P-value of statistics U
      U <- UEst(dis = dis, label = data[, 1], p = p,  
                fit = fitNB, y = y)
      UDis <- replicate(nper, UPer(dis = dis, label = data[, 1], p = p,  
                                   fit = fitNB, y = y))
      pval <- sum(as.numeric(U) <= UDis)/nper
    }
  }
  return(list(U = as.numeric(U), pval = pval))
}
