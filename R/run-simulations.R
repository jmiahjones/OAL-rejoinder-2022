# inside function
# modified by Jeremiah Jones, 2021 to investigate positivity violations

### define some functions for generating data, ATE estimates, and the wAMD,
expit = plogis
ATE_est = function(fY,fw,fA){
  t_ATE = fY*fw
  tt_ATE = ( ( sum(t_ATE[fA==1]) / sum(fw[fA==1]) ) - ( sum(t_ATE[fA==0]) /  sum(fw[fA==0]) ) )
  return(tt_ATE) 
}
create_weights = function(fp,fA,fw){
  fw = (fp)^(-1)
  fw[fA==0] = (1 - fp[fA==0])^(-1)
  return(fw)
}
wAMD_function = function(DataM,varlist,trt.var,wgt,beta){
  trt = untrt = diff_vec = rep(NA,length(beta)) 
  names(trt) = names(untrt) = names(diff_vec) = varlist
  for(jj in 1:length(varlist)){ 
    this.var = paste("w",varlist[jj],sep="") 
    DataM[,this.var] = DataM[,varlist[jj]] * DataM[,wgt] 
    trt[jj] = sum( DataM[DataM[,trt.var]==1, this.var ]) / sum(DataM[DataM[,trt.var]==1, wgt]) 
    untrt[jj] = sum(DataM[DataM[,trt.var]==0, this.var]) / sum(DataM[DataM[,trt.var]==0, wgt]) 
    diff_vec[jj] = abs( trt[jj] - untrt[jj] ) 
  } 
  wdiff_vec = diff_vec * abs(beta) 
  wAMD = c( sum(wdiff_vec))
  ret = list( diff_vec = diff_vec, wdiff_vec = wdiff_vec, wAMD = wAMD )
  return(ret) 
}

ates <- sapply(1L:num_simulations, function(sim_idx){
  
  ### simulate data
  Mean_x = rep(mean_x,length(var.list))
  Data = as.data.frame(mvrnorm(n = n,mu=Mean_x,Sigma = Sigma_x,empirical = FALSE))
  names(Data) = var.list
  gA_x = rowSums(Data[,var.list]*matrix(alpha_v,nrow=n,ncol=length(var.list),byrow=TRUE))
  pA = expit( gA_x )
  summary(pA)
  Data$A = as.numeric( runif(n=length(pA)) < pA) # simulate A 
  gY_xA = rowSums(Data[,var.list]*matrix(beta_v,nrow=n,ncol=length(var.list),byrow=TRUE))   
  Data$Y = gY_xA + rnorm(n=n)
  Data$Y = Data$Y + Data$A*bA
  
  # Normlize coviarates to have mean 0 and standard deviation 1
  temp.mean = colMeans(Data[,var.list])
  Temp.mean = matrix(temp.mean,ncol=length(var.list),nrow=nrow(Data),byrow=TRUE)
  Data[,var.list] = Data[,var.list] - Temp.mean
  temp.sd = apply(Data[var.list],FUN=sd,MARGIN=2)
  Temp.sd = matrix(temp.sd,ncol=length(var.list),nrow=nrow(Data),byrow=TRUE)
  Data[var.list] = Data[,var.list] / Temp.sd
  rm(list=c("temp.mean","Temp.mean","temp.sd","Temp.sd"))
  
  # estimate outcome model
  y.form = formula(paste("Y~A+",paste(var.list,collapse="+")))
  lm.Y = lm(y.form,data=Data)
  betaXY = coef(lm.Y)[var.list] 
  
  ## Want to save ATE, wAMD and propensity score coefficients for each lambda value
  ATE = wAMD_vec = rep(NA, length(lambda_vec))
  names(ATE) = names(wAMD_vec) = names(lambda_vec)
  coeff_XA = as.data.frame(matrix(NA,nrow=length(var.list),ncol=length(lambda_vec)))
  names(coeff_XA) = names(lambda_vec)
  rownames(coeff_XA) = var.list
  
  ######################################################################################
  #####  Run outcome adaptive lasso for each lambda value 
  ######################################################################################
  # weight model with all possible covariates included, this is passed into lasso function
  w.full.form = formula(paste("A~",paste(var.list,collapse="+")))
  lil = names(lambda_vec)[1]
  for( lil in names(lambda_vec) ){
    il = lambda_vec[lil]
    ig = gamma_vals[lil]
    
    pen_w <- abs(betaXY)^(-ig) # gets internally rescaled to sum to p
    pen_w_c <- sum(pen_w)/p
    
    X <- as.matrix(Data[,1:p])
    logit_oal <- glmnet(X, Data$A, family=binomial(), intercept=FALSE,
                        penalty.factor=abs(betaXY)^(-ig))
    
    # save propensity score coefficients
    coeff_XA[var.list,lil] = as.numeric(coef(
      logit_oal, s=pen_w_c*(n^(il)), exact=TRUE,
      x=X, y=Data$A, family=binomial(), intercept=FALSE,
      penalty.factor=abs(betaXY)^(-ig)
    ))[-1]
    
    # generate propensity score
    Data[,paste("f.pA",lil,sep="")] = plogis(X%*%coeff_XA[var.list,lil])
    # create inverse probability of treatment weights
    Data[,paste("w",lil,sep="")] = create_weights(fp=Data[,paste("f.pA",lil,sep="")],fA=Data$A)
    
    # estimate weighted absolute mean different over all covaraites using this lambda to generate weights
    wAMD_vec[lil] = wAMD_function(DataM=Data,varlist=var.list,trt.var="A",
                                  wgt=paste("w",lil,sep=""),beta=betaXY)$wAMD
    
    # save ATE estimate for this lambda value
    ATE[lil] = ATE_est(fY=Data$Y,fw=Data[,paste("w",lil,sep="")],fA=Data$A)
  } # close loop through lambda values
  
  # print out wAMD for all the lambda values tried
  wAMD_vec
  # find the lambda value that creates the smallest wAMD
  tt = which.min(wAMD_vec)
  # print out ATE corresponding to smallest wAMD value
  ATE[tt]
})

mean(ates)
# print out the coefficients for the propensity score that corresponds with smalles wAMD value 
# coeff_XA[,tt]