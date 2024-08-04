
#library(survival)
#######################
### CA randomization ###
########################


urn_design<-function(n_within_strata,omega,gamma,beta){ # Wei's urn design (Wei, 1978, JASA, Vol.73, No. 363)
  I<-numeric(n_within_strata)
  n_balls_0<-omega
  n_balls_1<-omega
  for(i in 1:n_within_strata){
    p<-n_balls_1/(n_balls_1+n_balls_0)
    I[i]<-rbinom(1,1,p)
    n_balls_0<-n_balls_0+(1-I[i])*gamma+I[i]*beta
    n_balls_1<-n_balls_1+(1-I[i])*beta+I[i]*gamma
  }
  return(I)
}

strata_urn_design<-function(z){
  omega<-1
  gamma<-0
  beta<-1
  n_strata<-dim(z)[2]
  n<-dim(z)[1]
  n_each_strata<-apply(z,2,sum)
  I<-numeric(n)
  for(i in 1:n_strata){
    if(n_each_strata[i]==0) {next}
    else if (n_each_strata[i]==1) {I[z[,i]==1]<-rbinom(1,1,0.5)}
    else{
      I[z[,i]==1]<-urn_design(n_each_strata[i],omega, gamma, beta)
    }
  }
  return(I)
}

discretize_z<-function(n,z_mean,z_sd,n_category){ # only normal
  Z<-rnorm(n,z_mean,z_sd)
  p_cut<-seq(0,1, length.out = (n_category+1))
  # p_cut = c(0, 0.5, 0.8, 1) #for case3_K=3
  # p_cut = c(0, 0.4, 0.7, 0.9, 1) #for case3_K=4
  # p_cut = c(0, 0.26, 0.74, 1) #for caseR01_3
  # p_cut = c(0, 0.25, 0.6, 0.8, 1) #for caseL_3
  q_cut<-sapply(p_cut,function(x) qnorm(x,z_mean,z_sd))
  Z_dis_scale<-as.factor(sapply(Z,function(x) max(which(q_cut<x))))
  Z_model_matrix<-model.matrix(~0+Z_dis_scale)
  return(list(Z=Z,Z_dis=Z_model_matrix))
}

# simple random-----------
# SR <- function(n, p_trt){
#   I<-rbinom(n,1,p_trt)
#   return(I)
# }

SR <- function(n, p_trt){
  if((n %% 2) == 0) I <- sample(c(rep(0,n/2),rep(1,n/2)))
  if((n %% 2) != 0) I <- sample( c(rep(0,n/2),rep(1,n/2),rbinom(1,1,p_trt)))
  return(I)
}

# covariate adaptive----- #BCD_p: biased proabaility
BC<-function(n_within_strata,BCD_p=2/3){ # Efron's biased coin randomization (Efron, 1971)
  stopifnot(BCD_p>1/2)
  I<-numeric(n_within_strata)
  I[1]<-rbinom(1,1,1/2)
  D<-2*I[1]-1
  for(i in 2:n_within_strata){
    if(D>0) I[i]<-rbinom(1,1,(1-BCD_p))
    else if(D<0) I[i]<-rbinom(1,1,BCD_p)
    else if(D==0) I[i]<-rbinom(1,1,1/2)
    D<-D+2*I[i]-1
  }
  return(I)
}


CABC <- function(z, BCD_p=2/3){ 
  n_strata<-dim(z)[2]
  n<-dim(z)[1]
  n_each_strata<-apply(z,2,sum)
  I<-numeric(n)
  for(s in 1:n_strata){
    if(n_each_strata[s]==0) {next}
    else if (n_each_strata[s]==1) {I[z[,s]==1]<-rbinom(1,1,0.5)}
    else{
      I[z[,s]==1]<-BC(n_each_strata[s],BCD_p)
    }
  }
  return(I)
}

# permuted block -----------  
permuted_block <- function(z, blocksize, p_trt){
  n_trt_blk<-blocksize*p_trt
  n_ctl_blk<-blocksize*(1-p_trt)
  if(!(floor(n_ctl_blk)==n_ctl_blk & floor(n_trt_blk)==n_trt_blk)){stop("number of treatment and control within block are not integers!")}
  if(blocksize)
    n_strata<-dim(z)[2]
  n<-dim(z)[1]
  n_each_strata<-apply(z,2,sum)
  I<-numeric(n)
  for(i in 1:n_strata){
    I_tmp<-0
    n_block_tmp<-ceiling(n_each_strata[i]/blocksize)
    for(b in 1:n_block_tmp){
      block_permuted_tmp<-sample(c(rep(0,n_ctl_blk),rep(1,n_trt_blk)))
      I_tmp<-c(I_tmp,block_permuted_tmp)
    }
    I_tmp<-I_tmp[-1]
    I[z[,i]==1]<-I_tmp[1:n_each_strata[i]]
  }
  return(I)
}


# treatment_assignment -----------  
treatment_assignment <- function(n,strata_z,minimization_z,randomization,p_trt){
  if(randomization=="SR"){
    I<-SR(n,p_trt)
  }else if(randomization=="CABC"){
    if(p_trt!=1/2){stop("For now, the proportion of treatment under CABC has to be 1/2.")}
    I<-CABC(strata_z)
  }else if(randomization=="permuted_block"){
    if(p_trt!=1/2){stop("For now, the proportion of treatment under permuted_block has to be 1/2.")}
    I<-permuted_block(strata_z,4,p_trt)
    # I<-permuted_block(strata_z,10,p_trt)
  }else if(randomization=="minimization"){
    if(p_trt!=1/2){stop("For now, the proportion of treatment under CABC has to be 1/2.")}
    I<-minimization(minimization_z)
  }else if(randomization=="urn"){
    if(p_trt!=1/2){stop("For now, the proportion of treatment under CABC has to be 1/2.")}
    I<-strata_urn_design(strata_z)
  }
  return(I)
}


#########################
######### Tests #########
#########################

data.sort<-function(data.simu){
  n<-dim(data.simu)[1]
  data.simu.order<- data.simu[order(data.simu$t),]
  data.rev<- data.frame(data.simu.order, Y=n:1, Y1=cumsum(data.simu.order$I1[n:1])[n:1], Y0=cumsum(data.simu.order$I0[n:1])[n:1])
  return(data.rev)
}

wlogrank<-function(data.simu){ 
  n<-dim(data.simu)[1]
  data.rev<-data.sort(data.simu)
  T.seq<-data.rev$t
  T.rep<-c(0,T.seq)[1:n]
  T.diff<-T.seq-T.rep
  same.ind<-which(T.diff==0)
  n_col<-dim(data.rev)[2]
  for(ind in same.ind){
    data.rev[ind,(n_col-2):n_col]<-data.rev[ind-1,(n_col-2):n_col]
  }
  S.alt<-sum(data.rev$delta*(data.rev$I1*data.rev$Y0-data.rev$I0*data.rev$Y1)/data.rev$Y)
  var.alt<-sum(data.rev$delta*data.rev$Y0*data.rev$Y1/data.rev$Y^2)
  res<-S.alt/sqrt(var.alt)
  return(list(S.alt=S.alt, res=res,var.alt=var.alt))
}

ind_to_factor<-function(data.simu){
  z<-data.simu[,grepl("strata",names(data.simu))] 
  n_strata<-dim(z)[2]
  n<-dim(z)[1]
  fz<-numeric(n)
  for(i in 1:n_strata){
    fz<-fz+z[,i]*i
  }
  return(as.factor(fz))
}

# model based Wald test
regular_wald<-function(data.simu){
  x<-data.simu[,grepl("model",names(data.simu))] #only those has model
  data_sub<-data.frame(t=data.simu$t,delta=data.simu$delta,I1=data.simu$I1,x)
  fit<-coxph(Surv(t,delta)~.,data=data_sub)
  wald_res<-fit$coefficients[1]/sqrt(fit$var[1,1])
  return(wald_res)
}

# randomization test to construct reference distribution (based on regular T_LR)
wlogrank_rand_distribution<-function(data.simu,randomization,p_trt,alpha=0.05,boot_n=200){
  n<-dim(data.simu)[1]
  T_logrank<-wlogrank(data.simu)$S.alt
  strata_z<-data.simu[,grepl("strata",names(data.simu))] 
  minimization_z<-data.simu[,grepl("minimization",names(data.simu))] 
  LR.rand<-numeric(boot_n)
  data.rand<-data.simu
  for(boot in 1:boot_n){
    I.rand<-treatment_assignment(n,strata_z,minimization_z,randomization,p_trt)
    data.rand$I1<-I.rand
    data.rand$I0<-1-I.rand
    LR.rand[boot]<-wlogrank(data.rand)$S.alt
  }
  T_logrank.rand.ordered<-sort(LR.rand)
  c1<-T_logrank.rand.ordered[ceiling(boot_n*alpha/2)]
  c2<-T_logrank.rand.ordered[ceiling(boot_n*(1-alpha/2))+1]
  T_rand_rej<-1*(T_logrank<c1 | T_logrank>c2)
  return(T_rand_rej)
}


# proposed T_RL in Ye and Shao, 2020, jrssb.
wlogrank_robust<-function(data.simu,randomization,p_trt){ # try new version 19 May
  n<-dim(data.simu)[1]
  S.alt<-wlogrank(data.simu)$S.alt
  data.rev<-data.sort(data.simu)
  strata_z<-data.rev[,grepl("strata",names(data.rev))]
  mean_at_risk<-data.rev$delta/data.rev$Y
  cumsum_at_risk<-cumsum(mean_at_risk)
  mean_at_risk_2<-data.rev$delta*data.rev$Y1/(data.rev$Y)^2
  cumsum_at_risk_2<-cumsum(mean_at_risk_2)
  mu_t<-data.rev$Y1/data.rev$Y
  O_i<-data.rev$delta*(data.rev$I1-data.rev$Y1/data.rev$Y)-data.rev$I1*cumsum_at_risk+cumsum_at_risk_2
  O_i1<-data.rev$I1*O_i
  O_i0<- -data.rev$I0*O_i
  strata_z<-data.rev[,grepl("strata",names(data.rev))] 
  n_strata<-dim(strata_z)[2]
  z<-strata_z
  cond_var_1<-numeric(n_strata)
  cond_var_0<-numeric(n_strata)
  prob_z<-numeric(n_strata)
  cond_exp_1<-numeric(n_strata)
  cond_exp_0<-numeric(n_strata)
  for(i in 1:n_strata){
    cond_var_1[i]<-var(O_i1[z[,i]==1 & data.rev$I1==1])
    cond_var_0[i]<-var(O_i0[z[,i]==1 & data.rev$I0==1])
    prob_z[i]<-mean(z[,i])
    cond_exp_1[i]<-mean(O_i1[z[,i]==1 & data.rev$I1==1])
    cond_exp_0[i]<-mean(O_i0[z[,i]==1 & data.rev$I0==1])
  }
  na.ind<-is.na(cond_var_1+cond_var_0)
  if(length(which(na.ind))!=0) warning("conditional variance in log-rank equals 0")
  cond_var_1<-cond_var_1[!na.ind]
  cond_var_0<-cond_var_0[!na.ind]
  prob_z<-prob_z[!na.ind]
  cond_exp_1<-cond_exp_1[!na.ind]
  cond_exp_0<-cond_exp_0[!na.ind]
  sai_1<-(p_trt*(cond_var_1%*%prob_z)+(1-p_trt)*(cond_var_0%*%prob_z))
  sai_2<-0 # for cabc and permuted_block
  if(randomization=="SR"){
    if(p_trt!=1/2)  stop("For now, SR is only designed for pi=1/2")
    sai_2<-1/4*sum((cond_exp_0+cond_exp_1)^2 * prob_z)
    
  }
  if(randomization =="urn"){
    if(p_trt!=1/2)  stop("For now, urn design is only designed for pi=1/2")
    sai_2<-1/12*sum((cond_exp_0+cond_exp_1)^2 * prob_z)
  }
  var_cal<-n*(sai_1+sai_2)
  T_cal<-S.alt/sqrt(var_cal)
  return(list(T_cal,var_cal))
}


# score test T_S in Ye and Shao, 2020, jrssb.
score<-function(data.simu){
  data.rev<-data.sort(data.simu)
  x<-data.rev[,grepl("model",names(data.rev))] #only those has model
  n<-dim(data.rev)[1]
  data_sub<-data.frame(t=data.rev$t,delta=data.rev$delta, x)
  fit<-coxph(Surv(t,delta)~.,data=data_sub)
  S_0_seq<-exp(fit$coefficients %*% t(x))  ##exp(beta*Wi)
  S_1_seq<-exp(fit$coefficients %*% t(x))*data.rev$I1
  S_0<-cumsum(S_0_seq[n:1])[n:1]/n
  S_1<-cumsum(S_1_seq[n:1])[n:1]/n
  U<-sum(data.rev$delta*(data.rev$I1-S_1/S_0))
  #mu_t<-data.rev$Y1/data.rev$Y # two choices, both valid
  mu_t<-S_1/S_0
  mean_at_risk<-data.rev$delta/(n*S_0)
  cumsum_at_risk<-cumsum(mean_at_risk)
  mean_at_risk_2<-data.rev$delta/(n*S_0)*mu_t
  cumsum_at_risk_2<-cumsum(mean_at_risk_2)
  O_i<-data.rev$delta*(data.rev$I1-mu_t)-S_1_seq*cumsum_at_risk+S_0_seq*cumsum_at_risk_2
  var_cal<-mean(O_i^2)*n
  T_WR<-U/sqrt(var_cal)
  return(list(T_WR, U, var_cal))
}

# proposed T_RS in Ye and Shao, 2020, jrssb.
score_robust<-function(data.simu,randomization,p_trt){
  data.rev <- data.sort(data.simu)
  x <- data.rev[,grepl("model",names(data.rev))] #only those has model
  n<-dim(data.rev)[1]
  data_sub<-data.frame(t=data.rev$t,delta=data.rev$delta,x)
  fit<-coxph(Surv(t,delta)~.,data=data_sub)
  S_0_seq<-exp(fit$coefficients %*% t(x))
  S_1_seq<-exp(fit$coefficients %*% t(x))*data.rev$I1
  S_0<-cumsum(S_0_seq[n:1])[n:1]/n
  S_1<-cumsum(S_1_seq[n:1])[n:1]/n
  U<-sum(data.rev$delta*(data.rev$I1-S_1/S_0))
  #mu_t<-data.rev$Y1/data.rev$Y # two choices, both valid
  mu_t<-S_1/S_0
  mean_at_risk<-data.rev$delta/(n*S_0)
  cumsum_at_risk<-cumsum(mean_at_risk)
  mean_at_risk_2<-data.rev$delta/(n*S_0)*mu_t
  cumsum_at_risk_2<-cumsum(mean_at_risk_2)
  O_i<-data.rev$delta*(data.rev$I1-mu_t)-S_1_seq*cumsum_at_risk+S_0_seq*cumsum_at_risk_2
  
  O_i1<-data.rev$I1*O_i
  O_i0<- -data.rev$I0*O_i
  
  strata_z<-data.rev[,grepl("strata",names(data.rev))] 
  n_strata<-dim(strata_z)[2]
  z<-strata_z
  cond_var_1<-numeric(n_strata)
  cond_var_0<-numeric(n_strata)
  prob_z<-numeric(n_strata)
  cond_exp_1<-numeric(n_strata)
  cond_exp_0<-numeric(n_strata)
  for(i in 1:n_strata){
    cond_var_1[i]<-var(O_i1[z[,i]==1 & data.rev$I1==1])
    cond_var_0[i]<-var(O_i0[z[,i]==1 & data.rev$I0==1])
    prob_z[i]<-mean(z[,i])
    cond_exp_1[i]<-mean(O_i1[z[,i]==1 & data.rev$I1==1])
    cond_exp_0[i]<-mean(O_i0[z[,i]==1 & data.rev$I0==1])
  }
  na.ind<-is.na(cond_var_1+cond_var_0)
  if(length(which(na.ind))!=0) warning("conditional variance in log-rank equals 0")
  cond_var_1<-cond_var_1[!na.ind]
  cond_var_0<-cond_var_0[!na.ind]
  prob_z<-prob_z[!na.ind]
  cond_exp_1<-cond_exp_1[!na.ind]
  cond_exp_0<-cond_exp_0[!na.ind]
  sai_1<-(p_trt*(cond_var_1%*%prob_z)+(1-p_trt)*(cond_var_0%*%prob_z))
  sai_2<-0 # for cabc and permuted_block
  if(randomization=="SR"){
    if(p_trt!=1/2)  stop("For now, SR is only designed for pi=1/2")
    sai_2<-1/4*sum((cond_exp_0+cond_exp_1)^2 * prob_z)
  }
  if(randomization =="urn"){
    if(p_trt!=1/2)  stop("For now, urn design is only designed for pi=1/2")
    sai_2<-1/12*sum((cond_exp_0+cond_exp_1)^2 * prob_z)
  }
  var_cal<-n*(sai_1+sai_2)
  T_SR<-U/sqrt(var_cal)
  return(list(T_SR, U, var_cal)) 
}

# numerator of T_S
score_numerator<-function(data.simu){
  data.rev<-data.sort(data.simu)
  x<-data.rev[,grepl("model",names(data.rev))] #only those has model
  n<-dim(data.rev)[1]
  data_sub<-data.frame(t=data.rev$t,delta=data.rev$delta,x)
  fit<-coxph(Surv(t,delta)~.,data=data_sub)
  S_0_seq<-exp(fit$coefficients %*% t(x))
  S_1_seq<-exp(fit$coefficients %*% t(x))*data.rev$I1
  S_0<-cumsum(S_0_seq[n:1])[n:1]/n
  S_1<-cumsum(S_1_seq[n:1])[n:1]/n
  U<-sum(data.rev$delta*(data.rev$I1-S_1/S_0))
  return(U)
}

# # Shao's bootstrap method + score test
score_boot<-function(data.simu,randomization,p_trt,alpha=0.05,boot_n=200){
  n<-dim(data.simu)[1]
  U<-score_numerator(data.simu)
  strata_z<-data.simu[,grepl("strata",names(data.simu))] 
  minimization_z<-data.simu[,grepl("minimization",names(data.simu))] 
  U.boot<-numeric(boot_n)
  #boostrap variance
  boot<-1
  while(boot<=boot_n){
    boot_ind<-sample(1:n, n, replace=TRUE)
    data.boot<-data.simu[boot_ind,]
    strata_z.boot<-strata_z[boot_ind,]
    minimization_z.boot<-minimization_z[boot_ind,]
    I.boot<-treatment_assignment(n,strata_z.boot,minimization_z.boot,randomization,p_trt)
    data.boot$I1<-I.boot
    data.boot$I0<-1-I.boot
    tt <- tryCatch(score_numerator(data.boot),error=function(e) e, warning=function(w) w)
    if(is(tt,"warning")) {next} # if warning exists, skip to the next run
    U.boot[boot]<-tt
    boot<-boot+1
  }
  var.boot<-var(U.boot)
  T_SB<-U/sqrt(var.boot)
  return(T_SB)
}

# bootstrap to construct reference distribution
score_boot_distribution<-function(data.simu,randomization,p_trt,alpha=0.05,boot_n=200){
  n<-dim(data.simu)[1]
  U<-score_numerator(data.simu)
  strata_z<-data.simu[,grepl("strata",names(data.simu))] 
  minimization_z<-data.simu[,grepl("minimization",names(data.simu))] 
  U.boot<-numeric(boot_n)
  #boostrap variance
  boot<-1
  while(boot<=boot_n){
    boot_ind<-sample(1:n, n, replace=TRUE)
    data.boot<-data.simu[boot_ind,]
    strata_z.boot<-strata_z[boot_ind,]
    minimization_z.boot<-minimization_z[boot_ind,]
    I.boot<-treatment_assignment(n,strata_z.boot,minimization_z.boot,randomization,p_trt)
    data.boot$I1<-I.boot
    data.boot$I0<-1-I.boot
    tt <- tryCatch(score_numerator(data.boot),error=function(e) e, warning=function(w) w)
    if(is(tt,"warning")) {next} # if warning exists, skip to the next run
    U.boot[boot]<-tt
    boot<-boot+1
  }
  U.boot.ordered<-sort(U.boot)
  c1<-U.boot.ordered[ceiling(boot_n*alpha/2)]
  c2<-U.boot.ordered[ceiling(boot_n*(1-alpha/2))+1]
  TS_rej<-1*(U<c1 | U>c2)
  return(TS_rej)
}

### Stratified permutation + score test (Bugni et al., 2018)
score_permutation<-function(data.simu,alpha=0.05,boot_n=200){
  U<-score_numerator(data.simu)
  strata_z<-data.simu[,grepl("strata",names(data.simu))] 
  z<-strata_z
  score.p<-numeric(boot_n)
  n_strata<-dim(z)[2]
  n<-dim(z)[1]
  data.tmp<-data.simu
  for(boot in 1:boot_n){
    for(s in 1:n_strata){
      data.tmp$I1[strata_z[,s]==1]<-sample(data.simu$I1[strata_z[,s]==1]) # permute within strata
    }
    data.tmp$I0<-1-data.tmp$I1
    score.p[boot]<-score_numerator(data.tmp)
  }
  score.p.ordered<-sort(score.p)
  c1<-score.p.ordered[ceiling(boot_n*alpha/2)]
  c2<-score.p.ordered[ceiling(boot_n*(1-alpha/2))+1]
  score_SP_rej<-1*(U<c1 | U>c2)
  return(score_SP_rej)
}

# randomization + score test 
score_rand<-function(data.simu,randomization,p_trt,boot_n=200){
  n<-dim(data.simu)[1]
  U<-score_numerator(data.simu)
  strata_z<-data.simu[,grepl("strata",names(data.simu))] 
  minimization_z<-data.simu[,grepl("minimization",names(data.simu))] 
  S.rand<-numeric(boot_n)
  data.rand<-data.simu
  for(boot in 1:boot_n){
    I.rand<-treatment_assignment(n,strata_z,minimization_z,randomization,p_trt)
    data.rand$I1<-I.rand
    data.rand$I0<-1-I.rand
    S.rand[boot]<-score_numerator(data.rand)
  }
  var.rand<-var(S.rand)
  score_rand<-U/sqrt(var.rand)
  return(score_rand)
}

# randomization test to construct reference distribution for score test 
score_rand_distribution<-function(data.simu,randomization,p_trt,alpha=0.05,boot_n=200){
  n<-dim(data.simu)[1]
  U<-score_numerator(data.simu)
  strata_z<-data.simu[,grepl("strata",names(data.simu))] 
  minimization_z<-data.simu[,grepl("minimization",names(data.simu))] 
  U.rand<-numeric(boot_n)
  data.rand<-data.simu
  for(boot in 1:boot_n){
    I.rand<-treatment_assignment(n,strata_z,minimization_z,randomization,p_trt)
    data.rand$I1<-I.rand
    data.rand$I0<-1-I.rand
    U.rand[boot]<-score_numerator(data.rand)
  }
  U.ordered<-sort(U.rand)
  c1<-U.ordered[ceiling(boot_n*alpha/2)]
  c2<-U.ordered[ceiling(boot_n*(1-alpha/2))+1]
  TW_rej<-1*(U<c1 | U>c2)
  return(TW_rej)
}


