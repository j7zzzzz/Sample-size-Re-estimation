library(survival)
source('E:/Dropbox/Adaptive design/Covariate_adjusted/code/simulation/all_functions.R')



########################
### Main Program     ###
########################
data_gen <- function(n,theta,randomization,p_trt,lambda0,prb,beta_c,tau,case=case){
  if(case=="caseR01"){ # z~multinomial(0.5,0.3,0.2)
    z <- rmultinom(n,1,prb)
    minimization_z <- t(z)
    n_strata <- dim(z)[1]
    strata_z <- matrix(nrow=n, ncol=n_strata)
    strata_z <- t(z)
    I <- treatment_assignment(n, strata_z, minimization_z, randomization, p_trt)
    lambda0 <- lambda0
    cumuhazard <- rexp(n)
    HR <- exp(I*theta + colSums(beta_c * z[2:3,]))
    #HR <- exp( I*theta + colSums(beta_c*z[1:length(beta_c),]) )
    t.star <- cumuhazard/lambda0/HR
    C <- runif(n,0,tau)
    t <- pmin(t.star, C)
    delta <- 1*(t.star<=C)
    #data.simu<-data.frame(t, C, delta, I, 1-I, z[1,], z[2,], strata_z, minimization_z)[1:n,]
    data.simu<-data.frame(t, C, delta, I, 1-I, z[2,], z[3,], strata_z, minimization_z)[1:n,]
    names(data.simu)<- c("t", "C", "delta", "I1", "I0", "model_Z1","model_Z2",
                         paste0("strata",1:n_strata), paste0("minimization",1:n_strata))
    return(data.simu)
  } #Real case:RE01
  else if (case=="case3_R01"){ # coutinuous z #incorrect model
    n_category<-3
    tmp<-discretize_z(n,0,1,n_category) #調整all_function p_cut選caseR01_3
    z2_continous<-tmp$Z
    z2<-t(tmp$Z_dis)
    
    minimization_z<-t(rbind(z2))
    minimization_n_col<-ncol(minimization_z)
    
    n_strata<-dim(z2)[1]
    strata_z<-matrix(nrow=n,ncol = n_strata)
    strata_z <- t(z2)
    I<-treatment_assignment(n,strata_z,minimization_z,randomization,p_trt)
    lambda0<-lambda0 #log(2)/12
    cumuhazard<-rexp(n)
    HR<-exp( I*theta + beta_c[1]*z2_continous )
    t.star<-cumuhazard/lambda0/HR
    C<-runif(n,0,tau)
    t<-pmin(t.star,C)
    delta<-1*(t.star<=C)
    data.simu<-data.frame(t, delta, I, 1-I,z2_continous,strata_z,minimization_z)[1:n,]
    names(data.simu)<- c("t", "delta", "I1", "I0","model_z2",paste0("strata",1:n_strata),paste("minimization",1:minimization_n_col))
    return(data.simu)
  }
  else if (case=="case4_R01"){ # z1~multinomial(0.5,0.3,0.2)+z2~N(0,1) K=3 # z1~multinomial(0.4,0.3,0.2,0.1)+z2~N(0,1) K=3
    z1<-rmultinom(n,1,prb) 
    row.names(z1)<-paste0("z1_",1:nrow(z1))
    n_category<-3
    tmp<-discretize_z(n,0,1,n_category)
    z2_continous<-tmp$Z
    z2<-t(tmp$Z_dis)
    
    minimization_z<-t(rbind(z1,z2))
    minimization_n_col<-ncol(minimization_z)
    
    n_strata<-dim(z1)[1]*dim(z2)[1]
    strata_z<-matrix(nrow=n,ncol = n_strata)
    ind<-1
    for(i in 1:dim(z1)[1]){
      for(j in 1:dim(z2)[1]){
        strata_z[,ind]<-z1[i,]*z2[j,]
        ind<-ind+1
      }
    }
    I<-treatment_assignment(n,strata_z,minimization_z,randomization,p_trt)
    lambda0<-lambda0 #log(2)/12
    cumuhazard<-rexp(n)
    HR<-exp( I*theta + colSums(beta_c * z1[2:3,]) + 1*z2_continous^2)
    t.star<-cumuhazard/lambda0/HR
    C<-runif(n,0,tau)
    t<-pmin(t.star,C)
    delta<-1*(t.star<=C)
    data.simu<-data.frame(t, delta, I, 1-I,z1[1,],z1[2,],z2_continous,strata_z,minimization_z)[1:n,]
    names(data.simu)<- c("t", "delta", "I1", "I0",paste0("model_z1",1:2),"model_z2",paste0("strata",1:n_strata),paste("minimization",1:minimization_n_col))
    return(data.simu)
  }
  else if (case=="case5_R01"){ # non-PH (mis), C unif
    z <- rmultinom(n,1,prb)
    minimization_z <- t(z)
    n_strata <- dim(z)[1]
    strata_z <- matrix(nrow=n, ncol=n_strata)
    strata_z <- t(z)
    minimization_n_col<-ncol(minimization_z)
    
    I<-treatment_assignment(n,strata_z,minimization_z,randomization,p_trt)
    t.star<-exp( theta*I+ colSums(beta_c*z[1:length(beta_c),]) )+rexp(n,1)
    C<-runif(n,0,tau)
    t<-pmin(t.star,C)
    delta<-1*(t.star<=C)
    data.simu<-data.frame(t, delta, I, 1-I, z[1,], z[2,],strata_z,minimization_z)[1:n,]
    names(data.simu)<- c("t", "delta", "I1", "I0","model_Z1","model_Z2",paste0("strata",1:n_strata),paste0("minimization",1:minimization_n_col))
    return(data.simu)
  }
  if(case=="caseL"){ # z~multinomial(0.25,0.35,0.2,0.2)
    z <- rmultinom(n,1,prb)
    minimization_z <- t(z)
    n_strata <- dim(z)[1]
    strata_z <- matrix(nrow=n, ncol=n_strata)
    strata_z <- t(z)
    I <- treatment_assignment(n, strata_z, minimization_z, randomization, p_trt)
    lambda0 <- lambda0
    cumuhazard <- rexp(n)
    HR <- exp(I*theta + colSums(beta_c * z[2:4,]))
    #HR <- exp( I*theta + colSums(beta_c*z[1:length(beta_c),]) )
    t.star <- cumuhazard/lambda0/HR
    C <- runif(n,0,tau)
    t <- pmin(t.star, C)
    delta <- 1*(t.star<=C)
    #data.simu<-data.frame(t, C, delta, I, 1-I, z[1,], z[2,], strata_z, minimization_z)[1:n,]
    data.simu<-data.frame(t, C, delta, I, 1-I, z[2,], z[3,], z[4,], strata_z, minimization_z)[1:n,]
    names(data.simu)<- c("t", "C", "delta", "I1", "I0", "model_Z1","model_Z2","model_Z3",
                         paste0("strata",1:n_strata), paste0("minimization",1:n_strata))
    return(data.simu)
  }   #Real case:Lung cancer
  else if (case=="case3_L"){ # coutinuous z #incorrect model
    n_category<-4
    tmp<-discretize_z(n,0,1,n_category) #調整all_function p_cut選caseL_3
    z2_continous<-tmp$Z
    z2<-t(tmp$Z_dis)
    
    minimization_z<-t(rbind(z2))
    minimization_n_col<-ncol(minimization_z)
    
    n_strata<-dim(z2)[1]
    strata_z<-matrix(nrow=n,ncol = n_strata)
    strata_z <- t(z2)
    I<-treatment_assignment(n,strata_z,minimization_z,randomization,p_trt)
    lambda0<-lambda0 #log(2)/12
    cumuhazard<-rexp(n)
    HR<-exp( I*theta + beta_c[1]*z2_continous )
    t.star<-cumuhazard/lambda0/HR
    C<-runif(n,0,tau)
    t<-pmin(t.star,C)
    delta<-1*(t.star<=C)
    data.simu<-data.frame(t, delta, I, 1-I,z2_continous,strata_z,minimization_z)[1:n,]
    names(data.simu)<- c("t", "delta", "I1", "I0","model_z2",paste0("strata",1:n_strata),paste("minimization",1:minimization_n_col))
    return(data.simu)
  }
  else if (case=="case4_L"){ # z1~multinomial(0.5,0.3,0.2)+z2~N(0,1) K=3 # z1~multinomial(0.4,0.3,0.2,0.1)+z2~N(0,1) K=3
    z1<-rmultinom(n,1,prb) 
    row.names(z1)<-paste0("z1_",1:nrow(z1))
    n_category<-3
    tmp<-discretize_z(n,0,1,n_category)
    z2_continous<-tmp$Z
    z2<-t(tmp$Z_dis)
    
    minimization_z<-t(rbind(z1,z2))
    minimization_n_col<-ncol(minimization_z)
    
    n_strata<-dim(z1)[1]*dim(z2)[1]
    strata_z<-matrix(nrow=n,ncol = n_strata)
    ind<-1
    for(i in 1:dim(z1)[1]){
      for(j in 1:dim(z2)[1]){
        strata_z[,ind]<-z1[i,]*z2[j,]
        ind<-ind+1
      }
    }
    I<-treatment_assignment(n,strata_z,minimization_z,randomization,p_trt)
    lambda0<-lambda0 #log(2)/12
    cumuhazard<-rexp(n)
    HR<-exp( I*theta + colSums(beta_c * z1[2:4,]) + 1*z2_continous^2)
    t.star<-cumuhazard/lambda0/HR
    C<-runif(n,0,tau)
    t<-pmin(t.star,C)
    delta<-1*(t.star<=C)
    data.simu<-data.frame(t, delta, I, 1-I,z1[1,],z1[2,],z1[3,],z2_continous,strata_z,minimization_z)[1:n,]
    names(data.simu)<- c("t", "delta", "I1", "I0",paste0("model_z1",1:3),"model_z2",paste0("strata",1:n_strata),paste("minimization",1:minimization_n_col))
    return(data.simu)
  }
  else if (case=="case5_L"){ # non-PH (mis), C unif
    z <- rmultinom(n,1,prb)
    minimization_z <- t(z)
    n_strata <- dim(z)[1]
    strata_z <- matrix(nrow=n, ncol=n_strata)
    strata_z <- t(z)
    minimization_n_col<-ncol(minimization_z)
    
    I<-treatment_assignment(n,strata_z,minimization_z,randomization,p_trt)
    t.star<-exp( theta*I+ colSums(beta_c*z[1:length(beta_c),]) )+rexp(n,1)
    C<-runif(n,0,tau)
    t<-pmin(t.star,C)
    delta<-1*(t.star<=C)
    data.simu<-data.frame(t, delta, I, 1-I, z[1,], z[2,], z[3,],strata_z,minimization_z)[1:n,]
    names(data.simu)<- c("t", "delta", "I1", "I0","model_Z1","model_Z2","model_Z3",paste0("strata",1:n_strata),paste0("minimization",1:minimization_n_col))
    return(data.simu)
  }
  if(case=="case1"){ # z~multinomial(0.5,0.3,0.2)
    z <- rmultinom(n,1,prb)
    minimization_z <- t(z)
    n_strata <- dim(z)[1]
    strata_z <- matrix(nrow=n, ncol=n_strata)
    strata_z <- t(z)
    I <- treatment_assignment(n, strata_z, minimization_z, randomization, p_trt)
    lambda0 <- lambda0
    cumuhazard <- rexp(n)
    #HR <- exp(I*theta+1*z[1,]+1*z[2,])
    HR <- exp( I*theta + colSums(beta_c*z[1:length(beta_c),]) )
    t.star <- cumuhazard/lambda0/HR
    C <- runif(n,0,tau)
    t <- pmin(t.star, C)
    delta <- 1*(t.star<=C)
    data.simu<-data.frame(t, C, delta, I, 1-I, z[1,], z[2,], strata_z, minimization_z)[1:n,]
    names(data.simu)<- c("t", "C", "delta", "I1", "I0", "model_Z1","model_Z2",
                         paste0("strata",1:n_strata), paste0("minimization",1:n_strata))
    return(data.simu)
  }
  else if (case=="case2"){ # z~multinomial(0.4,0.3,0.2,0.1)
    z <- rmultinom(n,1,prb)
    minimization_z <- t(z)
    n_strata <- dim(z)[1]
    strata_z <- matrix(nrow=n, ncol=n_strata)
    strata_z <- t(z)
    I <- treatment_assignment(n, strata_z, minimization_z, randomization, p_trt)
    lambda0 <- lambda0
    cumuhazard <- rexp(n)
    #HR <- exp(I*theta+1*z[1,]+1*z[2,])
    HR <- exp( I*theta + colSums(beta_c*z[1:length(beta_c),]) )
    t.star <- cumuhazard/lambda0/HR
    C <- runif(n,0,tau)
    t <- pmin(t.star, C)
    delta <- 1*(t.star<=C)
    data.simu<-data.frame(t, C, delta, I, 1-I, z[1,], z[2,], z[3,], strata_z, minimization_z)[1:n,]
    names(data.simu)<- c("t", "C", "delta", "I1", "I0", "model_Z1", "model_Z2", "model_Z3",
                         paste0("strata",1:n_strata), paste0("minimization",1:n_strata))
    return(data.simu)
  }
  else if (case=="case3_K=3"){ # coutinuous z #incorrect model
    n_category<-3
    tmp<-discretize_z(n,0,1,n_category) #all_funtion的p_cut要選case3_K=3
    z2_continous<-tmp$Z
    z2<-t(tmp$Z_dis)
    
    minimization_z<-t(rbind(z2))
    minimization_n_col<-ncol(minimization_z)
    
    n_strata<-dim(z2)[1]
    strata_z<-matrix(nrow=n,ncol = n_strata)
    strata_z <- t(z2)
    I<-treatment_assignment(n,strata_z,minimization_z,randomization,p_trt)
    lambda0<-lambda0 #log(2)/12
    cumuhazard<-rexp(n)
    HR<-exp( I*theta + beta_c[1]*z2_continous )
    t.star<-cumuhazard/lambda0/HR
    C<-runif(n,0,tau)
    t<-pmin(t.star,C)
    delta<-1*(t.star<=C)
    data.simu<-data.frame(t, delta, I, 1-I,z2_continous,strata_z,minimization_z)[1:n,]
    names(data.simu)<- c("t", "delta", "I1", "I0","model_z2",paste0("strata",1:n_strata),paste("minimization",1:minimization_n_col))
    return(data.simu)
  }
  else if (case=="case3_K=4"){ # coutinuous z #incorrect model
    n_category<-4
    tmp<-discretize_z(n,0,1,n_category) #all_funtion的p_cut要選case3_K=4
    z2_continous<-tmp$Z
    z2<-t(tmp$Z_dis)
    
    minimization_z<-t(rbind(z2))
    minimization_n_col<-ncol(minimization_z)
    
    n_strata<-dim(z2)[1]
    strata_z<-matrix(nrow=n,ncol = n_strata)
    strata_z <- t(z2)
    I<-treatment_assignment(n,strata_z,minimization_z,randomization,p_trt)
    lambda0<-lambda0 #log(2)/12
    cumuhazard<-rexp(n)
    HR<-exp( I*theta + beta_c[1]*z2_continous )
    t.star<-cumuhazard/lambda0/HR
    C<-runif(n,0,tau)
    t<-pmin(t.star,C)
    delta<-1*(t.star<=C)
    data.simu<-data.frame(t, delta, I, 1-I,z2_continous,strata_z,minimization_z)[1:n,]
    names(data.simu)<- c("t", "delta", "I1", "I0","model_z2",paste0("strata",1:n_strata),paste("minimization",1:minimization_n_col))
    return(data.simu)
  }
  else if (case=="case4_K=3"){ # z1~multinomial(0.5,0.3,0.2)+z2~N(0,1) K=3 # z1~multinomial(0.4,0.3,0.2,0.1)+z2~N(0,1) K=3
    z1<-rmultinom(n,1,prb) 
    row.names(z1)<-paste0("z1_",1:nrow(z1))
    n_category<-3
    tmp<-discretize_z(n,0,1,n_category)
    z2_continous<-tmp$Z
    z2<-t(tmp$Z_dis)
    
    minimization_z<-t(rbind(z1,z2))
    minimization_n_col<-ncol(minimization_z)
    
    n_strata<-dim(z1)[1]*dim(z2)[1]
    strata_z<-matrix(nrow=n,ncol = n_strata)
    ind<-1
    for(i in 1:dim(z1)[1]){
      for(j in 1:dim(z2)[1]){
        strata_z[,ind]<-z1[i,]*z2[j,]
        ind<-ind+1
      }
    }
    I<-treatment_assignment(n,strata_z,minimization_z,randomization,p_trt)
    lambda0<-lambda0 #log(2)/12
    cumuhazard<-rexp(n)
    HR<-exp( I*theta + colSums(beta_c*z1[1:length(beta_c),]) + 1*z2_continous^2)
    t.star<-cumuhazard/lambda0/HR
    C<-runif(n,0,tau)
    t<-pmin(t.star,C)
    delta<-1*(t.star<=C)
    data.simu<-data.frame(t, delta, I, 1-I,z1[1,],z1[2,],z2_continous,strata_z,minimization_z)[1:n,]
    names(data.simu)<- c("t", "delta", "I1", "I0",paste0("model_z1",1:2),"model_z2",paste0("strata",1:n_strata),paste("minimization",1:minimization_n_col))
    return(data.simu)
  }
  else if (case=="case4_K=4"){ # z1~multinomial(0.5,0.3,0.2)+z2~N(0,1) K=3 # z1~multinomial(0.4,0.3,0.2,0.1)+z2~N(0,1) K=3
    z1<-rmultinom(n,1,prb) 
    row.names(z1)<-paste0("z1_",1:nrow(z1))
    n_category<-3
    tmp<-discretize_z(n,0,1,n_category)
    z2_continous<-tmp$Z
    z2<-t(tmp$Z_dis)
    
    minimization_z<-t(rbind(z1,z2))
    minimization_n_col<-ncol(minimization_z)
    
    n_strata<-dim(z1)[1]*dim(z2)[1]
    strata_z<-matrix(nrow=n,ncol = n_strata)
    ind<-1
    for(i in 1:dim(z1)[1]){
      for(j in 1:dim(z2)[1]){
        strata_z[,ind]<-z1[i,]*z2[j,]
        ind<-ind+1
      }
    }
    I<-treatment_assignment(n,strata_z,minimization_z,randomization,p_trt)
    lambda0<-lambda0 #log(2)/12
    cumuhazard<-rexp(n)
    HR<-exp( I*theta + colSums(beta_c*z1[1:length(beta_c),]) + 1*z2_continous^2)
    t.star<-cumuhazard/lambda0/HR
    C<-runif(n,0,tau)
    t<-pmin(t.star,C)
    delta<-1*(t.star<=C)
    data.simu<-data.frame(t, delta, I, 1-I,z1[1,],z1[2,],z1[3,],z2_continous,strata_z,minimization_z)[1:n,]
    names(data.simu)<- c("t", "delta", "I1", "I0",paste0("model_z1",1:3),"model_z2",paste0("strata",1:n_strata),paste("minimization",1:minimization_n_col))
    return(data.simu)
  }
  else if (case=="case5_K=3"){ # non-PH (mis), C unif
    z <- rmultinom(n,1,prb)
    minimization_z <- t(z)
    n_strata <- dim(z)[1]
    strata_z <- matrix(nrow=n, ncol=n_strata)
    strata_z <- t(z)
    minimization_n_col<-ncol(minimization_z)
    
    I<-treatment_assignment(n,strata_z,minimization_z,randomization,p_trt)
    t.star<-exp( theta*I+ colSums(beta_c*z[1:length(beta_c),]) )+rexp(n,1)
    C<-runif(n,0,tau)
    t<-pmin(t.star,C)
    delta<-1*(t.star<=C)
    data.simu<-data.frame(t, delta, I, 1-I, z[1,], z[2,],strata_z,minimization_z)[1:n,]
    names(data.simu)<- c("t", "delta", "I1", "I0","model_Z1","model_Z2",paste0("strata",1:n_strata),paste0("minimization",1:minimization_n_col))
    return(data.simu)
  }
  else if (case=="case5_K=4"){ # non-PH (mis), C unif
    z <- rmultinom(n,1,prb)
    minimization_z <- t(z)
    n_strata <- dim(z)[1]
    strata_z <- matrix(nrow=n, ncol=n_strata)
    strata_z <- t(z)
    minimization_n_col<-ncol(minimization_z)
    
    I<-treatment_assignment(n,strata_z,minimization_z,randomization,p_trt)
    t.star<-exp( theta*I+ colSums(beta_c*z[1:length(beta_c),]) )+rexp(n,1)
    C<-runif(n,0,tau)
    t<-pmin(t.star,C)
    delta<-1*(t.star<=C)
    data.simu<-data.frame(t, delta, I, 1-I, z[1,], z[2,], z[3,],strata_z,minimization_z)[1:n,]
    names(data.simu)<- c("t", "delta", "I1", "I0","model_Z1","model_Z2","model_Z3",paste0("strata",1:n_strata),paste0("minimization",1:minimization_n_col))
    return(data.simu)
  }
}

##############################################################################
##############################        Parameter settings         #############
##############################################################################
#------------------------------------Example----------------------------------
# n = 428
# theta = 0.0
# Pow = 0.8
# randomization=c("SR") #"SR","CABC","permuted_block","minimization","urn"
# 
# lambda0 = 0.05 
# tau_a = 25; tau = 30  # C~U(0,tau); W~U(0, tau_a)
# p_trt <- 1/2 
# 
# 
# # prb <- c(0.4, 0.3, 0.2, 0.1)  # multinomial distribution
# # beta_c <- c(1, 1, 1)
# 
# prb <- c(0.5, 0.3, 0.2)  # multinomial distribution
# beta_c <- c(1, 1)
# 
# 
# alpha = 0.05
# 
# case <- "case1"
# simu = 1000

###############################################################################
#####################################        Main Program         #############
###############################################################################
SIMU_all <- function(n, theta, Pow, randomization, lambda0, tau_a, tau, p_trt, prb, beta_c, case, simu){

alpha = 0.05

p_Wd <- p_LR <- p_StrLR <- p_stdLR <- p_RoLR <- p_TS <- p_RS <- p_TLR <- p_stdLR1 <- k <- rep(NA, simu)
C_alpha <- exp(-qchisq(alpha, df=4, lower.tail = F)/2)
  
for (i in 1:simu){
    
  data.simu <- Calender_date <- D_all <- NA
  data.simu <- data_gen(n, theta=theta, randomization=randomization, beta_c=beta_c,
                          lambda0=lambda0, prb=prb, tau=tau, p_trt=p_trt, case=case)
    
    
    ### ----------------- original data include  -------------
    
    Calender_date <- c(sort(runif(n, 0, tau_a)) )  # W: uniformly inclusion
    
    D_all <- data.simu[1:n,]
    D_all$delta <- D_all$delta * (D_all$t < tau - Calender_date)  
    D_all$t <- pmin(D_all$t, tau - Calender_date)
    
    k[i] = sum(D_all$delta)
    #------------other test----------------
    ### model based wald test ###
    Wd <- regular_wald(D_all)
    p_Wd[i] <- Wd > qnorm(alpha, lower.tail = F)
    
    ### simple log rank test by Survival Package ###   
    LR <- survdiff(Surv(t,delta)~I1,data=D_all)
    Z_LR = (LR$obs[2]-LR$exp[2]) / sqrt(LR$var[1])
    p_LR[i] = Z_LR > qnorm(alpha, lower.tail = F)

    ### stratified log-rank test T_SL ### 
    fz <- ind_to_factor(D_all)
    SLR <- survdiff(Surv(t,delta)~I1+strata(fz),data=D_all)
    Z_SLR = sum( (SLR$obs-SLR$exp)[2,] ) / sqrt(SLR$var[1])
    p_StrLR[i] = Z_SLR > qnorm(alpha, lower.tail = F)

    ### standard log-rank test statistic ###
    stdLR = wlogrank(D_all)
    p_stdLR[i] = unlist(stdLR[2]) > qnorm(alpha, lower.tail = F)
  
    ### proposed robust log-rank test ###
    logrank.test_robust <- wlogrank_robust(D_all, randomization=randomization, p_trt=p_trt)
    p_RoLR[i] <- unlist(logrank.test_robust)[1] > qnorm(alpha, lower.tail = F)

    ### score test in Lin and Wei, 1989. Ts with robust variance ###
    score.test <- score(D_all)
    p_TS[i] <- unlist(score.test)[1] > qnorm(alpha, lower.tail = F)

    ### proposed T_RS in Ye and Shao, 2020, jrssb. ###
    score.test_robust <- score_robust (D_all,randomization,p_trt)
    p_RS[i] <- unlist(score_robust (D_all,randomization,p_trt))[1] > qnorm(alpha, lower.tail = F)
  }
  
  
  Tab1 = data.frame(death = round(mean(k,na.rm = TRUE), 0),
                      RS = round(mean(p_RS,na.rm = TRUE),4),
                      Wd = round(mean(p_Wd,na.rm = TRUE),4),
                      LR = round(mean(p_LR,na.rm = TRUE),4),
                      StrLR = round(mean(p_StrLR,na.rm = TRUE),4),
                      stdLR = round(mean(p_stdLR,na.rm = TRUE),4),
                      RoLR = round(mean(p_RoLR,na.rm = TRUE),4), 
                      TS = round(mean(p_TS,na.rm = TRUE),4)
  )
  
  print(paste("death.rate=", round(round(mean(k,na.rm = TRUE), 0)/n, 4) , "|death = ", round(mean(k,na.rm = TRUE), 0)))
  print(paste("| n=",n , "| theta=", theta, "| Power=", Pow, "| Design=", randomization))
  return(Tab1)

}

########################################################
#################        Function         ##############
########################################################
SIMU(n=428, theta=0.3, Pow=0.8, randomization=c("permuted_block"), lambda0=0.05, tau_a=25, tau=30, 
     p_trt=1/2, prb=c(0.5, 0.3, 0.2), beta_c=c(1, 1), case="case1", simu=1000)



