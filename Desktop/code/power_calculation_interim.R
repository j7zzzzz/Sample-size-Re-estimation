library(survival)
source('Desktop/code/all_functions.R')

########################
### Main Program     ###
########################

data_gen <- function(n,theta,randomization,p_trt,lambda0,prb,beta_c,tau,case=case){
  if(case=="caseR01"){ # z~multinomial(0.26,0.48,0.26)
    z <- rmultinom(n,1,prb)
    minimization_z <- t(z)
    n_strata <- dim(z)[1]
    strata_z <- matrix(nrow=n, ncol=n_strata)
    strata_z <- t(z)
    I <- treatment_assignment(n, strata_z, minimization_z, randomization, p_trt)
    lambda0 <- lambda0
    cumuhazard <- rexp(n)
    HR <- exp(I*theta + colSums(beta_c * z[2:3,]))
    t.star <- cumuhazard/lambda0/HR
    C <- runif(n,0,tau)
    t <- pmin(t.star, C)
    delta <- 1*(t.star<=C)
    data.simu<-data.frame(t, C, delta, I, 1-I, z[2,], z[3,], strata_z, minimization_z)[1:n,]
    names(data.simu)<- c("t", "C", "delta", "I1", "I0", "model_Z1","model_Z2",
                         paste0("strata",1:n_strata), paste0("minimization",1:n_strata))
    return(data.simu)
  } #Real case:RE01
  else if (case=="case3_R01"){ # coutinuous z #incorrect model
    n_category<-3
    tmp<-discretize_z(n,0,1,n_category) # To implement this case, kindly modify the p_cut in the discretize_z function within all_function to match the settings specified in case3_R01.
    z2_continous<-tmp$Z
    z2<-t(tmp$Z_dis)
    
    minimization_z<-t(rbind(z2))
    minimization_n_col<-ncol(minimization_z)
    
    n_strata<-dim(z2)[1]
    strata_z<-matrix(nrow=n,ncol = n_strata)
    strata_z <- t(z2)
    I<-treatment_assignment(n,strata_z,minimization_z,randomization,p_trt)
    lambda0<-lambda0
    cumuhazard<-rexp(n)
    HR<-exp( I*theta + beta_c[1]*z2_continous )
    t.star<-cumuhazard/lambda0/HR
    C<-runif(n,0,tau)
    t<-pmin(t.star,C)
    delta<-1*(t.star<=C)
    data.simu<-data.frame(t, delta, I, 1-I, z2_continous, strata_z, minimization_z)[1:n,]
    names(data.simu)<- c("t", "delta", "I1", "I0","model_z2",paste0("strata",1:n_strata),paste("minimization",1:minimization_n_col))
    return(data.simu)
  }
  else if (case=="case4_R01"){ # z1~multinomial(0.26,0.48,0.26)+z2~N(0,1) with K=3
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
    lambda0<-lambda0
    cumuhazard<-rexp(n)
    HR<-exp( I*theta + colSums(beta_c * z1[2:3,]) + 1*z2_continous^2)
    t.star<-cumuhazard/lambda0/HR
    C<-runif(n,0,tau)
    t<-pmin(t.star,C)
    delta<-1*(t.star<=C)
    data.simu<-data.frame(t, delta, I, 1-I, z1[1,], z1[2,], z2_continous, strata_z, minimization_z)[1:n,]
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
    data.simu<-data.frame(t, delta, I, 1-I, z[1,], z[2,], strata_z, minimization_z)[1:n,]
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
    t.star <- cumuhazard/lambda0/HR
    C <- runif(n,0,tau)
    t <- pmin(t.star, C)
    delta <- 1*(t.star<=C)
    data.simu<-data.frame(t, C, delta, I, 1-I, z[2,], z[3,], z[4,], strata_z, minimization_z)[1:n,]
    names(data.simu)<- c("t", "C", "delta", "I1", "I0", "model_Z1","model_Z2","model_Z3",
                         paste0("strata",1:n_strata), paste0("minimization",1:n_strata))
    return(data.simu)
  }   #Real case:Lung cancer
  else if (case=="case3_L"){ # coutinuous z #incorrect model
    n_category<-4
    tmp<-discretize_z(n,0,1,n_category) # To implement this case, kindly modify the p_cut in the discretize_z function within all_function to match the settings specified in case3_L.
    z2_continous<-tmp$Z
    z2<-t(tmp$Z_dis)
    
    minimization_z<-t(rbind(z2))
    minimization_n_col<-ncol(minimization_z)
    
    n_strata<-dim(z2)[1]
    strata_z<-matrix(nrow=n,ncol = n_strata)
    strata_z <- t(z2)
    I<-treatment_assignment(n,strata_z,minimization_z,randomization,p_trt)
    lambda0<-lambda0
    cumuhazard<-rexp(n)
    HR<-exp( I*theta + beta_c[1]*z2_continous )
    t.star<-cumuhazard/lambda0/HR
    C<-runif(n,0,tau)
    t<-pmin(t.star,C)
    delta<-1*(t.star<=C)
    data.simu<-data.frame(t, delta, I, 1-I, z2_continous, strata_z, minimization_z)[1:n,]
    names(data.simu)<- c("t", "delta", "I1", "I0","model_z2",paste0("strata",1:n_strata),paste("minimization",1:minimization_n_col))
    return(data.simu)
  }
  else if (case=="case4_L"){ # z1~multinomial(0.25,0.35,0.2,0.2)+z2~N(0,1) with K=3
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
    lambda0<-lambda0
    cumuhazard<-rexp(n)
    HR<-exp( I*theta + colSums(beta_c * z1[2:4,]) + 1*z2_continous^2)
    t.star<-cumuhazard/lambda0/HR
    C<-runif(n,0,tau)
    t<-pmin(t.star,C)
    delta<-1*(t.star<=C)
    data.simu<-data.frame(t, delta, I, 1-I, z1[1,], z1[2,], z1[3,], z2_continous,strata_z, minimization_z)[1:n,]
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
    data.simu<-data.frame(t, delta, I, 1-I, z[1,], z[2,], z[3,], strata_z,minimization_z)[1:n,]
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
    tmp<-discretize_z(n,0,1,n_category) # To implement this case, kindly modify the p_cut in the discretize_z function within all_function to match the settings specified in case3_K=3.
    z2_continous<-tmp$Z
    z2<-t(tmp$Z_dis)
    
    minimization_z<-t(rbind(z2))
    minimization_n_col<-ncol(minimization_z)
    
    n_strata<-dim(z2)[1]
    strata_z<-matrix(nrow=n,ncol = n_strata)
    strata_z <- t(z2)
    I<-treatment_assignment(n,strata_z,minimization_z,randomization,p_trt)
    lambda0<-lambda0
    cumuhazard<-rexp(n)
    HR<-exp( I*theta + beta_c[1]*z2_continous )
    t.star<-cumuhazard/lambda0/HR
    C<-runif(n,0,tau)
    t<-pmin(t.star,C)
    delta<-1*(t.star<=C)
    data.simu<-data.frame(t, delta, I, 1-I, z2_continous, strata_z, minimization_z)[1:n,]
    names(data.simu)<- c("t", "delta", "I1", "I0","model_z2",paste0("strata",1:n_strata),paste("minimization",1:minimization_n_col))
    return(data.simu)
  }
  else if (case=="case3_K=4"){ # coutinuous z #incorrect model
    n_category<-4
    tmp<-discretize_z(n,0,1,n_category) # To implement this case, kindly modify the p_cut in the discretize_z function within all_function to match the settings specified in case3_K=4.
    z2_continous<-tmp$Z
    z2<-t(tmp$Z_dis)
    
    minimization_z<-t(rbind(z2))
    minimization_n_col<-ncol(minimization_z)
    
    n_strata<-dim(z2)[1]
    strata_z<-matrix(nrow=n,ncol = n_strata)
    strata_z <- t(z2)
    I<-treatment_assignment(n,strata_z,minimization_z,randomization,p_trt)
    lambda0<-lambda0
    cumuhazard<-rexp(n)
    HR<-exp( I*theta + beta_c[1]*z2_continous )
    t.star<-cumuhazard/lambda0/HR
    C<-runif(n,0,tau)
    t<-pmin(t.star,C)
    delta<-1*(t.star<=C)
    data.simu<-data.frame(t, delta, I, 1-I, z2_continous, strata_z, minimization_z)[1:n,]
    names(data.simu)<- c("t", "delta", "I1", "I0","model_z2",paste0("strata",1:n_strata),paste("minimization",1:minimization_n_col))
    return(data.simu)
  }
  else if (case=="case4_K=3"){ # z1~multinomial(0.5,0.3,0.2)+z2~N(0,1) with K=3
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
    lambda0<-lambda0
    cumuhazard<-rexp(n)
    HR<-exp( I*theta + colSums(beta_c*z1[1:length(beta_c),]) + 1*z2_continous^2)
    t.star<-cumuhazard/lambda0/HR
    C<-runif(n,0,tau)
    t<-pmin(t.star,C)
    delta<-1*(t.star<=C)
    data.simu<-data.frame(t, delta, I, 1-I, z1[1,], z1[2,], z2_continous, strata_z, minimization_z)[1:n,]
    names(data.simu)<- c("t", "delta", "I1", "I0",paste0("model_z1",1:2),"model_z2",paste0("strata",1:n_strata),paste("minimization",1:minimization_n_col))
    return(data.simu)
  }
  else if (case=="case4_K=4"){ # z1~multinomial(0.4,0.3,0.2,0.1)+z2~N(0,1) with K=3
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
    lambda0<-lambda0
    cumuhazard<-rexp(n)
    HR<-exp( I*theta + colSums(beta_c*z1[1:length(beta_c),]) + 1*z2_continous^2)
    t.star<-cumuhazard/lambda0/HR
    C<-runif(n,0,tau)
    t<-pmin(t.star,C)
    delta<-1*(t.star<=C)
    data.simu<-data.frame(t, delta, I, 1-I, z1[1,], z1[2,], z1[3,], z2_continous, strata_z, minimization_z)[1:n,]
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
    data.simu<-data.frame(t, delta, I, 1-I, z[1,], z[2,], strata_z, minimization_z)[1:n,]
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
    data.simu<-data.frame(t, delta, I, 1-I, z[1,], z[2,], z[3,], strata_z, minimization_z)[1:n,]
    names(data.simu)<- c("t", "delta", "I1", "I0","model_Z1","model_Z2","model_Z3",paste0("strata",1:n_strata),paste0("minimization",1:minimization_n_col))
    return(data.simu)
  }
}


##############################################################################
##############################        Parameter settings         #############
##############################################################################
#------------------------------------Example----------------------------------
# n1 <- 290
# n <- 428
# theta <- 0.0
# theta_Ha <- 0.3
# Pow <- 0.8
# randomization <- c("SR") #"SR","CABC","permuted_block"
# 
# lambda0 <- 0.05 
# p_trt <- 1/2; alpha <- 0.05; alpha0 <- 0.3; alpha1 <- 0.0299
# tau1 <- 20; tau_a <- 25; tau <- 30 # C~U(0, tau); W~U(0, tau_a)
# 
#
# #-- 4 groups 
# # prb <- c(0.4, 0.3, 0.2, 0.1)  # multinomial distribution
# # beta_c <- c(1, 1, 1)
# #-- 3 groups 
# prb <- c(0.5, 0.3, 0.2)  # multinomial distribution
# beta_c <- c(1, 1)
# 
# 
# case <- "case1" # Simulation:"case1","case2","case3_K=3","case3_K=4","case4_K=3","case4_K=4","case5_K=3","case5_K=4"
                  # Real case-RE01:"caseR01,"case3_R01","case4_R01","case5_R01"
                  # Real case-Lung cancer:"caseL,"case3_L","case4_L","case5_L"
# simu <- 1000

#################################################################################
#####################################        Main Program         #############
#################################################################################
SIMU_interim <- function(n, n1, theta, theta_Ha, Pow, randomization, lambda0, tau1, tau_a, tau, p_trt, alpha, alpha0, alpha1, prb, beta_c, case, simu){
  
  E_stop <- Rej <- O.Rej <- k1 <- k2 <- p1 <- p2 <- p_all <- k <- rep(NA, simu)
  E_stop_Rej <- rep(NA, simu) # E_stop_Rej=1:In early stopping, representing Rej[i]=1
  n2 <- rep(0, simu)
  D <- rep(NA, simu)
  T1 <- rep(NA, simu)
  C_alpha <- exp(-qchisq(alpha, df=4, lower.tail = F)/2)

  for (i in 1:simu){
  
    data.simu <- Calender_date <-  D1 <- D2 <- D_all <- TS1 <- TS2 <- TS <- TS3 <- Scor1 <- Scor2 <- NA
    data.simu <- data_gen(n+n1+2000, theta=theta, randomization=randomization, beta_c=beta_c,
                        lambda0=lambda0, prb=prb, tau=tau, p_trt=p_trt, case=case)
  
  ### ----------------- original data include  -------------
  
    Calender_date <- c(sort(runif(n, 0, tau_a)) )  # W: uniformly inclusion
  
    D_all <- data.simu[1:n,]
    D_all$delta <- D_all$delta * (D_all$t < tau - Calender_date)  
    D_all$t <- pmin(D_all$t, tau - Calender_date)
    TS <- as.numeric( score_robust (D_all,randomization,p_trt))  #return(list(T_SR, U, var_cal))
    p_all[i]  <- pnorm( TS[1], lower.tail = FALSE)
  
    k[i] <- sum(D_all$delta)
  
  ### ----------------- D1  for first stage-------------
  
    Calender_date_1 <- c(sort(runif(n1, 0, tau1)))  # W: uniformly inclusion
  
    D1 <- data.simu[1:n1,]
    D1$delta <- D1$delta * (D1$t < tau1 - Calender_date_1)  
    D1$t <- pmin(D1$t, tau1 - Calender_date_1)
  
  
    ### proposed T_RS in Ye and Shao, 2020, jrssb. ###
    TS1 <- as.numeric( score_robust (D1,randomization,p_trt))  #return(list(T_SR, U, var_cal))
    Scor1 <- as.numeric(score(D1))
  
    p1[i] <- pnorm( TS1[1], lower.tail = FALSE )
    sigma_s2_D1 <- TS1[3]/n1
    k1[i] <- sum(D1$delta)
  
    if (!is.na(p1[i]) & p1[i] >= alpha0) { E_stop[i] <- 1; Rej[i] <- 0; E_stop_Rej[i] <- 0}
    if (!is.na(p1[i]) & p1[i] <= alpha1) { E_stop[i] <- 1; Rej[i] <- 1; E_stop_Rej[i] <- 1}
    if (is.na(p1[i])) { E_stop[i] <- NA; Rej[i] <- NA; E_stop_Rej[i] <- NA}
  
  
  ### -------------PART 2 ----------interim analysis  ------------- 
  
    # calculate n2 -------
    if(theta!=0 & !is.na(p1[i]) & alpha1<p1[i] & p1[i]<alpha0){
      C_alpha_new <- C_alpha/p1[i]
      if(C_alpha/p1[i]>=1) C_alpha_new=0.9999999
      if(C_alpha/p1[i]<=0) C_alpha_new=0.0000001 
      n2[i] <- ((qnorm(Pow, lower.tail = T) + qnorm(1-C_alpha_new, lower.tail = T))/theta_Ha)^2/(sigma_s2_D1)
      n2[i] <- ceiling(n2[i])
      if(n2[i] < 0) {n2[i] = 0}
      }
    if(theta==0 & alpha1<p1[i] & p1[i]<alpha0) {
      C_alpha_new <- C_alpha/p1[i]
      if(C_alpha/p1[i]>=1) C_alpha_new=0.9999999
      if(C_alpha/p1[i]<=0) C_alpha_new=0.0000001 
      n2[i] <- ((qnorm(Pow, lower.tail = T) + qnorm(1-C_alpha_new, lower.tail = T))/theta_Ha)^2/sigma_s2_D1
      if(n2[i] < 0) {n2[i] = 0}
      }
  
  
  
    if(n2[i] > 0 & !is.na(p1[i]) & alpha1<p1[i] & p1[i]<alpha0) {
      D2 <- data.simu[1:(n1+n2[i]),]
      Calender_date <- c(Calender_date_1, sort(runif(n2[i], tau1, tau_a)) )  # W: uniformly inclusion
      D2$delta <- D2$delta * (D2$t < tau - Calender_date[1:(n1+n2[i])])  
      D2$t <- pmin(D2$t, tau - Calender_date[1:(n1+n2[i])])
      k2[i] <- sum(D2$delta)
      
      TS2 <- as.numeric(score_robust(D2,randomization,p_trt))  
      Scor2 <- as.numeric(score(D2))

      ### T2=T_RS-T1 ### 
      TS3 = (TS2[2] - TS1[2]) / sqrt( Scor2[3]-Scor1[3] )
      
      # ------- Confirming asymptotic independence of two-stage test statistics -------
      T1[i] <- sqrt(1/n1)*TS1[2] # test statistic of stage1 
      D[i] <- sqrt(1/(n1+n2[i]))*TS2[2] - sqrt(1/n1)*TS1[2] # test statistic of stage2
      # -------------------------------------------------------------------------------
      
      p2[i] <- pnorm(TS3, lower.tail = FALSE)
      
      E_stop[i] <- 0
      if (!is.na(p2[i]) & p1[i]*p2[i] <= C_alpha) Rej[i] <- 1
      if (!is.na(p2[i]) & p1[i]*p2[i]  > C_alpha) Rej[i] <- 0
      if (is.na(p2[i])) Rej[i] <- NA
    
      }
  
  
  ##-------- original p values ------
  O.Rej[i] <- 0
  if(p_all[i] <= alpha)  {O.Rej[i] <- 1}
  
}

id <- which(E_stop == 0)
Result <- data.frame(Design = randomization,
                     n = n, 
                     n1 = n1,
                     All_Power = round(mean(O.Rej, na.rm = TRUE), 4),   # without interim analysis 
                     death = round(mean(k, na.rm = TRUE), 0),           # # of event for all trial
                     K1 = round(mean(k1, na.rm = TRUE), 0),             # # of event for first stage
                     interim_Power = round(mean(Rej, na.rm = TRUE), 4), # with interim analysis 
                     ASN = ceiling(mean(n1+n2, na.rm = TRUE)),
                     S.D._all = round(sd(n1+n2, na.rm = TRUE), 0),   
                     ASN_stage2 = ceiling(mean(n2[id], na.rm = TRUE)),  # Calculate the average of the sample size for n2 when E_stop=0 only：The n2 for trials that entered stage2 divided by the number of trials that entered stage2
                     S.D._stage2 = round(sd(n2[id], na.rm = TRUE), 0),  # Calculate the standard deviation of the sample size for n2 when E_stop=0 only
                     Early = round(mean(E_stop, na.rm = TRUE), 4),      # the proportion of trials terminated early among the 3000 simulations
                     Early_RejH0 = round(mean(E_stop_Rej, na.rm = TRUE), 4), # Power of early termination tests
                     stage2_Power = round(mean(Rej[id], na.rm = TRUE), 4))   # Power of stage2


print(paste("death.rate=", round(round(mean(k, na.rm = TRUE), 0)/n, 4) , "|death = ", round(mean(k, na.rm = TRUE), 0)))
print(paste("n1=", n1,"| n=",n , "| theta=", theta, "| Power=", Pow, "| Design=", randomization))
return(Result)

}


########################################################
#################        Function         ##############
########################################################
SIMU_interim(n=428, n1=290, theta=0.3, theta_Ha=0.3, Pow=0.8, randomization=c("CABC"), lambda0=0.05, tau1=20, tau_a=25, tau=30, 
     p_trt=1/2, alpha=0.05, alpha0=0.3, alpha1=0.0299, prb=c(0.5, 0.3, 0.2), beta_c=c(1, 1), case="case1", simu=1000)
