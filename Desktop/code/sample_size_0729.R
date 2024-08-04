
# library(survival)


##############################################################################
##############################        Parameter settings         #############
##############################################################################
#------------------------------------Example----------------------------------
# #theta=0.3; Pow=0.8
# #randomization=c("CABC") #"SR","CABC","permuted_block","minimization","urn"
# lambda0 <- 0.05 
# p_trt <- 1/2; alpha=0.05
# tau1=20; tau_a=25; tau=35  # C~U(0,tau); W~U(0, tau_a) 56 70 112
# 
# #-- 4 groups 
# # prb <- c(0.4, 0.3, 0.2, 0.1)  # multinomial distribution
# # beta_c <- c(1, 1, 1)
# #-- 3 groups 
# prb <- c(0.5, 0.3, 0.2)  # multinomial distribution
# beta_c <- c(1, 1)

######################################################
####          Sample size  (original)             ####
######################################################

#---------------calculate original n*_0 for Score test ------

SIZE <- function(theta, Pow ,lambda0, p_trt, alpha, tau1, tau_a, tau, prb, beta_c){
  
    
    p <- length(beta_c)
    sigma <- rep(NA,p)
    bb <- c(beta_c,0)
    
    for (i in 1:(p+1)){
      EYt1 <- function(t) prb[i]* lambda0 *(p_trt) * (1-p_trt) * exp(bb[i]) *
        exp( -lambda0*t*exp(bb[i]) ) *  (tau-t)/(tau_a)* (tau-t)/tau
      EYt2 <- function(t) prb[i]* lambda0 *(p_trt) * (1-p_trt) * exp(bb[i]) *
        exp( -lambda0*t*exp(bb[i]) ) *  1* (tau-t)/tau
      sigma[i] <- integrate(EYt2, lower = 0, upper = tau-tau_a)$value +
        integrate(EYt1, lower = tau-tau_a, upper = tau)$value 
      
    }
    sigma_s2 <- sum(sigma)
    
    n.score <- ((qnorm(Pow, lower.tail = T) + qnorm(1-alpha, lower.tail = T))/theta)^2/sigma_s2

    
    #--------------- # of event ------------
    
    death <- rep(NA,p+1)
    for (i in 1:(p+1)){
      death[i] <- prb[i] * (1/(tau*lambda0*exp(bb[i]))) * (exp(-lambda0*tau*exp(bb[i])) - 1)
    }
    death.rate <- 1 + sum(death)
    n.death <- n.score * death.rate   #number of events
    
    
    
    ######################################################
    ####          Sample size Re-calculation          ####
    ######################################################
    
    
    alpha0 <- 0.3
    alpha1 <- 0.0299
    B <- 1-Pow
    
    
    sigma <- event <- q1 <-rep(NA,p)
    bb <- c(beta_c,0)
    for (i in 1:(p+1)){
      EYt3 <- function(t) prb[i]* lambda0 *(p_trt) * (1-p_trt) * exp(bb[i]) *
        exp( -lambda0*t*exp(bb[i]) ) *  ( ((tau1-t)/(tau1)) * ((tau1-t)/tau1) )
      sigma[i] <- integrate(EYt3, lower = 0, upper = tau1)$value 
    }
    sigma_s2 <- sum(sigma)
    
    #---------  stage 1 -----------
    # sample size
    equ <- function(xx){
      1 - pnorm(qnorm(1-alpha1,lower.tail = T), xx, sd=1) - (1-B)/(B)* pnorm(qnorm(1-alpha0,lower.tail = T), xx, sd=1)
      #1 - pnorm(qnorm(alpha1), xx, 1) - (1-B)/(B)* pnorm(qnorm(alpha0), xx, 1)
    } 
    xi <- uniroot(equ, interval = c(-100, 100))$root
    n1 <- (xi)^2/(theta)^2/(sigma_s2)
    
    # event
    for (i in 1:(p+1)){
      Ft <- function(t) prb[i]*(1-exp(-lambda0*t*exp(bb[i])))
      event[i] <- integrate(Ft, lower = 0, upper = tau1)$value
    }
    k1 <- n1/tau1*sum(event)
    
    size <- data.frame(n = ceiling(n.score) ,
                       n1 = ceiling(n1) ,
                       Death = ceiling( n.death),
                       K1 = ceiling(k1))
  return(size)
}


########################################################
#################        Function         ##############
########################################################
SIZE(theta=0.3, Pow=0.8, lambda0=0.05, p_trt=1/2, alpha=0.05, 
     tau1=20, tau_a=25, tau=30, prb=c(0.5,0.3,0.2), beta_c=c(1,1))








