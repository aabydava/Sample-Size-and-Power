
######################################################################################
# Logistic Regression - Power function
# Adapated from PASS software
    #Hsieh, F.Y., Block, D.A., and Larsen, M.D. 1998. 
    #'A Simple Method of Sample Size Calculation for Linear and Logistic Regression', 
    #'Statistics in Medicine, Volume 17, pages 1623-1634.
# Author: David Aaby
# Updated: March 27 2017
######################################################################################


# Sample size Calculator for Logistic Regression #
# can solve for alpha, power, OR, N

LogisticSampleSize <- function(alpha=NULL, power=NULL, P0=NULL, OR=NULL, R=NULL, N=NULL) {
  
  #####     Inputs     #####
  # alpha = alpha level (usually 0.05)
  # power = 1 - Pr(type II error) (usually .80)
  # P0 = baseline probabliity that Y=1
  # OR = odds ratio (odds1 / odds0)
  # R = percent of N with X1=1
  # N = sample size
  
  
  ##################
  # error messages #
  ##################
  if(any(alpha < 0 | alpha > 1)) stop('alpha not between 0 and 1')
  if(any(power < 0 | power > 1)) stop('power not between 0 and 1')
  if(any(P0 < 0 | P0 > 1))       stop('P0 not between 0 and 1')
  if(any(R < 0 | R > 1))         stop('R not between 0 and 1')
  if(any(OR < 0))                stop('OR not a positive value')
  
  mylist = list(alpha, power, P0, OR, R, N)
  l.mylist = lengths(mylist)
  if(length(l.mylist[l.mylist>=2]) > 1) stop('Only vary one parameter at a time')
  
  
  
  ############################
  # sample size calculations #
  ############################
  
  
  if(is.null(alpha)) {
    beta = 1 - power
    P1 = (P0*OR) / (OR*P0 + 1 - P0)
    Pbar = (1-R)*P0 + R*P1
    
    A1 = sqrt(N*((P0-P1)^2)*(1-R))
    A2 = qnorm(1-beta)*sqrt(P0*(1-P0) + ((P1*(1-P1)*(1-R))/R))
    A3 = sqrt(Pbar*(1-Pbar)/R)
    A4 = (A1 - A2) / A3
    
    alpha = 2*(1-pnorm(A4))
    
    # output results #
    results = NULL
    
    if(length(power) > 1 | length(P0) > 1 | length(OR) > 1 | length(R) > 1 | length(N) > 1) {
      rownum = max(length(alpha), length(power), length(P0), length(OR), length(R))
      for(i in 1:rownum) {
        results = cbind(alpha, power, P0, OR, R, N)
      }
    }
    
    else {
      results = c(alpha, power, P0, OR, R, N)
      results = data.frame(results)
      results = t(results)
      colnames(results) = c("alpha", "power", "P0", "OR", "R" ,"N")
    }
    
    #return(results)
    
    
  }
  
  
  
  if(is.null(power)) {
    
    P1 = (P0*OR) / (OR*P0 + 1 - P0)
    Pbar = (1-R)*P0 + R*P1
    
    A1 = sqrt(N*((P0-P1)^2)*(1-R))
    A2 = qnorm(1-alpha/2)* sqrt(Pbar*(1-Pbar)/R)
    A3 = sqrt(P0*(1-P0) + ((P1*(1-P1)*(1-R))/R))
    A4 = (A1 - A2) / A3
    
    power = pnorm(A4)
    
    # output results #
    results = NULL
    
    if(length(alpha) > 1 | length(P0) > 1 | length(OR) > 1 | length(R) > 1 | length(N) > 1) {
      rownum = max(length(alpha), length(power), length(P0), length(OR), length(R))
      for(i in 1:rownum) {
        results = cbind(alpha, power, P0, OR, R, N)
      }
    }
    
    else {
      results = c(alpha, power, P0, OR, R, N)
      results = data.frame(results)
      results = t(results)
      colnames(results) = c("alpha", "power", "P0", "OR", "R" ,"N")
    }
    
    #return(results)
    
    
  }
  
  
  
  
  
  if(is.null(OR)) {
    
    beta = 1 - power
    P1 = seq(0, .90, .00001)
    N1 = N
    OR = NULL
    for(i in 1:length(P0)) {
      Pbar = (1-R)*P0[i] + R*P1
      n = (((qnorm(1-alpha/2)* sqrt(Pbar*(1-Pbar)/R)) + 
              (qnorm(1-beta)*sqrt(P0[i]*(1-P0[i]) + ((P1*(1-P1)*(1-R))/R))))^2) / 
        (((P0[i]-P1)^2) * (1-R))
      
      or = (P1*(1-P0[i])) / (P0[i]*(1-P1))
      N.N1 = abs(n - N1)
      foo = cbind(P1,N,or, N.N1)
      foo = foo[which(foo[,3] >= 1),]
      x = foo[which.min(foo[,4]),]
      or = x[[3]]
      #print(or)
      OR = c(OR, or)
      OR = round(OR,3)
    }
    
    # output results #
    results = NULL
    
    if(length(alpha) > 1 | length(power) > 1 | length(P0) > 1 | length(R) > 1 | length(N) > 1) {
      rownum = max(length(alpha), length(power), length(P0), length(OR), length(R))
      for(i in 1:rownum) {
        results = cbind(alpha, power, P0, OR, R, N)
      }
    }
    
    else {
      results = c(alpha, power, P0, OR, R, N)
      results = data.frame(results)
      results = t(results)
      colnames(results) = c("alpha", "power", "P0", "OR", "R" ,"N")
    }
    
    #return(results)
    
    
  }
  
  
  
  if(is.null(N)) {
    
    beta = 1 - power
    
    # first, solve for P1 #
    #OR = P1(1-P0) / P0(1-P1)
    P1 = (P0*OR) / (OR*P0 + 1 - P0)
    
    #calculate Pbar #
    Pbar = (1-R)*P0 + R*P1
    
    A1 = qnorm(1-alpha/2)* sqrt(Pbar*(1-Pbar)/R)
    A2 = qnorm(1-beta)*sqrt(P0*(1-P0) + ((P1*(1-P1)*(1-R))/R))
    A3 = ((P0-P1)^2) * (1-R)
    
    N = ((A1 + A2)^2) / A3
    N = round(N)
    
    # output results #
    results = NULL
    
    if(length(alpha) > 1 | length(power) > 1 | length(P0) > 1 | length(OR) > 1 | length(R) > 1) {
      rownum = max(length(alpha), length(power), length(P0), length(OR), length(R))
      for(i in 1:rownum) {
        results = cbind(alpha, power, P0, OR, R, N)
      }
    }
    
    else {
      results = c(alpha, power, P0, OR, R, N)
      results = data.frame(results)
      results = t(results)
      colnames(results) = c("alpha", "power", "P0", "OR", "R" ,"N")
    }
    
    #return(results)
    
    
  }
  
  return(results)
  
}


LogisticSampleSize(alpha=0.05, power=.90, P0=0.07, OR=1.5, R=.50)            #outputs sample size
LogisticSampleSize(alpha=0.05, power=.90, P0=0.07, OR=c(1.5,2.0), R=.50)     #outputs sample size N
LogisticSampleSize(alpha=0.05, power=.90, P0=0.07, R=.50, N=3326)            #outputs OR
LogisticSampleSize(power=.90, P0=0.07, OR=1.5, R=.50, N=3326)                #outputs alpha
LogisticSampleSize(alpha=0.05, P0=0.07, OR=1.5, R=.50, N=3326)               #outputs beta

LogisticSampleSize(alpha=c(0.05,0.10,0.15,0.20), power=.90, P0=0.07, OR=1.5, R=.50) 
LogisticSampleSize(alpha=0.05, power=.90, P0=c(0.07,0.1), R=.50, N=3300) 




LogisticSampleSize(alpha=5, power=.90, P0=0.07, OR=1.5, R=.50)       #produces error message (it is supposed to)
LogisticSampleSize(alpha=c(.05,.10,.15), power=.90, P0=0.07, OR=c(1.5,2.0,2.5), R=.50)       #produces error message (it is supposed to)

# try using numbers from NIJ Parental Incar Grant #
LogisticSampleSize(alpha=.05, power=.80, P0=seq(.05,.50,.05), R=.54, N=300)
LogisticSampleSize(alpha=.05, power=.80, P0=.05, R=.54, N=300)




alpha=.05
power=.80
P0=seq(.05,.25,.05)
R=.54
N=300
beta = 1 - power
P1 = seq(0, .90, .00001)
N1 = N
OR = NULL
for(i in 1:length(P0)) {
  Pbar = (1-R)*P0[1] + R*P1
  n = (((qnorm(1-alpha/2)* sqrt(Pbar*(1-Pbar)/R)) + 
          (qnorm(1-beta)*sqrt(P0[1]*(1-P0[1]) + ((P1*(1-P1)*(1-R))/R))))^2) / 
    (((P0[1]-P1)^2) * (1-R))
  
  or = (P1*(1-P0[1])) / (P0[1]*(1-P1))
  N.N1 = abs(n - N1)
  foo = cbind(P1,N,or, N.N1)
  foo2 = foo[which(foo[,3] >= 1),]
  x = foo2[which.min(foo2[,4]),]
  or = x[[3]]
  #print(or)
  OR = c(OR, or)
  OR = round(OR,3)
}
OR



