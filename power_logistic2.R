
######################################################################################
# Logistic Regression - Power function
# Adapated from PASS software
    #Hsieh, F.Y., Block, D.A., and Larsen, M.D. 1998. 
    #'A Simple Method of Sample Size Calculation for Linear and Logistic Regression', 
    #'Statistics in Medicine, Volume 17, pages 1623-1634.
# Author: David Aaby
# Updated: March 27 2017
######################################################################################


library(powerMediation)

alpha = .05
beta = .10
power = 1 - beta
p0 = 0.07
OR = 1.5
B = log(OR)

Pstar = .07


N = (qnorm(1-alpha/2) + qnorm(1-beta))^2 / (Pstar*(1-Pstar)*B^2)



alpha = 0.05
beta = .10
power = 1 - beta
P0 = 0.07
P1 =  0.1014493           #NEED TO SOLVE FOR THIS
OR = 1.5
B = log(OR)
R = .50
Pbar = (1-R)*P0 + R*P1

p1 = 0.07
OR = p2*(1-p1) / p1*(1-p2)
p2 = (p1*OR) / (OR*p1 + 1 - p1)

A1 = qnorm(1-alpha/2)* sqrt(Pbar*(1-Pbar)/R)
A2 = qnorm(1-beta)*sqrt(P0*(1-P0) + ((P1*(1-P1)*(1-R))/R))
A3 = ((P0-P1)^2) * (1-R)

N = ((A1 + A2)^2) / A3



############
# function #
############

LogisticSampleSize <- function(alpha, power, P0, OR, R) {
  
  #####     Inputs     #####
  # alpha = alpha level (usually 0.05)
  # power = 1 - Pr(type II error (usually .80))
  # P0 = baseline probabliity that Y=1
  # OR = odds ratio (odds1 / odds0)
  # R = percent of N with X1=1
  
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
  return(round(N))
}

LogisticSampleSize(alpha=0.05, power=.90, P0=0.07, OR=c(1.5,2.0), R=.50)


LogisticSampleSize <- function(alpha=NULL, power=NULL, P0=NULL, OR=NULL, R=NULL, N=NULL) {
  
  #####     Inputs     #####
  # alpha = alpha level (usually 0.05)
  # power = 1 - Pr(type II error (usually .80))
  # P0 = baseline probabliity that Y=1
  # OR = odds ratio (odds1 / odds0)
  # R = percent of N with X1=1
  # N = sample size
  
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
    return(round(N))
  }
  
  if(is.null(OR)) {
    beta = 1 - power
    P1 = seq(0, .90, .00001)
    N1 = N
    Pbar = (1-R)*P0 + R*P1
    n = (((qnorm(1-alpha/2)* sqrt(Pbar*(1-Pbar)/R)) + 
            (qnorm(1-beta)*sqrt(P0*(1-P0) + ((P1*(1-P1)*(1-R))/R))))^2) / 
      (((P0-P1)^2) * (1-R))
    
    OR = (P1*(1-P0)) / (P0*(1-P1))
    N.N1 = abs(n - N1)
    foo = cbind(P1,N,OR, N.N1)
    #print(head(foo))
    #foo = data.frame(foo)
    #x = foo[which(foo$N.N1==min(foo$N.N1)),]
    x = foo[which.min(foo[,4]),]
    OR = x[[3]]
    return(OR)
  }
  
  if(is.null(alpha)) {
    beta = 1 - power
    P1 = (P0*OR) / (OR*P0 + 1 - P0)
    Pbar = (1-R)*P0 + R*P1
    
    A1 = sqrt(N*((P0-P1)^2)*(1-R))
    A2 = qnorm(1-beta)*sqrt(P0*(1-P0) + ((P1*(1-P1)*(1-R))/R))
    A3 = sqrt(Pbar*(1-Pbar)/R)
    A4 = (A1 - A2) / A3
    
    alpha = 2*(1-pnorm(A4))
    return(alpha)
    
  }
  
  
  if(is.null(power)) {
    #beta = 1 - power
    P1 = (P0*OR) / (OR*P0 + 1 - P0)
    Pbar = (1-R)*P0 + R*P1
    
    A1 = sqrt(N*((P0-P1)^2)*(1-R))
    A2 = qnorm(1-alpha/2)* sqrt(Pbar*(1-Pbar)/R)
    A3 = sqrt(P0*(1-P0) + ((P1*(1-P1)*(1-R))/R))
    A4 = (A1 - A2) / A3
    
    power = pnorm(A4)
    return(power)
    
  }
  
}

LogisticSampleSize(alpha=0.05, power=.90, P0=0.07, OR=c(1.5,2.0), R=.50)     #outputs sample size N
LogisticSampleSize(alpha=0.05, power=.90, P0=0.07, R=.50, N=3326)            #outputs OR
LogisticSampleSize(power=.90, P0=0.07, OR=1.5, R=.50, N=3326)                #outputs alpha
LogisticSampleSize(alpha=0.05, P0=0.07, OR=1.5, R=.50, N=3326)               #outputs beta


alpha=0.05
power=.90
beta = 1-power
P0=0.07
P1 = seq(0,.15,.0001)
R=.50
N1= 3326
Pbar = (1-R)*P0 + R*P1
N = (((qnorm(1-alpha/2)* sqrt(Pbar*(1-Pbar)/R)) + 
        (qnorm(1-beta)*sqrt(P0*(1-P0) + ((P1*(1-P1)*(1-R))/R))))^2) / 
    (((P0-P1)^2) * (1-R))
#N-N1
#solve(2,6)

OR = (P1*(1-P0)) / (P0*(1-P1))
N.N1 = abs(N - N1)
foo = cbind(P1,N,OR, N.N1)
foo = data.frame(foo)

x= foo[which(foo$N.N1==min(foo$N.N1)),]
x= foo$OR[which.min(foo$N.N1)]
(0.1014493*(1-.070)) / (.070*(1-0.1014493))

9.997582e+00
-9.421308e+00
-2.866189e+01







#####################################################
# Version 2: output vector instead of single values #
#####################################################

LogisticSampleSize <- function(alpha=NULL, power=NULL, P0=NULL, OR=NULL, R=NULL, N=NULL) {
  
  #####     Inputs     #####
  # alpha = alpha level (usually 0.05)
  # power = 1 - Pr(type II error (usually .80))
  # P0 = baseline probabliity that Y=1
  # OR = odds ratio (odds1 / odds0)
  # R = percent of N with X1=1
  # N = sample size
  
  
  # error messages #
  if(any(alpha < 0 | alpha > 1)) stop('alpha not between 0 and 1')
  if(any(power < 0 | power > 1)) stop('power not between 0 and 1')
  if(any(P0 < 0 | P0 > 1))       stop('P0 not between 0 and 1')
  if(any(R < 0 | R > 1))         stop('R not between 0 and 1')
  if(any(OR < 0))                stop('OR not a positive value')
  
  
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
    #return(round(N))
    N = round(N)
    
  }
  
  if(is.null(OR)) {
    beta = 1 - power
    P1 = seq(0, .90, .00001)
    N1 = N
    Pbar = (1-R)*P0 + R*P1
    n = (((qnorm(1-alpha/2)* sqrt(Pbar*(1-Pbar)/R)) + 
            (qnorm(1-beta)*sqrt(P0*(1-P0) + ((P1*(1-P1)*(1-R))/R))))^2) / 
      (((P0-P1)^2) * (1-R))
    
    OR = (P1*(1-P0)) / (P0*(1-P1))
    N.N1 = abs(n - N1)
    foo = cbind(P1,N,OR, N.N1)
    #print(head(foo))
    #foo = data.frame(foo)
    #x = foo[which(foo$N.N1==min(foo$N.N1)),]
    x = foo[which.min(foo[,4]),]
    OR = x[[3]]
    #print(OR)
    #return(OR)
  
  }
  
  if(is.null(alpha)) {
    beta = 1 - power
    P1 = (P0*OR) / (OR*P0 + 1 - P0)
    Pbar = (1-R)*P0 + R*P1
    
    A1 = sqrt(N*((P0-P1)^2)*(1-R))
    A2 = qnorm(1-beta)*sqrt(P0*(1-P0) + ((P1*(1-P1)*(1-R))/R))
    A3 = sqrt(Pbar*(1-Pbar)/R)
    A4 = (A1 - A2) / A3
    
    alpha = 2*(1-pnorm(A4))
    #return(alpha)
    
  }
  
  
  if(is.null(power)) {
    #beta = 1 - power
    P1 = (P0*OR) / (OR*P0 + 1 - P0)
    Pbar = (1-R)*P0 + R*P1
    
    A1 = sqrt(N*((P0-P1)^2)*(1-R))
    A2 = qnorm(1-alpha/2)* sqrt(Pbar*(1-Pbar)/R)
    A3 = sqrt(P0*(1-P0) + ((P1*(1-P1)*(1-R))/R))
    A4 = (A1 - A2) / A3
    
    power = pnorm(A4)
    #return(power)
    
  }
  
  #print(OR)
  
  #results = c(alpha, power, P0, OR, R, N)
  #results = data.frame(results)
  #results = t(results)
  #colnames(results) = c("alpha", "power", "P0", "OR", "R" ,"N")
  
  results = NULL
  
  if(length(OR) > 1) {
    for(i in 1:length(OR)) {
      results = cbind(alpha, power, P0, OR, R, N)
    }
  }
  
  else {
    results = c(alpha, power, P0, OR, R, N)
    results = data.frame(results)
    results = t(results)
    colnames(results) = c("alpha", "power", "P0", "OR", "R" ,"N")
  }
  
  #for(i in 1:length(OR)) {
  #  results = cbind(alpha, power, P0, OR, R, N)
  #}
  
  
  return(results)
  
}

LogisticSampleSize(alpha=0.05, power=.90, P0=0.07, OR=1.5, R=.50)            #outputs sample size
LogisticSampleSize(alpha=0.05, power=.90, P0=0.07, OR=c(1.5,2.0), R=.50)     #outputs sample size N
LogisticSampleSize(alpha=0.05, power=.90, P0=0.07, R=.50, N=3326)            #outputs OR
LogisticSampleSize(power=.90, P0=0.07, OR=1.5, R=.50, N=3326)                #outputs alpha
LogisticSampleSize(alpha=0.05, P0=0.07, OR=1.5, R=.50, N=3326)               #outputs beta


LogisticSampleSize(alpha=5, power=.90, P0=0.07, OR=1.5, R=.50)       #produces error message (it is supposed to)




#####################################################
# Version 3: accept one input as a vector of values #
#####################################################

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
  if(length(l.mylist[l.mylist==2]) > 1) stop('Only vary one parameter at a time')
  
  
 
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
    Pbar = (1-R)*P0 + R*P1
    n = (((qnorm(1-alpha/2)* sqrt(Pbar*(1-Pbar)/R)) + 
            (qnorm(1-beta)*sqrt(P0*(1-P0) + ((P1*(1-P1)*(1-R))/R))))^2) / 
      (((P0-P1)^2) * (1-R))
    
    OR = (P1*(1-P0)) / (P0*(1-P1))
    N.N1 = abs(n - N1)
    foo = cbind(P1,N,OR, N.N1)
    x = foo[which.min(foo[,4]),]
    OR = x[[3]]
    
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
LogisticSampleSize(alpha=c(.05,.10), power=.90, P0=0.07, OR=c(1.5,2.0), R=.50)       #produces error message (it is supposed to)


# try using numbers from NIJ Parental Incar Grant #
LogisticSampleSize(alpha=.05, power=.80, P0=seq(.05,.50,.05), R=.54, N=300)


alpha=.05
power=.80
P0=seq(.05,.25,.05)
R=.54
N=300


beta = 1 - power
P1 = seq(0, .90, .01)
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
  x = foo[which.min(foo[,4]),]
  or = x[[3]]
  print(or)
  OR = c(OR, or)
}
OR








###############################################################################
# Version 4: correct for when P0 is input as a vector instead of single value #
###############################################################################

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



