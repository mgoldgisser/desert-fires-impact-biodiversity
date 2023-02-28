## OTHER CODE FOR TEST
##
PanEstFun.Sam <-
  function(y1, y2) {
    n1 <- y1[1]
    n2 <- y2[1]
    x1 <- y1[-1]
    x2 <- y2[-1]
    D12 <- sum(x1 > 0 & x2 > 0)
    f11 <- sum(x1 == 1 & x2 == 1)
    f22 <- sum(x1 == 2 & x2 == 2)
    f1p <- sum(x1 == 1 & x2 >= 1)
    fp1 <- sum(x1 >= 1 & x2 == 1)
    f2p <- sum(x1 == 2 & x2 >= 1)
    fp2 <- sum(x1 >= 1 & x2 == 2)
    K1 <- (n1 - 1) / n1
    K2 <- (n2 - 1) / n2
    if (f2p == 0 || fp2 == 0 || f22 == 0) {
      est <- D12 + K1 * f1p * (f1p - 1) / 2 / (f2p + 1) + 
        K2 * fp1 * (fp1 - 1) / 2 / (fp2 + 1) + 
        K1 * K2 * f11 * (f11 - 1) / 4 / (f22 + 1)
    } else {
      est <- D12 + K1 * f1p^2 / 2 / f2p + K2 * fp1^2 / 2 / fp2 + 
        K1 * K2 * f11^2 / 4 / f22
    }
    return(est)
  }

Chat.Ind_Inc <- function(x, m)
{
  t <- x[1]
  x=x[-1]; x <- x[x>0]
  Q1 <- sum(x == 1)
  Q2 <- sum(x == 2)
  if(Q1>0 & Q2>0)
  {
    a=(t - 1) * Q1 / ((t - 1) * Q1 + 2 * Q2)
  }
  if(Q1>1 & Q2==0)
  {
    a=(t-1)*(Q1-1) / ( (t-1)*(Q1-1) + 2 )
  } 
  if(Q1==1 & Q2==0) {a=0}
  if(Q1==0) {a=0} 
  
  Sub <- function(m){
    if(m < t) out <- 1-sum(x / t * exp(lchoose(t - x, m)-lchoose(t - 1, m)))
    if(m == t) out <- 1-Q1/t*a
    if(m > t) out <- 1-Q1/t*a^(m-t+1)
    out
  }
  sapply(m, Sub)		
}

correct_obspi_Inc<- function(X)
{
  t = X[1]
  x = X[-1]
  Sobs <- sum(x > 0) 	
  
  f1 <- sum(x == 1) 	
  f2 <- sum(x == 2)
  if(f1>0 & f2>0){
    a=(t - 1) * f1 / ((t - 1) * f1 + 2 * f2) * f1 / t
  }
  if(f1>1 & f2==0){
    a=(t-1)*(f1-1) / ( (t-1)*(f1-1) + 2 )*f1/t
  } 
  if(f1==1 & f2==0) {a=0}
  if(f1==0 ) {a=0} 	
  b <- sum(x / t * (1 - x / t) ^ t)
  w <- a / b  			
  Prob.hat <- x / t * (1 - w * (1 - x / t) ^ t)	
  Prob.hat
}

U2_equ_Inc=function(X, Y)
{
  t1=X[1]
  t2=Y[1]
  X = X[-1]; Y = Y[-1]
  Q1.=sum(X==1 & Y>0)
  p1barhat_1=p1bar_equ_Inc(c(t1,X))
  out3=Q1./t1*(1-p1barhat_1)
  #############################
  Q11=sum(X==1 & Y==1)
  p1barhat_2=p1bar_equ_Inc(c(t2,Y))
  out4=Q11/t1/t2*(1-p1barhat_1)*(1-p1barhat_2)/p1barhat_2
  if(p1barhat_2==0){out4=0}
  output=out3+out4
  return(output) 
}

p1bar_equ_Inc=function(X)
{
  t = X[1]
  X = X[-1]
  
  Q1=sum(X==1)
  Q2=sum(X==2)
  a=1
  if(Q1>0 & Q2>0)
  {
    a=2*Q2/( (t-1)*Q1+2*Q2  )
  }
  if(Q1>1 & Q2==0)
  {
    a=2/( (t-1)*(Q1-1)+2      )
  }
  if(Q1==1 &  Q2==0){a=1}
  if(Q1==0){a=1}
  return(a)
}
