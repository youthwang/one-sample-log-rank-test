
alpha <-  0.05   # Type I error rate
side  <-  1      # number of side test
pw    <-  0.8    # Desired Power
M0    <-  3      # Median survival time for historical control
HR    <-  NULL   # Either M1 OR HR needed,not need for both
K     <-  1.05   # Assumed shape parameter for weibo survival distribution.
M1    <-  6      # Median survival time assumed for study cohort
t.a   <-  12     # total accrual time
t.f   <-  12     # total follow-up time

Onesamplelogrank<-function(alpha,side,pw,M0,M1,K,t.a,t.f,HR){
                  if (is.null(M1)) M1<-(M0^K/HR)^(1/K) else HR<-(M0/M1)^K
                  f1p0 <- function(x){exp(-log(2)*(x/M1)^K)*K*log(2)*x^(K-1)/M0^K}
                  f2p0 <- function(x){(t.a+t.f-x)*exp(-log(2)*(x/M1)^K)*K*log(2)*x^(K-1)/M0^K}
                  p0   <- integrate(f1p0,0,t.f)$value+(1/t.a)*integrate(f2p0,t.f,(t.f+t.a))$value
                  f1p1 <- function(x){exp(-log(2)*(x/M1)^K)*K*log(2)*x^(K-1)/M1^K}
                  f2p1 <- function(x){(t.a+t.f-x)*exp(-log(2)*(x/M1)^K)*K*log(2)*x^(K-1)/M1^K}
                  p1   <- integrate(f1p1,0,t.f)$value+(1/t.a)*integrate(f2p1,t.f,(t.f+t.a))$value
                  f1p00<- function(x){exp(-log(2)*(x/M1)^K)*K*log(2)^2*x^(2*K-1)/M0^(2*K)}
                  f2p00<- function(x){(t.a+t.f-x)*exp(-log(2)*(x/M1)^K)*K*log(2)^2*x^(2*K-1)/M0^(2*K)}
                  p00  <- integrate(f1p00,0,t.f)$value+(1/t.a)*integrate(f2p00,t.f,(t.f+t.a))$value
                  f1p01<- function(x){exp(-log(2)*(x/M1)^K)*K*log(2)^2*x^(2*K-1)/(M0*M1)^K}
                  f2p01<- function(x){(t.a+t.f-x)*exp(-log(2)*(x/M1)^K)*K*log(2)^2*x^(2*K-1)/(M0*M1)^K}
                  p01  <- integrate(f1p01,0,t.f)$value+(1/t.a)*integrate(f2p01,t.f,(t.f+t.a))$value
                  sigma.s<-p1-p1^2+2*p00-p0^2-2*p01+2*p0*p1
                  sigma0.s<-p0
                  sigma1.s<-p1
                  omiga<-sigma1.s-sigma0.s
                  N <-(sqrt(sigma0.s)*qnorm(1-alpha/side)+sqrt(sigma.s)*qnorm(pw))^2/omiga^2
                  n<-ceiling(N)
                  Power.actual<-pnorm(-sqrt(sigma0.s/sigma.s)*qnorm(1-alpha/side)-omiga*sqrt(n/sigma.s))
                  return(data.frame(N.total=n,Events=round(n*p1),Accr.time=t.a,FU.time=t.f,HR=HR,M.cntl=M0,M.new=round(M1,2),
                        Weibo.shape=K,Alpha=alpha,Actural.power=round(Power.actual,4),Event.prob=round(p1,4)))
                  }
Onesamplelogrank(alpha,side,pw,M0,M1,K,t.a,t.f,HR)
