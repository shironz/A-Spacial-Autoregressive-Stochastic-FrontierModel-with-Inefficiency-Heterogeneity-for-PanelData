#########################################################################
##  Here Try different initial values and different longer chains
##  to show that my MCMC estimation methods work well





############################################################################  
#############################  SARH  ####################################### 
############################################################################  
############################################################################  
##########################################################
#   
#   A Simple SAR SFM with exponential distr of inefficiency
#   Shirong Zhao
#
##########################################################

# y = rho*w*y + x*beta + v - u, panel data
# Model: (I_NT - rho*I_T*kronecker*W)*Y = x*beta + v - u
# u ~ exp(lambda), lambda = exp(z*gamma)
# beta ~ N(0,sigma2*R)
# sigma2 ~ gamma(as,bs)
# gamma ~ N(0,R1)
# u ~ exp(lambda)
# rho ~ Unif(1/r.min, 1)


library(MASS)
library(coda)
library(mvtnorm)
library(Matrix) 
library(mnormt)
library(lattice)
library(msm)
library(mnormt)
library(gdata) # use inside command write.fwf

# subfunctions
# source("./Functions/data_gen.ipynb") # simulate spacial data
##############################################################
####### This function generates spacial data with Y, D and W
data.gen<-function(P1,P2){
  
  V<-P1*P2
  Vid<-array(1:V,c(P1,P2)) # Arbitrarily identifies voxels   
  
  Y<-rep(-99,V) # Arbitrary response vector                 
  Dj<-rep(-99,V)
  Nj<-matrix(-99,nrow=V,ncol=8)  # This will list the neighbors for each voxel
  
  
  for(p1 in 1:P1){
    for(p2 in 1:P2){
      
      # Create the neighborhood matrix, the design matrix, and the segmentation
      # vectors
      l1<-max(c(1,p1-1))
      l2<-max(c(1,p2-1))
      u1<-min(c(P1,p1+1))
      u2<-min(c(P2,p2+1))
      
      neighbors<-as.vector(Vid[l1:u1,l2:u2])
      neighbors<-neighbors[neighbors!=Vid[p1,p2]]
      
      Y[Vid[p1,p2]]<- 10*dnorm((p1-p2),0,10) # rnorm(1,0,0.1 is the noise data, here no noise
      
      wjv<-length(neighbors)
      Dj[Vid[p1,p2]]<-wjv
      Nj[Vid[p1,p2],1:wjv]<-neighbors           
    }
  }
  # Creating the necessary sparse matrices
  i<-NULL
  j<-NULL
  
  for(v in 1:V){
    i<-c(i,rep(v,Dj[v]))
    j<-c(j,Nj[v,1:Dj[v]])
  }
  
  W<- sparseMatrix(i, j, x = 1)  
  D<-sparseMatrix(1:V,1:V,x=Dj)
  
  return(list("Y"=Y,"D"=D,"W"=W))
}

##########################################
# Generates D and W Matrix

P1 = 10
P2 = 10

data<-data.gen(P1=P1,P2=P2)
Y0<-matrix(data$Y, nrow = P1, ncol = P2) 
#levelplot(Y0)

D<-data$D
WW<-data$W 
W=WW/rowSums(WW) # weight matrix, row sums up to 1

##########################################
# Generates Simulation Data
#
set.seed(1111)
N = P1*P2
T = 5
NT = N*T
X = cbind(1,rnorm(NT,0,1))
Z = cbind(1,rnorm(NT,0,1))
beta.true = c(1,4)
sigma2.true = 1
gamma.true = c(-2,1)
lambda.true = exp(Z%*%gamma.true)
rho.true = 0.5
IW=kronecker(diag(T), W)
IRW = solve(diag(NT)-rho.true*IW) # IRW = (I_NT - rho*I_T*kronecker*W)^{-1}
v = rnorm(NT,0,sd = sqrt(sigma2.true))
u.true = rexp(NT, rate=lambda.true)
Y = IRW%*%(X%*%beta.true + v - u.true)

##########################################
# specify the priors

NT = dim(X)[1]
p = dim(X)[2] 
q = dim(Z)[2]
R = 1000*diag(p)
R1 = 100
R2 = 100
RI<-solve(R)
as = 1
bs = 1
al = 1
bl = 1
delta1 = 0.5 # need to be tuned, gamma1.p/gamma1 ~ N(gamma1, delta1)
delta2 = 0.5 # need to be tuned, gamma2.p/gamma2 ~ N(gamma2, delta2)

acc1 = 0 # calculate the accepting rate for gamma1
acc2 = 0 # calculate the accepting rate  for gamma2
thin = 5e1
iter = 1e3
G=iter/thin
Burn = iter/4

### Set the inital value
beta = c(0,0)
gamma = c(0,0)
gamma1 = gamma[1]
gamma2 = gamma[2]
sigma2 = 1
lambda = exp(Z%*%gamma)
rho = 0
u = rexp(NT, rate=lambda)


##########################################
# Create Vector/Matrix to save parameters

beta.save<-matrix(-99,nrow=iter/thin,ncol=p)
sigma2.save<-rep(-99,iter/thin)
gamma.save<-matrix(-99,nrow=iter/thin,ncol=q)
rho.save<-rep(-99,iter/thin)
u.save<-matrix(-99,nrow=iter/thin,ncol=NT)
beta.save[1,]<-beta
sigma2.save[1]<-sigma2
gamma.save[1,]<-gamma
rho.save[1]<-rho
u.save[1,]<-u


##################################################################
# calculate the range of rho, rho ~(1/r.min,1)
# where r.min is the most negative real characteristic root of W
e <- eigen(W)
r.min <- min(e$values)
r.minI <- 1/r.min


##############################################
# quantities only computed once
XTXRI<-solve(t(X)%*%X+RI)
YTWTWYI<-solve(t(Y)%*%t(IW)%*%IW%*%Y)


############################################################################  
############################################################################   
# Burn in loop

for(g in 1:Burn){
  
  
  ##############################################
  # Sample the coefficients beta
  
  Yh=(diag(NT)-rho*IW)%*%Y   # Yh = (I-rho*W)*Y
  
  cv.beta<-XTXRI*sigma2
  mu.beta<-XTXRI%*%t(X)%*%(Yh+u)
  
  
  beta<-as.vector(rmvnorm(1,mean=mu.beta,sigma=cv.beta)) 
  
  
  ##############################################
  # sample inefficiency components 
  
  
  Yt<-Yh-X%*%beta  # Yt = (I-rho*W)*Y - X%*%beta
  
  l1<-0
  u1<-Inf
  mu.u<-as.vector(-(Yt+lambda*sigma2*rep(1,NT)))
  cv.u<-as.vector(sigma2*diag(NT))
  u<-rtnorm(NT, mean=mu.u, sd=sqrt(sigma2), lower=l1, upper=u1)
  u<-as.vector(u)
  
  
  ###############################################
  # Sample sigma2
  
  
  asp<-as+NT/2+p/2
  bsp<-as.vector(bs+sum((Yh-X%*%beta+u)^2)/2+t(beta)%*%RI%*%beta/2)
  
  samp<-rgamma(1,asp,bsp)
  
  sigma2<-1/samp
  
  
  ###############################################
  # Sample gamma1
  
  gamma1.p<-rnorm(1,gamma1,sd=sqrt(delta1))
  gamma.p<-c(gamma1.p,gamma2)
  lambda.p<-exp(Z%*%gamma.p)
  llik=sum(dexp(u,lambda,log=TRUE))
  gamma1.prior<-dnorm(gamma1,0,sd=sqrt(R1),log=TRUE)
  llik.p=sum(dexp(u,lambda.p,log=TRUE))
  gamma1.p.prior<-dnorm(gamma1.p,0,sd=sqrt(R1),log=TRUE)
  
  r=exp(llik.p + gamma1.p.prior - llik - gamma1.prior)
  z<-rbinom(1,1,min(r,1))  
  
  if(z==1){
    gamma1 =gamma1.p
    gamma =gamma.p
    lambda =lambda.p
    acc1 = acc1 + 1 	
  }
  
  
  
  ###############################################
  # Sample gamma2
  
  gamma2.p<-rnorm(1,gamma2,sd=sqrt(delta2))
  gamma.p<-c(gamma1,gamma2.p)
  lambda.p<-exp(Z%*%gamma.p)
  llik=sum(dexp(u,lambda,log=TRUE))
  gamma2.prior<-dnorm(gamma2,0,sd=sqrt(R2),log=TRUE)
  llik.p=sum(dexp(u,lambda.p,log=TRUE))
  gamma2.p.prior<-dnorm(gamma2.p,0,sd=sqrt(R2),log=TRUE)
  
  r=exp(llik.p + gamma2.p.prior - llik - gamma2.prior)
  z<-rbinom(1,1,min(r,1))  
  
  if(z==1){
    gamma2 =gamma2.p
    gamma =gamma.p
    lambda =lambda.p
    acc2 = acc2 + 1 	
  }
  
  
  
  ###############################################
  # Sample rho
  
  
  l2<--1
  u2<-1
  
  mu.rho <- as.vector(YTWTWYI%*%t(Y)%*%t(IW)%*%(Y-X%*%beta+u))  
  cv.rho <- as.vector(sigma2*YTWTWYI)
  
  rho<-rtnorm(1, mean=mu.rho, sd=sqrt(cv.rho), lower=l2, upper=u2)
  
  
  ###############################################
  # Tune delta, delta1, delta2
  
  if(g%%1000==0){   
    
    delta1<- delta1 + (acc1/1000 >0.55)*0.75*delta1 - (acc1/1000 < 0.35)*0.75*delta1
    delta2<- delta2 + (acc2/1000 >0.55)*0.75*delta2 - (acc2/1000 < 0.35)*0.75*delta2
    print(c(acc1,acc2))
    
    acc1<-0
    acc2<-0
    
  }
  
}


############################################################################  
############################################################################  
# Smpling loop


for(g in (thin+1):iter){
  
  ##############################################
  # Sample the coefficients beta
  
  Yh=(diag(NT)-rho*IW)%*%Y   # Yh = (I-rho*W)*Y
  
  cv.beta<-XTXRI*sigma2
  mu.beta<-XTXRI%*%t(X)%*%(Yh+u)
  
  
  beta<-as.vector(rmvnorm(1,mean=mu.beta,sigma=cv.beta)) 
  
  
  ##############################################
  # sample inefficiency components 
  
  
  Yt<-Yh-X%*%beta  # Yt = (I-rho*W)*Y - X%*%beta
  
  l1<-0
  u1<-Inf
  mu.u<-as.vector(-(Yt+lambda*sigma2*rep(1,NT)))
  cv.u<-as.vector(sigma2*diag(NT))
  u<-rtnorm(NT, mean=mu.u, sd=sqrt(sigma2), lower=l1, upper=u1)
  u<-as.vector(u)
  
  ###############################################
  # Sample sigma2
  
  
  asp<-as+NT/2+p/2
  bsp<-as.vector(bs+sum((Yh-X%*%beta+u)^2)/2+t(beta)%*%RI%*%beta/2)
  
  samp<-rgamma(1,asp,bsp)
  
  sigma2<-1/samp
  
  
  ###############################################
  # Sample gamma1
  
  gamma1.p<-rnorm(1,gamma1,sd=sqrt(delta1))
  gamma.p<-c(gamma1.p,gamma2)
  lambda.p<-exp(Z%*%gamma.p)
  llik=sum(dexp(u,lambda,log=TRUE))
  gamma1.prior<-dnorm(gamma1,0,sd=sqrt(R1),log=TRUE)
  llik.p=sum(dexp(u,lambda.p,log=TRUE))
  gamma1.p.prior<-dnorm(gamma1.p,0,sd=sqrt(R1),log=TRUE)
  
  r=exp(llik.p + gamma1.p.prior - llik - gamma1.prior)
  z<-rbinom(1,1,min(r,1))  
  
  if(z==1){
    gamma1 =gamma1.p
    gamma =gamma.p
    lambda =lambda.p
    acc1 = acc1 + 1 	
  }
  
  
  
  ###############################################
  # Sample gamma2
  
  gamma2.p<-rnorm(1,gamma2,sd=sqrt(delta2))
  gamma.p<-c(gamma1,gamma2.p)
  lambda.p<-exp(Z%*%gamma.p)
  llik=sum(dexp(u,lambda,log=TRUE))
  gamma2.prior<-dnorm(gamma2,0,sd=sqrt(R2),log=TRUE)
  llik.p=sum(dexp(u,lambda.p,log=TRUE))
  gamma2.p.prior<-dnorm(gamma2.p,0,sd=sqrt(R2),log=TRUE)
  
  r=exp(llik.p + gamma2.p.prior - llik - gamma2.prior)
  z<-rbinom(1,1,min(r,1))  
  
  if(z==1){
    gamma2 =gamma2.p
    gamma =gamma.p
    lambda =lambda.p
    acc2 = acc2 + 1 	
  }
  
  
  
  ###############################################
  # Sample rho
  
  
  l2<--1
  u2<-1
  
  mu.rho <- as.vector(YTWTWYI%*%t(Y)%*%t(IW)%*%(Y-X%*%beta+u))  
  cv.rho <- as.vector(sigma2*YTWTWYI)
  
  rho<-rtnorm(1, mean=mu.rho, sd=sqrt(cv.rho), lower=l2, upper=u2)
  
  
  
  ###############################################
  # Save the parameters
  
  if(g%%thin==0){   
    
    print(g)
    
    beta.save[g/thin,]<-beta
    sigma2.save[g/thin]<-sigma2
    gamma.save[g/thin,]<-gamma
    rho.save[g/thin]<-rho
    u.save[g/thin,]<-u
    
  }
}


###############################################
# calculate the accepting rate


acc1/(iter-thin)
acc2/(iter-thin)



###############################################
# Summarizing results 




beta.mcmc = as.mcmc(beta.save[(G/2):G,])
print(paste0("Estimate Mean of beta:  ", apply(beta.mcmc,2,mean)))
print(paste0("True beta:  ", beta.true))
HPDinterval(beta.mcmc)
print(paste0("Effective Sample Size:  ", effectiveSize(beta.mcmc)))
plot(beta.mcmc)
autocorr.plot(beta.mcmc)




sigma2.mcmc = as.mcmc(sigma2.save[(G/2):G])
print(paste0("Estimate Mean of sigma2:  ", mean(sigma2.mcmc)))
print(paste0("True sigma2:  ", sigma2.true))
HPDinterval(sigma2.mcmc)
print(paste0("Effective Sample Size:  ", effectiveSize(sigma2.mcmc)))
plot(sigma2.mcmc)
autocorr.plot(sigma2.mcmc)




gamma.mcmc = as.mcmc(gamma.save[(G/2):G,])
print(paste0("Estimate Mean of gamma:  ", apply(gamma.mcmc,2,mean)))
print(paste0("True gamma:  ", gamma.true))
HPDinterval(gamma.mcmc)
print(paste0("Effective Sample Size:  ", effectiveSize(gamma.mcmc)))
plot(gamma.mcmc)
autocorr.plot(gamma.mcmc)




rho.mcmc = as.mcmc(rho.save[(G/2):G])
print(paste0("Estimate Mean of rho:  ", mean(rho.mcmc)))
print(paste0("True rho:  ", rho.true))
HPDinterval(rho.mcmc)
print(paste0("Effective Sample Size:  ", effectiveSize(rho.mcmc)))
plot(rho.mcmc)
autocorr.plot(rho.mcmc)





##########################################
# Save the data for SARH

SARH1=cbind(beta.save[,1], beta.save[,2], sigma2.save, gamma.save[,1], gamma.save[,2], rho.save)
formatInfoSARH1<-write.fwf(x=SARH1, file="./Data/SARH1.txt", formatInfo=TRUE)
write.csv(formatInfoSARH1, "./Data/formatInfoSARH1.csv")

############################################
# Save the graph


pdf('./Output/SARH.pdf')

beta.mcmc = as.mcmc(beta.save[(G/2):G,])
print(paste0("Estimate Mean of beta:  ", apply(beta.mcmc,2,mean)))
print(paste0("True beta:  ", beta.true))
HPDinterval(beta.mcmc)
print(paste0("Effective Sample Size:  ", effectiveSize(beta.mcmc)))
plot(beta.mcmc)
autocorr.plot(beta.mcmc)


sigma2.mcmc = as.mcmc(sigma2.save[(G/2):G])
print(paste0("Estimate Mean of sigma2:  ", mean(sigma2.mcmc)))
print(paste0("True sigma2:  ", sigma2.true))
HPDinterval(sigma2.mcmc)
print(paste0("Effective Sample Size:  ", effectiveSize(sigma2.mcmc)))
plot(sigma2.mcmc)
autocorr.plot(sigma2.mcmc)


gamma.mcmc = as.mcmc(gamma.save[(G/2):G,])
print(paste0("Estimate Mean of gamma:  ", apply(gamma.mcmc,2,mean)))
print(paste0("True gamma:  ", gamma.true))
HPDinterval(gamma.mcmc)
print(paste0("Effective Sample Size:  ", effectiveSize(gamma.mcmc)))
plot(gamma.mcmc)
autocorr.plot(gamma.mcmc)


rho.mcmc = as.mcmc(rho.save[(G/2):G])
print(paste0("Estimate Mean of rho:  ", mean(rho.mcmc)))
print(paste0("True rho:  ", rho.true))
HPDinterval(rho.mcmc)
print(paste0("Effective Sample Size:  ", effectiveSize(rho.mcmc)))
plot(rho.mcmc)
autocorr.plot(rho.mcmc)


dev.off()


################################################################
# make latex code for tables:

# true values
tab=matrix(nrow=6,ncol=5)
tab[1:2,1] = beta.true
tab[3,1] = rho.true
tab[4,1] = sigma2.true
tab[5:6,1] = gamma.true

# initial values
tab[1:2,2] = beta.save[1,]
tab[3,2] = rho.save[1]
tab[4,2] = sigma2.save[1]
tab[5:6,2] = gamma.save[1,]


# estimate
tab[1:2,3] = round(apply(beta.mcmc,2,mean),2)
tab[3,3] = round(mean(rho.mcmc),2)
tab[4,3] = round(mean(sigma2.mcmc),2)
tab[5:6,3] = round(apply(gamma.mcmc,2,mean),2)



# 95% HPD: left
tab[1:2,4] = round(HPDinterval(beta.mcmc)[1:2,1],2)
tab[3,4] = round(HPDinterval(rho.mcmc)[1],2)
tab[4,4] = round(HPDinterval(sigma2.mcmc)[1],2)
tab[5:6,4] = round(HPDinterval(gamma.mcmc)[1:2,1],2)


# 95% HPD: right
tab[1:2,5] = round(HPDinterval(beta.mcmc)[1:2,2],2)
tab[3,5] = round(HPDinterval(rho.mcmc)[2],2)
tab[4,5] = round(HPDinterval(sigma2.mcmc)[2],2)
tab[5:6,5] = round(HPDinterval(gamma.mcmc)[1:2,2],2)

print(tab)

res=formatC(tab,width=7,digits=2,format="f")

print(res)

# \\ means backslash  \ 
param = c("$\\beta_1 $",  
          "$\\beta_2 $",    
          "$\\rho    $",
          "$\\sigma^2$",   
          "$\\gamma_1$",
          "$\\gamma_2$") 


tex=paste( 
          param, "&",
          res[,1],"&",
          res[,2],"&",
          res[,3],"&","[",
          res[,4],",",
          res[,5],"]","\\\\")

print(tex)

write(tex,file="./Output/SARH_test1_inital.tex")










