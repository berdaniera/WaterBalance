library(data.table)
library(mvtnorm)
library(msm)
a.n <- function(x) as.numeric(x)

############
# Bring in data
dat <- read.csv("data.csv")
design <- read.csv("design.csv")

# assign leaf areas and then drop them from the design matrix
# (I put them in that file so I didn't need to store them in a separate one)
Aleaf <- design$Aleaf
design <- design[,-which(colnames(design)=="Aleaf")]

dimat <- table(dat[,1:2]) # trees by droughts
ni <- length(unique(dat$i)) # number of individuals
nd <- max(dat$d) # number of days

############
# Functions for Gibbs sampling
getmuall <- function(pet,et,ao,x0){
  unlist(lapply(1:ni, function(i){ # across individuals
    unlist(lapply(1:nd, function(d){ # across droughts
      wt <- which(dat$d==d&dat$i==i)
      nt <- length(wt)
      SETtmp <- c(x0[i,d],cumsum(et[wt[-nt]]))
      ao[i,1]*pet[wt]*(1-SETtmp/ao[i,2])
    }))
  }))
}

updateOmega <- function(){
  prop <- cbind(aoG[,1],rtnorm(length(aoG[,2]),aoG[,2],jumpw,lower=minw))

  munow <- getmuall(piG,etG,exp(aoG),x0)
  munew <- getmuall(piG,etG,exp(prop),x0)

  pnow <- as.vector(tapply(dnorm(etG,munow,sqrt(sigG),log=T),dat$i,sum)) + # like
    dnorm(aoG[,2],as.matrix(design)%*%betaG,sqrt(rhoG),log=T) # prior
  pnew <- as.vector(tapply(dnorm(etG,munew,sqrt(sigG),log=T),dat$i,sum)) +
    dnorm(prop[,2],as.matrix(design)%*%betaG,sqrt(rhoG),log=T)
  a <- exp(pnew - pnow)
  z <- runif(nrow(prop),0,1)
  keep <- which(z<a,arr.ind=T)
  aoG[keep,2] <- prop[keep,2]
  list(x=aoG,accept=length(keep))
}

getVvPET <- function(pet,et,ao,x0){
  no <- length(pet)
  V <- ( (ao[1]*(1-cumsum(c(x0,et[-no]))/ao[2]))^2/sigG + 1/phiG )^(-1)
  v <- ao[1]*(1-cumsum(c(x0,et[-no]))/ao[2])*et/sigG + pet/phiG
  list(V=V,v=v)
}


updatePET <- function(){
  ao <- exp(aoG)
  Vv <- rbindlist(lapply(1:ni,function(i){
    rbindlist(lapply(which(dimat[i,]!=0),function(d){
      wg <- which(dat$d==d&dat$i==i)
      getVvPET(dat$pet[wg],
               etG[wg],
               exp(aoG[i,]),
               x0[i,d])
    }))
  }))
#  rnorm(nrow(Vv),Vv$V*Vv$v,sqrt(Vv$V))
  rtnorm(nrow(Vv),Vv$V*Vv$v,sqrt(Vv$V),etG/exp(aoG[dat$i,1])) #also try 0 lower limit
}

getVv <- function(pet,et,y,ao,x0){
  no <- length(pet)
  Va <- 1/( (1+(ao[1]*pet[-1]/ao[2])^2)/sigG + 1/tauG )
  va <- ( ao[1]*pet[-no]*(1-cumsum(c(x0,et[-c(no-1,no)]))/ao[2]) -
            (ao[1]*pet[-1]/ao[2])*(et[-1]-ao[1]*pet[-1]+ao[1]*pet[-1]*cumsum(c(x0,et[-c(no-1,no)]))/ao[2]) )/sigG + y[-no]/tauG
  # sample last value
  Vl <- 1/( 1/sigG + 1/tauG )
  vl <- (ao[1]*pet[no]*(1-sum(et[-no])/ao[2]))/sigG + y[no]/tauG
  V <- c(Va,Vl)
  v <- c(va,vl)
  list(V=V,v=v)
}


updateETall <- function(){
  Vv <- rbindlist(lapply(1:ni,function(i){
    rbindlist(lapply(which(dimat[i,]!=0),function(d){
      wg <- which(dat$d==d&dat$i==i)
      getVv(piG[wg],
            etG[wg],
            dat$y[wg],
            exp(aoG[i,]),
            x0[i,d]
      )})
    )})
  )
#  rnorm(nrow(dat),Vv$V*Vv$v,sqrt(Vv$V))
  rtnorm(nrow(dat),Vv$V*Vv$v,sqrt(Vv$V),0)
}

updateX0 <- function(){
  ao <- exp(aoG)
  Vv <- rbindlist(lapply(1:ni,function(i){
    rbindlist(lapply(which(dimat[i,]!=0),function(d){
      wg <- which(dat$d==d&dat$i==i&dat$t==1)
      V <- ( 1/(500) + ((ao[i,1]*piG[wg]/ao[i,2])^2)/sigG )^(-1)
      v <- (ao[i,1]*piG[wg] - etG[wg])*(ao[i,1]*piG[wg]/ao[i,2])/sigG
      list(V=V,v=v)
    }))
  }))
  xos <- rtnorm(nrow(Vv),Vv$V*Vv$v,sqrt(Vv$V),0)
  tmpmat <- matrix(NA,nrow=ni,ncol=nd)
  tmpmat[which(dimat!=0)] <- xos
  tmpmat
}

updateBeta <- function(x,w,Sig){
  V <- solve(crossprod(x)/Sig + priorIV)
  v <- crossprod(x,w)/Sig + priorIV%*%priorB
  t(rmvnorm(1,V%*%v,V))
}

updateSigmaIG <- function(y,mu,s1,s2){
  u1 <- s1 + length(y)/2
  u2 <- s2 + .5*sum( (y - mu)^2 )
  1/rgamma(1,u1,u2)
}


############
# Initial values
# REMEMBER TO SORT DATA BY INDIVIDUAL AND DROUGHT IDs
etG <- dat$y # best estimate is observed value
piG <- dat$pet # latent PET
betaG <- rep(0,ncol(design))
minw <- log(sapply(1:ni,function(x) max(tapply(dat$y[which(dat$i==x)],dat$d[which(dat$i==x)],sum)))*0.85) # try a fraction too....
winit <- minw + log(1.1)
aoG <- cbind(log(Aleaf),winit)

x0prior <- x0 <- x0post <- matrix(0,nrow=ni,ncol=nd)

phiG <- 1 # uncertainty on PET
rhoG <- 0.5 # regression var
sigG <- 5 # process var
tauG <- 5 # obs var
# jump parameters
jumpw <- rep(0.2,nrow(design))
# acceptance counters
ao <- 0

############
# PRIORS
priorB <- as.vector(numeric(ncol(design)))
priorIV <- solve(diag(1000,ncol(design)))
# prior on PET uncertainty
q1 <- 5000
q2 <- (0.25^2)*(q1-1)
#q1 <- q2 <- 0.001
#hist(sqrt(1/rgamma(10000,q1,q2)))
# regression error
r1 <- r2 <- 0.001
# process error - sig
s1 <- 5000#nrow(dat)#100#r1#
s2 <- (2^2)*(s1-1)#(10^2)*(nt-1)#r2#
#s1 <- s2 <- 0.001
#hist(sqrt(1/rgamma(10000,s1,s2)))
#quantile(sqrt(1/rgamma(10000,s1,s2)))
# obs error - tau
v1 <- 5000
v2 <- (3^2)*(v1-1)
#v1 <- v2 <- 0.001
#hist(sqrt(1/rgamma(10000,v1,v2)))
#mean(sqrt(1/rgamma(10000,v1,v2)))

############
# Number of samples
# This is a small number of iterations for testing, not for production use...
ng <- 1000
burns <- 100

# Store data
# Parameter names match those in Equation 5 in the manuscript
etmugib <- etgib <- petgib <- matrix(NA,ng,nrow(dat))
wgib <- matrix(NA,ng,ni)
bgib <- matrix(NA,ng,ncol(design))
taugib <- siggib <- rhogib <- phigib <- numeric(ng)

x0post <- array(0,c(ni,nd,ng))

############
# Gibbs sampler
prog <- txtProgressBar(min=0, max=ng, char="*", style=3)
for(g in 1:ng){

  piG <- updatePET()

  x0 <- updateX0()

  etG <- updateETall() # prediction of latent et

  omg <- updateOmega()
  aoG <- omg$x
  ao <- ao + omg$accept

  betaG <- updateBeta(as.matrix(design),aoG[,2],rhoG)

  mutmp <- getmuall(piG,etG,exp(aoG),x0) # model prediction for et

  tauG <- updateSigmaIG(dat$y,etG,v1,v2) # obs err
  sigG <- updateSigmaIG(etG[-which(dat$t==1)],mutmp[-which(dat$t==1)],s1,s2) # proc err
  rhoG <- updateSigmaIG(aoG[,2],as.matrix(design)%*%betaG,r1,r2) # regression error
  phiG <- updateSigmaIG(dat$pet,piG,q1,q2) # pet error

  etgib[g,] <- etG
  etmugib[g,] <- mutmp
  petgib[g,] <- piG
  wgib[g,] <- exp(aoG[,2])
  bgib[g,] <- betaG
  x0post[,,g] <- x0

  taugib[g] <- tauG
  siggib[g] <- sigG
  rhogib[g] <- rhoG
  phigib[g] <- phiG
  #adaptive jump for w?
  if(g == burns){
    jumpw <- apply(log(wgib[20:g,]),2,sd)*2.38
    if(any(jumpw==0)){
      print("Error")
      break
    }
  }
  setTxtProgressBar(prog, g)
}
close(prog)
