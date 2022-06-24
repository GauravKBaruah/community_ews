# predator-prey eco-evolutionary dynamics function for Early warning system measurements in a simple ecological food-web!

library(entropy)


require(MASS)


# error functon
erf<-function(x) {return(2 * pnorm(x * sqrt(2)) - 1)}


# eco-evolutionary dynamic function for predator-prey module that returns species abundance time series, trait timeseries, trait variance
foodweb.model.ews<-function(tt, parms = parameters,
                            adjacency.mat,forcing.species=0)
{
  
  
  
  adjacent_matrix<-adjacency.mat
  diagonal.elements <- diag(adjacency.mat)
  n.middle.sp<- 2 
  n.basal.sp<-length(which(diagonal.elements == -1)) 
  n.consumers<-6-n.basal.sp-n.middle.sp
  
  
  middle.sp <- 2
  tt<-start.time
  a=2
  
  
  if(forcing.species==0){ 
    forcing.middle.sp <- 0
    forcing.consumer<-0
    forcing.n.basal.sp<-0
  } else if(forcing.species==4 || forcing.species==5 ){ 
    forcing.middle.sp <- forcing.species -n.basal.sp
    forcing.consumer<-0
    forcing.n.basal.sp<-0
  
    } else if(forcing.species==6){
      forcing.middle.sp <- 0  
      forcing.consumer<-forcing.species - n.middle.sp-n.basal.sp
      forcing.n.basal.sp<-0
    }  else if (forcing.species ==1 || forcing.species == 2 ||forcing.species==3){
    forcing.middle.sp <- 0
    forcing.consumer<-0
    forcing.n.basal.sp<-forcing.species
    
  }
  
  
  sa<- array(dim=c(1,(n.basal.sp+n.consumers+n.middle.sp)))
  
  gammA<-array(dim=c(tt,n.basal.sp,n.middle.sp))
  betaA<-array(dim=c(tt,n.basal.sp,n.middle.sp))
  t.gammA<-array(dim=c(tt,n.middle.sp,n.basal.sp))
  t.betaA<-array(dim=c(tt,n.middle.sp,n.basal.sp))
  
  M.gammA<-array(dim=c(tt,n.middle.sp,n.consumers))  # mutual species to predator matrix
  M.betaA<-array(dim=c(tt,n.middle.sp,n.consumers)) # mutual species to predator matrix
  dM.gammA<-array(dim=c(tt,n.consumers,n.middle.sp))  # mutual species to predator matrix
  dM.betaA<-array(dim=c(tt,n.consumers,n.middle.sp)) # mutual species to predator matrix
  
  
  det_bip<-array(dim=c(tt,n.consumers))
  z<-array(dim=c(tt,(n.basal.sp+n.consumers+n.middle.sp)))
  det_biB<-array(dim=c(tt, n.basal.sp))
  det_biM<-array(dim=c(tt,n.middle.sp))
  
  C<-cg<-array(dim=c(tt,n.consumers))
  R<-gA<-rg<-array(dim=c(tt,n.basal.sp))
  M<-mg<-array(dim=c(tt,n.middle.sp))
  gM<-array(dim=c(tt,n.middle.sp))
  
  
  c<-matrix(0, tt,n.consumers)
  ff<-matrix(0, tt,n.basal.sp)
  fc<-matrix(0, tt, n.middle.sp)
  
  b_bar<-array(dim=c(tt,n.consumers))
  h<-0.4
  gamma<-array(dim=c(tt, n.basal.sp))
  xeta<-array(dim=c(tt, n.consumers))
  
  
  R[1,]<-parms$InitialAB[1:n.basal.sp]
  M[1,]<-parms$InitialAB[(n.basal.sp+1):(n.basal.sp+n.middle.sp)]
  C[1,]<-parms$InitialAB[(n.basal.sp+n.middle.sp+1) :(n.middle.sp+n.consumers+n.basal.sp)]
  
  z[1,]<-parms$initialmu
  
  sa[1,]<-parms$variance
  
  gamma<-parms$gamma
  thetaB<-parms$theta[1]
  thetaP<-parms$theta[2]
  w<-parms$w
  dt<-0.1
  env.change<-env.function(name = parms$forcing.type)
  
  
  for (t in 1:(tt-1)) {
    
    #print(t)
    
    # growth function of the model for Plant species and animal species
    for ( j in 1: (n.basal.sp)){
      det_biB[t,j]<-0.5*(erf((thetaB-z[t,j])/(sqrt(2*sa[1,j])))+erf((thetaB+z[t,j])/(sqrt(2*sa[1,j]))))
      # this is growth rate
    }
    
    
    det_biM[t,]<-  parms$R[(n.basal.sp+1):(n.basal.sp+n.middle.sp)]
    det_bip[t,]<-  parms$R[(n.basal.sp+n.middle.sp+1):(n.basal.sp+n.consumers+n.middle.sp)]
    
    #basal species random competitive interactions
    
    
    alpha   <- parms$Pmatrix
    
    m.alpha <- parms$Amatrix
    
    
    for (i in 1: n.basal.sp){
      for (k in (n.basal.sp+1): (n.basal.sp+n.middle.sp)){
        gammA[t,i,(k-(n.basal.sp))] <- adjacency.mat[i,k]*w[3]/(sqrt(2*sa[i]+2*sa[k]+w[3]^2))*exp(-(z[t,i]-z[t,k])^2/(2*sa[i]+2*sa[k]+w[3]^2))
        betaA[t,i,(k-(n.basal.sp))] <- -adjacency.mat[i,k]*w[3]*2*sa[i]*(z[t,i]-z[t,k])/((2*sa[i]+2*sa[k]+w[3]^2)^1.5)*exp(-(z[t,i]-z[t,k])^2/(2*sa[i]+2*sa[k]+w[3]^2))
      }
    }
    
    
    
    
    t.gammA[t,,]<- -t(gammA[t,,])
    t.betaA[t,,]<- -t(betaA[t,,])
    
    
    
    
    for( i in (n.basal.sp+1): (n.basal.sp+n.middle.sp)){
      for (k in (n.basal.sp+n.middle.sp+1) :(n.basal.sp+n.consumers+n.middle.sp) )
      {
        M.gammA[t,(i-(n.middle.sp+1)),(k-(n.basal.sp+n.middle.sp))] <- adjacency.mat[i,k]*gamma/(sqrt(2*sa[i]+2*sa[k]+gamma^2))*exp(-(z[t,i]-z[t,k])^2/(2*sa[i]+2*sa[k]+gamma^2))
        M.betaA[t,(i-(n.middle.sp+1)),(k-(n.basal.sp+n.middle.sp))] <- -adjacency.mat[i,k]*gamma*2*sa[i]*(z[t,i]-z[t,k])/((2*sa[i]+2*sa[k]+gamma^2)^1.5)*exp(-(z[t,i]-z[t,k])^2/(2*sa[i]+2*sa[k]+gamma^2))
        
      }
    }
    
    dM.gammA[t,,]<- -t(M.gammA[t,,])
    dM.betaA[t,,]<- -t(M.betaA[t,,])
    
    
    
    for (j in 1: n.basal.sp){ 
      
      if (t < 500) {
        R[t+1,j]<- R[t,j] + R[t,j]*(det_biB[t,j]- sum(alpha[j,]*R[t,]) + sum(1*gammA[t,j,]*M[t,]))*dt +
          dt*rnorm(1,mean = 0,sd=0.1)*R[t,j]^(a/2) 
        
        # 
        rg[t,j]<- R[t,j]*(det_biB[t,j]- sum(alpha[j,]*R[t,]) + sum(1*gammA[t,j,]*M[t,]))*dt 
        
        
        gA[t,j]  <-  (sqrt(sa[1,j]))/(sqrt(2*pi))*exp(-(thetaB+z[t,j])^2/(2*sa[1,j]))*(1- exp(2*thetaB*z[t,j]/sa[1,j])) # this is the growth rate of the trait z
        
        z[t+1,j]<- z[t,j] + sa[1,j]*h*(gA[t,j]  + sum(betaA[t,j,]*M[t,]))+ dt*rnorm(1,0,0.005)
        
      }
      
      else if (t >= 500){
        if ( forcing.n.basal.sp == 0 ) {
          ff[t,]<- 0
        }
        else if ( forcing.n.basal.sp != 0 ){
          ff[t,forcing.n.basal.sp] <-env.change[t-499]
        }
        R[t+1,j]<- R[t,j] + R[t,j]*(det_biB[t,j]- sum(alpha[j,]*R[t,]) + sum(1*gammA[t,j,]*M[t,]))*dt +
          dt*rnorm(1,mean = 0,sd=0.05)*R[t,j]^(a/2) -ff[t,j]*R[t,j]^2/(1+R[t,j]^2)*dt 
        
    
        rg[t,j]<- R[t,j]*(det_biB[t,j]- sum(alpha[j,]*R[t,]) + sum(1*gammA[t,j,]*M[t,]))*dt 
        
        
        gA[t,j]  <-  (sqrt(sa[1,j]))/(sqrt(2*pi))*exp(-(thetaB+z[t,j])^2/(2*sa[1,j]))*(1- exp(2*thetaB*z[t,j]/sa[1,j])) # this is the growth rate of the trait z
        z[t+1,j]<-z[t,j] + sa[1,j]*h*(gA[t,j] + sum(betaA[t,j,]*M[t,]))+ dt*rnorm(1,0,0.005)
      }
      if (R[t+1,j] < 1e-5){ R[t+1,j] <- 0}
    } 
    
    
    
    
    for ( i in 1: middle.sp){
      if(t < 500 ){
        M[t+1,i]<- M[t,i]+ M[t,i]*( det_biM[t,i]- sum(m.alpha[i,]*M[t,])+ sum(1*t.gammA[t,i,]*R[t,]) + sum( M.gammA[t,i,]*C[t,]))*dt+
          dt*rnorm(1,mean = 0,sd=0.1)*M[t,i]^(a/2) 
        mg[t,i]<-M[t,i]*( det_biM[t,i]- sum(m.alpha[i,]*M[t,])+ sum(1*t.gammA[t,i,]*R[t,]) + sum( M.gammA[t,i,]*C[t,]))*dt
        
        z[t+1,(n.basal.sp+i)]<-z[t,(n.basal.sp+i)] + 
          sa[1,(n.basal.sp+i)]*h*( sum(t.betaA[t,i,]*R[t,]) + sum(M.betaA[t,i,]*C[t,])) + dt*rnorm(1,0,0.005)
      }
      
      else if ( t >= 500){
        if ( forcing.middle.sp == 0 ){
          fc[t,] <- 0
          
        }else if ( forcing.middle.sp != 0 ){
          
          fc[t,forcing.middle.sp] <- env.change[t-499]
          
        }
        
        M[t+1,i]<- M[t,i]+ M[t,i]*( det_biM[t,i]- sum(m.alpha[i,]*M[t,])+ sum(1*t.gammA[t,i,]*R[t,]) + sum( M.gammA[t,i,]*C[t,]))*dt+
          dt*rnorm(1,mean = 0,sd=0.1)*M[t,i]^(a/2) - fc[t,i]*dt*M[t,i]^2/(1+M[t,i]^2)
        
        mg[t,i]<-M[t,i]*( det_biM[t,i]- sum(m.alpha[i,]*M[t,])+ sum(1*t.gammA[t,i,]*R[t,]) + sum( M.gammA[t,i,]*C[t,]))*dt
        
        
        z[t+1,(n.basal.sp+i)]<-z[t,(n.basal.sp+i)] + 
          sa[1,(n.basal.sp+i)]*h*( sum(t.betaA[t,i,]*R[t,])+ sum(M.betaA[t,i,]*C[t,])) + dt*rnorm(1,0,0.005)
        
      }
      if (M[t+1,i] < 1e-5){ M[t+1,i] <- 0}
      
    }
    
    for ( i in 1: n.consumers){ 
      if (t < 500){
        C[t+1,i]<- C[t,i]+ C[t,i]*( det_bip[t,i] - 0.005*C[t ,i]  + sum(1*dM.gammA[t,i,]*M[t,]) )*dt+
          dt*rnorm(1,mean = 0,sd=0.05)*C[t,i]^(a/2)
        
        cg[t,i]<- C[t,i]*( det_bip[t,i] - 0.005*C[t ,i]  + sum(1*dM.gammA[t,i,]*M[t,]) )*dt
        
        
        z[t+1,(n.basal.sp+middle.sp+i)]<- z[t,(n.basal.sp+middle.sp+i)] + 
          sa[1,(n.basal.sp+middle.sp+i)]*h*(  sum(dM.betaA[t,i,]*M[t,]))+ dt*rnorm(1,0,0.005)
        
      }
      
      else if (t >= 500){
        if (forcing.consumer == 0 )
        {c[t+1, ] <- 0 }
        else if (forcing.consumer != 0 ){
          c[t+1,forcing.consumer] <-env.change[t-499]
        }
        C[t+1,i]<- C[t,i]+ C[t,i]*( det_bip[t,i] - 0.005*C[t ,i] + sum( 1*dM.gammA[t,i,]*M[t,]) )*dt+
          dt*rnorm(1,mean = 0,sd=0.1)*C[t,i]^(a/2) - c[t,i]*dt*C[t,i]^2/(1+C[t,i]^2)
        
        cg[t,i]<- C[t,i]*( det_bip[t,i] - 0.005*C[t ,i]  + sum(1*dM.gammA[t,i,]*M[t,]) )*dt
        
        z[t+1,(n.basal.sp+middle.sp+i)]<- z[t,(n.basal.sp+middle.sp+i)] + 
          sa[1,(n.basal.sp+middle.sp+i)]*h*(  sum(dM.betaA[t,i,]*M[t,])) + dt*rnorm(1,0,0.01)
      }
      
      if (C[t+1,i] < 1e-5){ C[t+1,i] <- 0}
      
    } 
    
    
    
    
    
    
    
  }
  
  
  output= list(Resource = R[1:start.time,], Consumer = C[1:start.time,], Middle.sp = M[1:start.time,], trait= z[1:start.time,],
               M.gammA=M.gammA[1:start.time,,], gammA=gammA[1:start.time,,] ,sa=sa,env=env.change,
               rg=rg, cg=cg, mg=mg)
  return(output)
}


Mcommunity = function(iter, tt, ...){
  set.seed(as.numeric(Sys.time())-Sys.getpid()) 
  init = tt
  replicate =foodweb.model.ews(tt=init, ...)
  replicate$start = init
  replicate$iter = iter
  return(replicate)
}
#

