
library(entropy)
library(pracma)
#library(igraph)
#library(bipartite)



set.seed(1234)



#mutualism community function
#eco-evolutionary dynamics of abundance and trait dynamics of two guilds of species
#returns: abundance of species over time, mean traits over time,variance over time
#just for sake of clarity, the two guilds of species in the mutualistic community are assummed to be plants and animals
Mutualistic.model.ews<-function(time=start.time,matrix,parameters=parameters,forcing.species=0)
{
  
  time<-start.time
  Aspecies<- dim(matrix)[2]
  Plantspecies<- dim(matrix)[1]
  
  if(forcing.species==0){ 
    forcing.plant.no <- 0
    forcing.animal.no<-0
    } else if(forcing.species==5 ){ 
    forcing.animal.no = (Aspecies+Plantspecies)-forcing.species 
    forcing.plant.no = 0 
    } else if(forcing.species==6){
      forcing.animal.no <- 2
      forcing.plant.no<-0
    }  else if (forcing.species ==1 || forcing.species == 2 ||forcing.species==3 || forcing.species==4){
    forcing.plant.no = forcing.species
    forcing.animal.no=0
  }
  
  
degree.animals<-degree.plants<-numeric()
  
  
  for(i in 1:Plantspecies){
    degree.plants[i]<-sum(matrix[i,])}
  for(j in 1:Aspecies){
    degree.animals[j]<-sum(matrix[,j])
  }
  
  
  
  a=2
  alphaP<- betaP <- array(dim=c( Plantspecies,Plantspecies))
  alphaA<- betaA <- array(dim=c(Aspecies,Aspecies))
  gamma_A<- two_betaA<-array(dim=c(time, Aspecies,Plantspecies))
  gamma_P<-two_betaP<-array(dim=c(time, Plantspecies,Aspecies))
  Np<-mup <-rp<- det_biP <- Ap <- gP <- MutP <- xetaP<-b_barP <-Kp<- fp<-array(dim=c(time,Plantspecies))
  Na<-muA<- ra<-det_biA<- Aa <- gA<- MutA<- xetaA<- b_barA<-Ka<-fa<-array(dim=c(time, Aspecies))
  ca<-array(dim=c(Aspecies,Plantspecies))
  cp<-array(dim=c(Plantspecies,Aspecies))
  sp<-hp<-array(dim=c(1, Plantspecies))
  sa<-ha<-array(dim=c(1, Aspecies))
  fp<-matrix(0,time,Plantspecies)
  fa<-matrix(0, time, Aspecies)
  Np[1,]<-parameters$InitialAB[1:Plantspecies]
  Na[1,]<- parameters$InitialAB[(Plantspecies+1) : (Plantspecies+Aspecies)]# initital species density
  mup[1,]<-parameters$initialmu[1:Plantspecies]
  muA[1,]<-parameters$initialmu[(Plantspecies+1) : (Plantspecies+Aspecies)]
  hp[1,]<- runif(Plantspecies,0.4,0.4)
  ha[1,]<- runif(Aspecies,0.4,0.4)
  
  sp[1,]<-parameters$Variance[1:Plantspecies]  # Intravar  #intraspecific trait variation 
  sa[1,]<-parameters$Variance[(Plantspecies+1) : (Plantspecies+Aspecies)]# Intravar  #intraspecific trait variation 
  dt<-parameters$dt
  w<-parameters$gamma
  omegaA<-parameters$Omega[1]
  omegaP<-parameters$Omega[2]
  thetaA<-parameters$Theta
  thetaP<-parameters$Theta  
  
  env.change.plants<-env.function(name = parameters$forcing.type)
  env.change.animals<-env.function(name = parameters$forcing.type)
  
  
  for (t in 1:(time-1)) {
    
    
    # growth function of the model for Plant species and animal species
    for ( j in 1: (Plantspecies)){
      
    if(j==2){
      det_biP[t,j]<- 0.5*(erf((thetaP-mup[t,j])/(sqrt(2*sp[1,j])))+erf((thetaP+mup[t,j])/(sqrt(2*sp[1,j]))))
    }else{
      det_biP[t,j]<- 0.5*(erf((thetaP-mup[t,j])/(sqrt(2*sp[1,j])))+erf((thetaP+mup[t,j])/(sqrt(2*sp[1,j]))))
      
    }
     
      
      
    }
    
    for ( j in 1: (Aspecies)){
     if(j==1){
        det_biA[t,j] <- 0.5*(erf((thetaA-muA[t,j])/(sqrt(2*sa[1,j])))+erf((thetaA+muA[t,j])/(sqrt(2*sa[1,j])))) 
     }else{
       det_biA[t,j] <- 0.5*(erf((thetaA-muA[t,j])/(sqrt(2*sa[1,j])))+erf((thetaA+muA[t,j])/(sqrt(2*sa[1,j])))) 
     }
    }
    
      #trait independent competition matrix
        alphaP <- parameters$Pmatrix
        alphaA <-   parameters$Amatrix
       
    #Mutualism term for animal and plant species
    for ( i in 1: Aspecies ) {
      for (k in 1: Plantspecies){
        
        
          ca[i,k]<-matrix[k,i]# /degree.animals[i]
          gamma_A[t,i,k] <-1*w/sqrt(2*sp[1,k]+2*sa[1,i]+w^2)*exp(-(muA[t,i]-mup[t,k])^2/(2*sp[1,k]+2*sa[1,i]+w^2))
          two_betaA[t,i,k] <-  -1*w*sa[1,i]*(muA[t,i]-mup[t,k])/(2*sp[1,k]+2*sa[1,i]+w^2)^1.5*exp(-(muA[t,i]-mup[t,k])^2/(2*sp[1,k]+2*sa[1,i]+w^2))
          
       
      
           }
      }
    
    
    # Mutualism term for plant and animal species  
    for ( i in 1: Plantspecies ) {
      for (k in 1: Aspecies){
        
        
        cp[i,k]<- matrix[i,k] #/degree.plants[i]
        gamma_P[t,i,k] <- 1*w/sqrt(2*sp[1,i]+2*sa[1,k]+w^2)*exp(-(mup[t,i]-muA[t,k])^2/(2*sp[1,i]+2*sa[1,k]+w^2)) 
        two_betaP[t,i,k] <- -1*w*sp[1,i]*(mup[t,i]-muA[t,k])/(2*sp[1,i]+2*sa[1,k]+w^2)^1.5*exp(-(mup[t,i]-muA[t,k])^2/(2*sp[1,i]+2*sa[1,k]+w^2)) 
        
      }
    }
    
    
    for ( i in 1: Aspecies){
      
      if (t<1000){ 
        gA[t,i]  <-  (sqrt(sa[1,i]))/(sqrt(2*pi))*exp(-(thetaA+muA[t,i])^2/(2*sa[1,i]))*(1- exp(2*thetaA*muA[t,i]/sa[1,i])) # this is the growth rate of the trait z
        xetaA[t,i] <- (sum(two_betaA[t,i,]*ca[i,]*Np[t,])) 
        muA[t+1,i]<- muA[t,i]  + sa[1,i]*ha[1,i]*(  gA[t,i] -   xetaA[t,i])*dt #+ dt*rnorm(1,mean = 0,sd=0.01)
        
        ra[t,i]<-(det_biA[t,i] -  sum(alphaA[i,]*Na[t,]) + sum(gamma_A[t,i,]*Np[t,]))*Na[t,i]*dt
        
        Na[t+1,i]<-Na[t,i]+ (det_biA[t,i] -  sum(alphaA[i,]*Na[t,]) + sum(gamma_A[t,i,]*Np[t,]))*Na[t,i]*dt+ 
          dt*rnorm(1,mean = 0,sd=0.1)*Na[t,i]^(a/2) 
        
      }
      else if (t>= 1000){
        if(forcing.animal.no == 0){
          
          fa[t+1,] <- fa[t,]
          
        }
        
        else if(forcing.animal.no != 0){
          
          fa[t+1,forcing.animal.no]<-env.change.animals[t-999]
          
          
        }
        
        Na[t+1,i]<-Na[t,i]+ (det_biA[t,i] -  sum(alphaA[i,]*Na[t,]) + sum(gamma_A[t,i,]*Np[t,]))*Na[t,i]*dt+
          dt*rnorm(1,mean = 0,sd=0.1)*Na[t,i]^(a/2) -fa[t,i]*Na[t,i]^2*dt/(1+Na[t,i]^2) 
        gA[t,i]  <-  (sqrt(sa[1,i]))/(sqrt(2*pi))*exp(-(thetaA+muA[t,i])^2/(2*sa[1,i]))*(1- exp(2*thetaA*muA[t,i]/sa[1,i])) # this is the growth rate of the trait z
        xetaA[t,i] <- (sum(two_betaA[t,i,]*ca[i,]*Np[t,])) 
        muA[t+1,i]<- muA[t,i]  + sa[1,i]*ha[1,i]*(  gA[t,i] +  xetaA[t,i])*dt+ dt*rnorm(1,mean = 0,sd=0.01)
        ra[t,i]<-(det_biA[t,i] -  sum(alphaA[i,]*Na[t,]) + sum(gamma_A[t,i,]*Np[t,]))*Na[t,i]*dt 
        
      }
      
      if (Na[t+1,i] < 1e-4) { Na[t+1,i]<- 0 } 
    }   
    #
    
    for ( i in 1: Plantspecies){
      
      
      if (t<1000){ 
        
        Np[t+1,i]<-Np[t,i]+ (det_biP[t,i] -   sum(alphaP[i,]*Np[t,])+sum(gamma_P[t,i,]*cp[i,]*Na[t,]))*Np[t,i]*dt +
          dt*rnorm(1,mean = 0,sd=0.1)*Np[t,i]^(a/2) # lotka volterra equation  1. det_bi is the growth rate and A[t,] is the pairwise competition coefficient matrix
        gP[t,i]  <-  (sqrt(sp[1,i]))/(sqrt(2*pi))*exp(-(thetaP+mup[t,i])^2/(2*sp[1,i]))*(1- exp(2*thetaA*mup[t,i]/sp[1,i])) # this is the growth rate of the trait z
        xetaP[t,i] <- (sum(two_betaP[t,i,]*cp[i,]*Na[t,])) 
        mup[t+1,i]<- mup[t,i] +  sp[1,i]*hp[1,i]*(gP[t,i] +  xetaP[t,i])*dt #+ dt*rnorm(1,mean = 0,sd=0.008)
        rp[t+1,i]<-(det_biP[t,i] -   sum(alphaP[i,]*Np[t,])+sum(gamma_P[t,i,]*cp[i,]*Na[t,]))*Np[t,i]*dt
        
      }
      else if (t>= 1000){
        if (forcing.plant.no == 0){
          fp[t+1,] <- 0 
        }
        else if(forcing.plant.no != 0){
          fp[t+1,forcing.plant.no] <-env.change.plants[t-999]
          #env.change[t-999]
        }
        
        Np[t+1,i]<-Np[t,i]+ (det_biP[t,i] -   sum(alphaP[i,]*Np[t,]) + sum(gamma_P[t,i,]*cp[i,]*Na[t,]))*Np[t,i]*dt +
          dt*rnorm(1,mean = 0,sd=0.1)*Np[t,i]^(a/2) - fp[t,i]*Np[t,i]^2*dt/(1+Np[t,i]^2)
        gP[t,i]  <-  (sqrt(sp[1,i]))/(sqrt(2*pi))*exp(-(thetaP+mup[t,i])^2/(2*sp[1,i]))*(1- exp(2*thetaA*mup[t,i]/sp[1,i])) # this is the growth rate of the trait z
        xetaP[t,i] <- (sum(two_betaP[t,i,]*cp[i,]*Na[t,])) 
        mup[t+1,i]<- mup[t,i] +  sp[1,i]*hp[1,i]*(gP[t,i]  +  xetaP[t,i])*dt # + dt*rnorm(1,mean = 0,sd=0.008)
        rp[t+1,i]<- (det_biP[t,i] -   sum(alphaP[i,]*Np[t,])+sum(gamma_P[t,i,]*cp[i,]*Na[t,]))*Np[t,i]*dt
        
      }
      
      
      if (Np[t+1,i] < 1e-4) { Np[t+1,i] <- 0 }
      
    }
  }
  
  output= list(Plants = Np[500:time,],Animals=Na[500:time,], Plant.trait = mup[500:time,], 
               Animal.trait=muA[500:time,], varP=sp[1,] , varA= sa[1,],fp=fp[500:time,],fa=fa[500:time,],
               Gamma=gamma_P[500:time,,],env.change=env.change.animals,ra=ra, rp=rp)
  return(output)
}

Mcommunity = function(iter, time, ...){
  set.seed(rnorm(1,as.numeric(Sys.time())-Sys.getpid(),10000)) 
  init = time
  replicate =try(Mutualistic.model.ews(time=init, ...))
  replicate$start = init
  replicate$iter = iter
  return(replicate)
}


