rm(list=ls())
#R script that simulates all the data of species collapse for the mutualism community



library(dplyr)


source("~/Dropbox/Zurich PhD Research/4_Chapter_4/Data archive/01 mutualism data/05_Mutualism_constant_compvalues.R")

load("~/Dropbox/Zurich PhD Research/4_Chapter_4/Data archive/01 mutualism data/DATA_MUTUALISM_CONSTANT_COMPETITION.RData")
source("~/Dropbox/Zurich PhD Research/4_Chapter_4/Data archive/01 mutualism data/01_functions.R", echo=F)





g<-adj.matrix<-matrix(c(1,1,
                        1,1,
                        1,1,
                        1,1), nrow=4, ncol=2, byrow = T)

start.time<-5000
#parameters for initialisation

parameters$forcing.type<-'stochastic'
parameters$Variance <- runif(6,0.005,0.06)
parameters$gamma<-0.1#0.01, 0.03, 0.05, 0.1
parameters$initialmu<- c(-0.5,0.6,0.12,-0.54,0.6,-0.5)
reps=1



#empty dataframes for EWS metrics
mut.data<-expand.grid(forcing_species=c(1,2,3,4,5,6),
                       `random_seed`=4327+(1:50)*100) %>%
  as_tibble %>%
  mutate(`community_Ar1`=0,
         `community_SD`=0,
         `community_trait`=0,
         `species_Ar1`=0,
         `species_SD`=0,
         `species_trait`=0,
         `Dominant.eigenvalue`=0)


#empty data frames for effective potential curve
 potential.data<-expand.grid(forcing_species=c(rep(1,501),rep(2,501),rep(3,501),
                                              rep(4,501),rep(5,501),rep(6,501)),
                            `random_seed`=4327+(1:50)*100) %>% 
  as_tibble %>% 
  mutate(V=0,
         density=0)

# empty data frame for metrics from effective potential curves                      
potential.curve.estimates<-expand.grid(forcing_species=c(1,2,3,4,5,6),
                            `random_seed`=4327+(1:50)*100) %>% 
  as_tibble %>% 
  mutate(slope=0,
         width.potential=0,
         potential.depth=0,
         diffsp=0)



# empty data frame for qusipotential wentzell friedlin potential
Wentzell_potential<-expand.grid(forcing_species=c(1,2,3,4,5,6),
                                `random_seed`=4327+(1:50)*100) %>% 
  as_tibble %>% 
  mutate(slope=0,
         width.potential=0,
         potential.depth=0,
         diffsp=0)


# simulations
for(r in 272:nrow(mut.data)){

  
  
  model.temp<-lapply(1, Mcommunity,time=start.time,matrix=g,
                              parameters=parameters,forcing.species=mut.data$forcing_species[r])
  model.temp[[1]]$Density<-cbind(model.temp[[1]]$Plants,model.temp[[1]]$Animals)
  model.temp[[1]]$trait<-cbind(model.temp[[1]]$Plant.trait, model.temp[[1]]$Animal.trait)
  model.temp[[1]]$r<-cbind(model.temp[[1]]$rp, model.temp[[1]]$ra)
  
 # ts.plot(model.temp[[1]]$Density,col=c("black","red","green","pink","orange","violet"))
  
  Uab<- Congugate_gradient_method(r = model.temp[[1]]$r,N =  model.temp[[1]]$Density,S = 6,
                                  species_index = mut.data$forcing_species[r])
  
  #plot(Uab$V) 
  xtip.stoch<-extent.decline.mutualism(model.temp[[1]],species =mut.data$forcing_species[r] ) 
  
  #Species AR1
  mut.data$species_Ar1[r]<-    genericEWS(model.temp[[1]]$Density[200:round(mean(xtip.stoch,na.rm=T),0),mut.data$forcing_species[r]],
                                           winsize = 50,detrending = "no",bandwidth = 50)$KTauAcf[1]

  #species SD
  mut.data$species_SD[r]<-    (genericEWS(model.temp[[1]]$Density[200:round(mean(xtip.stoch,na.rm=T),0),mut.data$forcing_species[r]],
                                           winsize = 50,detrending = "no",bandwidth = 50)$KTauSD[1])
  #species trait
  mut.data$species_trait[r] <- genericEWS(model.temp[[1]]$trait[200:round(mean(xtip.stoch,na.rm=T),0),mut.data$forcing_species[r]],
                                           winsize = 50,detrending = "no",bandwidth = 50)$Ktau.trait[1]

  #dominant eigenvalue
  mut.data$Dominant.eigenvalue[r] <-delta.shift.eigenvalue(var.cov.mat(15,model.temp[[1]], 6,
                                                                        tipping.point = round(mean(xtip.stoch,na.rm=T),0)))

  #community SD
  biomass.comm.stochastic<-community.measures(data =model.temp,
                                              tipping.point = round(mean(xtip.stoch,na.rm=T),0),
                                              harvest.species =mut.data$forcing_species[r]   )$total.abundance

  mut.data$community_SD[r]<- genericEWS(biomass.comm.stochastic[,1],winsize = 50,
                                         detrending = "no",bandwidth = 50)$KTauSD[1]

  #community AR1
  mut.data$community_Ar1[r]<-genericEWS(biomass.comm.stochastic[,1],winsize = 50,
                                         detrending = "no",bandwidth = 50)$KTauAcf[1]

  #community trait
  mean.comm.stochastic<-community.measures(data =model.temp,
                                           tipping.point = round(mean(xtip.stoch,na.rm=T),0),
                                           harvest.species = mut.data$forcing_species[r])$community.traitmean

  mut.data$community_trait[r]<-genericEWS(mean.comm.stochastic[,1],winsize = 50,
                                           detrending = "no",bandwidth = 50)$Ktau.trait[1]




  #potential data
 
  
  vv<-grates(model.temp,forcing.species=mut.data$forcing_species[r])
  potential.data$V[(r*501+1-501):(r*501)]<- vv$V.mean
  potential.data$density[(r*501+1-501):(r*501)]<-vv$N1
  
  potential.curve.estimates$slope[r]<- vv$slope.sp
  potential.curve.estimates$width.potential[r] <- vv$width.potential.sp
  potential.curve.estimates$potential.depth[r] <- vv$pot.depth.sp
  potential.curve.estimates$diffsp[r]<- vv$diffsp
  
  Wentzell_potential$width.potential[r]<-Uab$width.potential.sp
  Wentzell_potential$slope[r]<-Uab$slope.sp
  Wentzell_potential$potential.depth[r]<-Uab$pot.depth.sp
  
  #vv<-grates(model.temp,forcing.species=mut.data$forcing_species[r])$diffsp
  print(r)
  
  
}







#save(mut.data,file="Mutualism_data_results_1.RData")
#save(potential.data,file="Mutualism_data_potential_results_v3.RData")
#save(potential.curve.estimates,file="Mutualism_data_potential_curve_ests_v3.RData")
