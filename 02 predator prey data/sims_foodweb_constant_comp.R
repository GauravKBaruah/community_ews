# R script for species collapse simulations for predator prey food-web module

rm(list = ls())


library(dplyr)
library(mgcv)
library(vegan)
library(tidyr)
source("01_fweb_functions_analyses.R", echo=F)
source("food_web_constant_Competition.R", echo=F)


load('FOOD_WEB_CONSTANT_COMPETITION.RData')
parameters 


# incidence matrix for the predator-prey module
g<-adj.matrix<-matrix(c(-1,0,0,-1,-1,0,
                        0,-1,0,-1,-1,0,
                        0,0,-1, -1,-1,0,
                        1,1,1,0,0,-1,
                        1,1,1,0,0,-1,
                        0,0,0,1,1,0),  nrow=6, ncol=6, byrow = T)


parameters$R<-c(1,  1,  1, 
                -0.02500000, -0.010000, -0.008)
parameters$gamma<-0.005 #top predator and basal resource interaction width
parameters$w[3]<-0.005 #basal-resource interaction width
parameters$initialmu<-c(0.77509908,  0.09888011, -0.75051236, -0.74836734,  0.5858192, -0.83202347)


parameters$A_matrix<-matrix.competition(middle.sp = 2,basal.sp = 3)$Amatrix
parameters$P_matrix<-matrix.competition(middle.sp = 2,basal.sp = 3)$Pmatrix


parameters$Omega<- 0.05 #0.01, 0.03, 0.05
parameters$forcing.type<-'stochastic'
parameters$Variance <- runif(6, 0.007,0.02)
reps=1



#empty data frame for EWS
f.data<-expand.grid(forcing_species=c(1,2,3,4,5,6),
                       `random_seed`=4327+(1:50)*100) %>%
  as_tibble %>%
  mutate(`community_Ar1`=0,
         `community_SD`=0,
         `community_trait`=0,
         `species_Ar1`=0,
         `species_SD`=0,
         `species_trait`=0,
         `Dominant.eigenvalue`=0)


#empty data frame for epotential curves
potential.data<-expand.grid(forcing_species=c(rep(1,9001),rep(2,9001),rep(3,9001),
                                              rep(4,9001),rep(5,9001),rep(6,9001)),
                            `random_seed`=4327+(1:50)*100) %>% 
  as_tibble %>% 
  mutate(V=0,
         density=0)



#empty data frame for metrics from potential curves
potential.curve.estimates<-expand.grid(forcing_species=c(1,2,3,4,5,6),
                                       `random_seed`=4327+(1:50)*100) %>% 
  as_tibble %>% 
  mutate(slope=0,
         width.potential=0,
         potential.depth=0,
         diffsp=0)



#empty dataframe for metrics from wentzell potential function 
Wentzell_potential<-expand.grid(forcing_species=c(1,2,3,4,5,6),
                                       `random_seed`=4327+(1:50)*100) %>% 
  as_tibble %>% 
  mutate(slope=0,
         width.potential=0,
         potential.depth=0,
         diffsp=0)






parameters$forcing.type <- 'stochastic'

start.time <-5000




for(r in 1:nrow(f.data)){
        model.temp <-  lapply(1, Mcommunity, tt = start.time, 
                                                    parms = parameters, adjacency.mat = g,
                                                    forcing.species=f.data$forcing_species[r])
      
      
      model.temp[[1]]$Density<-cbind(model.temp[[1]]$Resource,model.temp[[1]]$Middle.sp, model.temp[[1]]$Consumer)
      model.temp[[1]]$r<-cbind(model.temp[[1]]$rg,model.temp[[1]]$mg, model.temp[[1]]$cg)

      xtip.stoch<-extent.decline.foodweb(model.temp[[1]],species =f.data$forcing_species[r] ) 

      
     Uab<- Congugate_gradient_method(r = model.temp[[1]]$r,N =  model.temp[[1]]$Density,S = 6,
                                     species_index = f.data$forcing_species[r])
      
      #Species AR1
      f.data$species_Ar1[r]<-    genericEWS(model.temp[[1]]$Density[200:round(mean(xtip.stoch,na.rm=T),0),f.data$forcing_species[r]],
                                              winsize = 50,detrending = "no",bandwidth = 50)$KTauAcf[1]

      #species SD
      f.data$species_SD[r]<-    (genericEWS(model.temp[[1]]$Density[200:round(mean(xtip.stoch,na.rm=T),0),f.data$forcing_species[r]],
                                              winsize = 50,detrending = "no",bandwidth = 50)$KTauSD[1])
      #species trait
      f.data$species_trait[r] <- genericEWS(model.temp[[1]]$trait[200:round(mean(xtip.stoch,na.rm=T),0),f.data$forcing_species[r]],
                                              winsize = 50,detrending = "no",bandwidth = 50)$Ktau.trait[1]

      #dominant eigenvalue
      f.data$Dominant.eigenvalue[r] <-delta.shift.eigenvalue(var.cov.mat(15,model.temp[[1]], 6,
                                                                           tipping.point = round(mean(xtip.stoch,na.rm=T),0)))

      #community SD
      biomass.comm.stochastic<-community.measures(data =model.temp,
                                                  tipping.point = round(mean(xtip.stoch,na.rm=T),0),
                                                  harvest.species =f.data$forcing_species[r]   )$total.abundance

      f.data$community_SD[r]<- genericEWS(biomass.comm.stochastic[,1],winsize = 50,
                                            detrending = "no",bandwidth = 50)$KTauSD[1]

      #community AR1
      f.data$community_Ar1[r]<-genericEWS(biomass.comm.stochastic[,1],winsize = 50,
                                            detrending = "no",bandwidth = 50)$KTauAcf[1]

      #community trait
      mean.comm.stochastic<-community.measures(data =model.temp,
                                               tipping.point = round(mean(xtip.stoch,na.rm=T),0),
                                               harvest.species = f.data$forcing_species[r])$community.traitmean

      f.data$community_trait[r]<-genericEWS(mean.comm.stochastic[,1],winsize = 50,
                                              detrending = "no",bandwidth = 50)$Ktau.trait[1]




      #potential data
      vv<-epotential(model.temp,forcing.species = f.data$forcing_species[r])
      
   
      
      potential.data$V[(r*9001+1-9001):(r*9001)]<- vv$V.mean
      potential.data$density[(r*9001+1-9001):(r*9001)]<-vv$N1
      
      
      Wentzell_potential$slope[r]<-Uab$slope.sp
      Wentzell_potential$width.potential[r]<-Uab$width.potential.sp
      Wentzell_potential$potential.depth[r]<-Uab$pot.depth.sp
      
            
      potential.curve.estimates$slope[r]<- vv$slope.sp
      potential.curve.estimates$width.potential[r] <- vv$width.potential.sp
      potential.curve.estimates$potential.depth[r] <- vv$pot.depth.sp
      potential.curve.estimates$diffsp[r]<- vv$diffsp
      print(r)
      
      
  }

#save(Wentzell_potential, file="Wentzell_potential_estimates.RData")
#save(f.data,file="foodweb_data_results.RData")
#save(potential.data,file="foodweb_data_potential_results_v2.RData")
#save(potential.curve.estimates,file="foodweb_data_potential_curve_ests_v2.RData")
