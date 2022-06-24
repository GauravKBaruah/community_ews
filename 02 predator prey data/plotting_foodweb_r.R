#R script for plotting predator-prey food-web results

rm(list=ls())



library(dplyr)
library(ggplot2)
library(gridBase)
library(gridExtra)
library(earlywarnings)
library(RColorBrewer)
library(wesanderson)
library(tidyr)
library(grid)
library(gridExtra)
library(cowplot)
load("foodweb_data_results.RData")
load("foodweb_data_potential_curve_ests_v2.RData")
load("foodweb_data_potential_results_v2.RData")

best_color_paletter<- c(wes_palettes$Darjeeling1, wes_palettes$Royal2,wes_palettes$Rushmore1)

wdata<-gather(f.data, metric, value, community_Ar1:Dominant.eigenvalue, factor_key=TRUE)





wdata$metric  <- factor(wdata$metric, levels =c("Dominant.eigenvalue","community_Ar1", 
                                                "community_SD",
                                                "community_trait", "species_Ar1",
                                                "species_SD", "species_trait"))

(ews_est<-ggplot(data= wdata, aes(x = value,y = factor(forcing_species), color= metric )) +
  geom_jitter(aes( color = metric),alpha=0.3,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.4),
              size = 0.5
  ) +
  geom_vline(xintercept=0, linetype ="dashed")+
  stat_summary(
    aes(color = metric),
    fun.data="mean_sdl",  fun.args = list(mult=1),
    geom = "pointrange",  size = 0.6,
    position = position_dodge(0.5)
  )+ggtitle("B")+
    theme_cowplot()+xlab("Kendall's tau")+ylab("Species forced to collapse")+
  scale_color_manual(values=best_color_paletter,
                     breaks = c( "species_trait","species_SD","species_Ar1","community_trait",
                                 "community_SD","community_Ar1","Dominant.eigenvalue"))+
  theme(plot.title = element_text(size = 15,   face = "bold"),text = element_text(size = 15 ),
        axis.title = element_text(face="bold"),legend.position = "right"))




  
species1<-potential.data %>% filter(forcing_species==1 )
species2<-potential.data %>% filter(forcing_species==2 )
species3<-potential.data %>% filter(forcing_species==3 )
species4<-potential.data %>% filter(forcing_species==4 )
species5<-potential.data %>% filter(forcing_species==5 )
species6<-potential.data %>% filter(forcing_species==6 )

d1<-v1<-d2<-v2<-d3<-v3<-d4<-v4<-d5<-v5<-d6<-v6<-matrix(NA,nrow=9001,ncol=50)
for(r in 1:50){
  

  d1[,r]<-species1$density[(r*9001+1-9001) : (r*9001)]
  v1[,r]<-species1$V[(r*9001+1-9001) : (r*9001)]
  
  d2[,r]<-species2$density[(r*9001+1-9001) : (r*9001)]
  v2[,r]<-species2$V[(r*9001+1-9001) : (r*9001)]
  
  d3[,r]<-species3$density[(r*9001+1-9001) : (r*9001)]
  v3[,r]<-species3$V[(r*9001+1-9001) : (r*9001)]
  
  d4[,r]<-species4$density[(r*9001+1-9001) : (r*9001)]
  v4[,r]<-species4$V[(r*9001+1-9001) : (r*9001)]
  
  d5[,r]<-species5$density[(r*9001+1-9001) : (r*9001)]
  v5[,r]<-species5$V[(r*9001+1-9001) : (r*9001)]
  
  d6[,r]<-species6$density[(r*9001+1-9001) : (r*9001)]
  v6[,r]<-species6$V[(r*9001+1-9001) : (r*9001)]
  
}

color.pal<-c("firebrick", "grey","darkgreen","darkblue","mediumpurple2","lightsalmon3")
potential.03<-data.frame( Potential=c(rowMeans(v1),rowMeans(v2),rowMeans(v3),rowMeans(v4),
                                      rowMeans(v5),rowMeans(v6)),
                          Density=c(rowMeans(d1),rowMeans(d2),rowMeans(d3),rowMeans(d4),
                                    rowMeans(d5),rowMeans(d6)),
                          Species = factor(c(rep("Species 1",each=9001),
                                             rep("Species 2", each=9001),
                                             rep("Species 3", each=9001),
                                             rep("Species 4", each=9001),
                                             rep("Species 5", each=9001),
                                             rep("Species 6", each=9001) ))) 



ggplot(data= potential.03, aes(y = Potential,x =Density , color = Species )) +geom_line(size=.5)+
  ylab("Effective potential V(N)")+xlab("Population density")+ggtitle("A")+
  theme_classic()+
  theme(plot.title = element_text(size = 10,   face = "bold"),
        text = element_text(size = 10 ),
        axis.title = element_text(face="bold"),
        legend.position = "right") +ylim(c(-0.5,4))+
  scale_color_manual(values=color.pal )+ labs(color="Species")



potential.curve<-ggplot(data= potential.03, aes(y = Potential,x =Density , color = Species )) +geom_line(size=1)+
  ylab("Effective potential V(N)")+xlab("Population density")+ggtitle("A")+
  theme_classic()+
  theme(plot.title = element_text(size = 10,   face = "bold"),
        text = element_text(size = 10 ),
        axis.title = element_text(face="bold"),
        legend.position = "right")+ylim(c(-1,10))+
  scale_color_manual(values=color.pal )+ labs(color="Species")


dd<-data.frame(images=list.files("~/Dropbox/Zurich PhD Research/4_Chapter_4/MultispeciesEWS/0_Codes & rmd/Potential plot figs/", full.names = T), stringsAsFactors = F)
dd$names = gsub("[a-zA-Z]|[[:punct:]]","",dd$images)
dd$values = sample(1:50, size=nrow(dd))
img = readPNG(dd$images[1])
g =  rasterGrob(img, interpolate=TRUE)


#g = list()

min_v<-potential.03 %>% dplyr::group_by(Species) %>% 
  dplyr::summarise(Min_v= min(Potential))

N1<-filter(potential.03, Species=="Species 1")$Density
N2<-filter(potential.03, Species=="Species 2")$Density
N3<-filter(potential.03, Species=="Species 3")$Density
N4<-filter(potential.03, Species=="Species 4")$Density
N5<-filter(potential.03, Species=="Species 5")$Density
N6<-filter(potential.03, Species=="Species 6")$Density

V1<-filter(potential.03, Species=="Species 1")$Potential
V2<-filter(potential.03, Species=="Species 2")$Potential
V3<-filter(potential.03, Species=="Species 3")$Potential
V4<-filter(potential.03, Species=="Species 4")$Potential
V5<-filter(potential.03, Species=="Species 5")$Potential
V6<-filter(potential.03, Species=="Species 6")$Potential

KK<-c(N1[which(V1== min(V1))],N2[which(V2 == min(V2))],
      N3[which(V3 == min(V3))],
      N4[which(V4 == min(V4))],
      N5[which(V5 == min(V5))],
      N6[which(V6 == min(V6))])
# this assembles the vector of carrying capacity or the point at which the population has the lowest effective potential value.


VV<-c(min(V1),min(V2),min(V3),min(V4),min(V5),min(V6)) # this vector assembles the lowest value of effective potential

# this for low below manages to put the .png pics of the ball at the lowest point in the effective potential curve
for(i in 1:(nrow(dd))){
  img = readPNG(dd$images[i])
  g[[i]] =  rasterGrob(img, interpolate=TRUE)
  
  potential.curve = potential.curve +
    annotation_custom(grob=g[[i]], xmin=KK[i]-0.45, xmax=KK[i]+0.45, ymin=VV[i], ymax=VV[i]+0.35)
}

potential.curve



# estimates of metrics from the potetnail curve
potential.curve.estimates$slope[is.nan(potential.curve.estimates$slope)] <- 0.000001
slope_est<-ggplot(data= potential.curve.estimates, 
                  aes(y = slope,x = factor(forcing_species), color= factor(forcing_species ))) +
  geom_boxplot(alpha=0,position= position_dodge(1))+
  geom_point(pch = 16,alpha=0.15, position = position_jitterdodge(1))+
  scale_color_manual(values = color.pal)+
  theme_cowplot()+xlab("Species")+ggtitle("D")+
  theme(plot.title = element_text(size = 10,
                                  face = "bold"),text = element_text(size = 10 ),
        axis.title = element_text(face="bold"),legend.position = "none") + 
  labs(color="Species")+ylab("Potential Slope")


dept_est<-ggplot(data= potential.curve.estimates, 
                 aes(y = -potential.depth,x = factor(forcing_species), color= factor(forcing_species ))) +
  geom_boxplot(alpha=0,position= position_dodge(1))+
  geom_point(pch = 16,alpha=0.15, position = position_jitterdodge(1))+
  scale_color_manual(values = color.pal)+
  theme_cowplot()+xlab("Species")+ggtitle("E")+
  theme(plot.title = element_text(size = 10,
                                  face = "bold"),text = element_text(size = 10 ),
        axis.title = element_text(face="bold"),legend.position = "none") + 
  labs(color="Species")+ylab("Potential depth") 

width_est<-ggplot(data= potential.curve.estimates, 
                  aes(y = width.potential,x = factor(forcing_species), color= factor(forcing_species ))) +
  geom_boxplot(alpha=0,position= position_dodge(1))+
  geom_point(pch = 16,alpha=0.15, position = position_jitterdodge(1))+
  scale_color_manual(values = color.pal)+
  theme_cowplot()+xlab("Species")+ggtitle("F")+
  theme(plot.title = element_text(size = 10,
                                  face = "bold"),text = element_text(size = 10 ),
        axis.title = element_text(face="bold"),legend.position = "none") + 
  labs(color="Species")+ylab("Potential width") 


dif_est<-ggplot(data= potential.curve.estimates, 
                  aes(y = diffsp^2,x = factor(forcing_species), color= factor(forcing_species ))) +
  geom_boxplot(alpha=0,position= position_dodge(1))+
  geom_point(pch = 16,alpha=0.15, position = position_jitterdodge(1))+
  scale_color_manual(values = color.pal)+
  theme_cowplot()+xlab("Species")+ggtitle("B")+
  theme(plot.title = element_text(size = 10,
                                  face = "bold"),text = element_text(size = 10 ),
        axis.title = element_text(face="bold"),legend.position = "none") + 
  labs(color="Species")+ylab("gains through predation/growth") 


source('~/Dropbox/Zurich PhD Research/1_Chapter_1/rolling_GAMs_methods.R', echo=F)

source("C:/Users/Baruah Gaurav/Dropbox/Zurich PhD Research/1_Chapter_1/rolling_GAMs_methods.R", echo=F)

lay_out(list(slope_est, 1, 1),
        list(dept_est,1,2),
        list(width_est,1,3) )



#wentzell potential estimates


slope_est_Wp<-ggplot(data= Wentzell_potential, 
                  aes(y = slope,x = factor(forcing_species), color= factor(forcing_species ))) +
  geom_boxplot(alpha=0,position= position_dodge(1))+
  geom_point(pch = 16,alpha=0.15,size=2, position = position_jitterdodge(1))+
  scale_color_manual(values = color.pal)+
  ylim(c(-0.1, 0.01))+
  theme_cowplot()+xlab("Species")+ggtitle("A")+
  theme(plot.title = element_text(size = 10,
                                  face = "bold"),text = element_text(size = 10 ),
        axis.title = element_text(face="bold"),legend.position = "none") + 
  labs(color="Species")+ylab("Potential Slope")


dept_est_wp<-ggplot(data= Wentzell_potential, 
                 aes(y = -potential.depth,x = factor(forcing_species), color= factor(forcing_species ))) +
  geom_boxplot(alpha=0,position= position_dodge(1))+
  geom_point(pch = 16,size=2,alpha=0.15, position = position_jitterdodge(1))+
  scale_color_manual(values = color.pal)+
  ylim(c(0.025,0))+
  theme_cowplot()+xlab("Species")+ggtitle("B")+
  theme(plot.title = element_text(size = 10,
                                  face = "bold"),text = element_text(size = 10 ),
        axis.title = element_text(face="bold"),legend.position = "none") + 
  labs(color="Species")+ylab("Potential depth") 

width_est_wp<-ggplot(data= Wentzell_potential, 
                  aes(y = abs(width.potential),x = factor(forcing_species), color= factor(forcing_species ))) +
  geom_boxplot(alpha=0,position= position_dodge(1))+
  geom_point(pch = 16,size=2,alpha=0.15, position = position_jitterdodge(1))+
  scale_color_manual(values = color.pal)+
  ylim(c(0,0.45))+
  theme_cowplot()+xlab("Species")+ggtitle("C")+
  theme(plot.title = element_text(size = 10,
                                  face = "bold"),text = element_text(size = 10 ),
        axis.title = element_text(face="bold"),legend.position = "none") + 
  labs(color="Species")+ylab("Potential width") 


lay_out(list(slope_est_Wp, 1, 1),
        list(dept_est_wp,1,2),
        list(width_est_wp,1,3) )

