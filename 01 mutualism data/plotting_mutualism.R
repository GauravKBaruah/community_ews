rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridBase)
library(gridExtra)
library(earlywarnings)
library(RColorBrewer)
library(wesanderson)
library(cowplot)

load("Mutualism_data_results_1.RData")
load("Mutualism_data_potential_results_v3.RData")
load("Mutualism_data_potential_curve_ests_v3.RData")

best_color_paletter<- c(wes_palettes$Darjeeling1, wes_palettes$Royal2,wes_palettes$Rushmore1)

wdata<-gather(mut.data, metric, value, community_Ar1:Dominant.eigenvalue, factor_key=TRUE)





wdata$metric  <- factor(wdata$metric, levels =c("Dominant.eigenvalue","community_Ar1", 
                                                "community_SD",
                                                "community_trait", "species_Ar1",
                                                "species_SD", "species_trait"))



(w1<-ggplot(data= wdata, aes(x = value,y = factor(forcing_species), color= metric )) +
  geom_jitter(aes(color = metric),alpha=0.3,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.4),
              size = 0.5
  ) +
  geom_vline(xintercept=0, linetype ="dashed")+
  stat_summary(
    aes(color = metric),
    fun.data="mean_sdl",  fun.args = list(mult=1),
    geom = "pointrange",  size = 0.6,
    position = position_dodge(0.5)
  )+ggtitle("A")+
  theme_cowplot()+xlab("Kendall's tau")+ylab("Species forced to collapse")+
  scale_color_manual(values=best_color_paletter,
                     breaks = c( "species_trait","species_SD","species_Ar1","community_trait",
                                 "community_SD","community_Ar1","Dominant.eigenvalue"))+
    theme(plot.title = element_text(size = 15,   face = "bold"),text = element_text(size = 15 ),
          axis.title = element_text(face="bold"),legend.position = "right") )




  
species1<-potential.data %>% filter(forcing_species==1 )
species2<-potential.data %>% filter(forcing_species==2 )
species3<-potential.data %>% filter(forcing_species==3 )
species4<-potential.data %>% filter(forcing_species==4 )
species5<-potential.data %>% filter(forcing_species==5 )
species6<-potential.data %>% filter(forcing_species==6 )

d1<-v1<-d2<-v2<-d3<-v3<-d4<-v4<-d5<-v5<-d6<-v6<-matrix(NA,nrow=501,ncol=50)
for(r in 1:50){
  

  d1[,r]<-species1$density[(r*501+1-501) : (r*501)]
  v1[,r]<-species1$V[(r*501+1-501) : (r*501)]
  
  d2[,r]<-species2$density[(r*501+1-501) : (r*501)]
  v2[,r]<-species2$V[(r*501+1-501) : (r*501)]
  
  d3[,r]<-species3$density[(r*501+1-501) : (r*501)]
  v3[,r]<-species3$V[(r*501+1-501) : (r*501)]
  
  d4[,r]<-species4$density[(r*501+1-501) : (r*501)]
  v4[,r]<-species4$V[(r*501+1-501) : (r*501)]
  
  d5[,r]<-species5$density[(r*501+1-501) : (r*501)]
  v5[,r]<-species5$V[(r*501+1-501) : (r*501)]
  
  d6[,r]<-species6$density[(r*501+1-501) : (r*501)]
  v6[,r]<-species6$V[(r*501+1-501) : (r*501)]
  
}

color.pal<-c("firebrick", "grey","darkgreen","darkblue","mediumpurple2","lightsalmon3")
potential.03<-data.frame( Potential=c(rowMeans(v1),rowMeans(v2),rowMeans(v3),rowMeans(v4),
                                      rowMeans(v5),rowMeans(v6)),
                          Density=c(rowMeans(d1),rowMeans(d2),rowMeans(d3),rowMeans(d4),
                                    rowMeans(d5),rowMeans(d6)),
                          Species = factor(c(rep("Species 1",each=501),
                                             rep("Species 2", each=501),
                                             rep("Species 3", each=501),
                                             rep("Species 4", each=501),
                                             rep("Species 5", each=501),
                                             rep("Species 6", each=501) ))) 



ggplot(data= potential.03, aes(y = Potential,x =Density , color = Species )) +geom_line(size=1.5)+
  ylab("Effective potential V(N)")+xlab("Population density")+ggtitle("A")+
  theme_cowplot()+
  theme(plot.title = element_text(size = 10,   face = "bold"),
        text = element_text(size = 10 ),
        axis.title = element_text(face="bold"),
        legend.position = "right") +
  scale_color_manual(values=color.pal )+ labs(color="Species")+ylim(c(-0.5,10))



# estimates of metrics from the potetnail curve
slope_est<-ggplot(data= potential.curve.estimates, 
                  aes(y = slope,x = factor(forcing_species), color= factor(forcing_species ))) +
  geom_boxplot(alpha=0,position= position_dodge(1))+
  geom_point(pch = 16,alpha=0.15, position = position_jitterdodge(1))+
  scale_color_manual(values = color.pal)+
  theme_cowplot()+xlab("Species")+ggtitle("A")+
  theme(plot.title = element_text(size = 10,
                                  face = "bold"),text = element_text(size = 10 ),
        axis.title = element_text(face="bold"),legend.position = "none") + 
  labs(color="Species")+ylab("Potential Slope")


dept_est<-ggplot(data= potential.curve.estimates, 
                 aes(y = -potential.depth,x = factor(forcing_species), color= factor(forcing_species ))) +
  geom_boxplot(alpha=0,position= position_dodge(1))+
  geom_point(pch = 16,alpha=0.15, position = position_jitterdodge(1))+
  scale_color_manual(values = color.pal)+
  theme_cowplot()+xlab("Species")+ggtitle("B")+
  theme(plot.title = element_text(size = 10,
                                  face = "bold"),text = element_text(size = 10 ),
        axis.title = element_text(face="bold"),legend.position = "none") + 
  labs(color="Species")+ylab("Potential depth") 
 
width_est<-ggplot(data= potential.curve.estimates, 
                  aes(y = width.potential,x = factor(forcing_species), color= factor(forcing_species ))) +
  geom_boxplot(alpha=0,position= position_dodge(1))+
  geom_point(pch = 16,alpha=0.15, position = position_jitterdodge(1))+
  scale_color_manual(values = color.pal)+
  theme_cowplot()+xlab("Species")+ggtitle("C")+
  theme(plot.title = element_text(size = 10,
                                  face = "bold"),text = element_text(size = 10 ),
        axis.title = element_text(face="bold"),legend.position = "none") + 
  labs(color="Species")+ylab("Potential width") 


diffsp<-ggplot(data= potential.curve.estimates, 
                  aes(y = diffsp,x = factor(forcing_species), color= factor(forcing_species ))) +
  geom_boxplot(alpha=0,position= position_dodge(1))+
  geom_point(pch = 16,alpha=0.15, position = position_jitterdodge(1))+
  scale_color_manual(values = color.pal)+
  theme_cowplot()+xlab("Species")+ggtitle("A")+
  theme(plot.title = element_text(size = 10,
                                  face = "bold"),text = element_text(size = 10 ),
        axis.title = element_text(face="bold"),legend.position = "none") + 
  labs(color="Species")+ylab("gains in growth through mutualism") 

source('~/Dropbox/Zurich PhD Research/1_Chapter_1/rolling_GAMs_methods.R', echo=F)

lay_out(list(slope_est, 1, 1),
        list(dept_est,1,2),
        list(width_est,1,3))



#wentzell potential

load("Wentzell_potential.RData")
(slope_est<-ggplot(data= Wentzell_potential, 
                  aes(y = slope,x = factor(forcing_species), color= factor(forcing_species ))) +
  geom_boxplot(alpha=0,position= position_dodge(1))+
  geom_point(pch = 16,alpha=0.15, position = position_jitterdodge(1))+
  scale_color_manual(values = color.pal)+
 
  theme_cowplot()+xlab("Species")+ggtitle("A")+
  theme(plot.title = element_text(size = 10,
                                  face = "bold"),text = element_text(size = 10 ),
        axis.title = element_text(face="bold"),legend.position = "none") + 
  labs(color="Species")+ylab("Potential Slope"))


(dept_est<-ggplot(data= Wentzell_potential, 
                 aes(y = -potential.depth,x = factor(forcing_species), color= factor(forcing_species ))) +
  geom_boxplot(alpha=0,position= position_dodge(1))+
  geom_point(pch = 16,alpha=0.15, position = position_jitterdodge(1))+
  scale_color_manual(values = color.pal)+
  theme_cowplot()+xlab("Species")+ggtitle("B")+
  theme(plot.title = element_text(size = 10,
                                  face = "bold"),text = element_text(size = 10 ),
        axis.title = element_text(face="bold"),legend.position = "none") + 
  labs(color="Species")+ylab("Potential depth") )

(width_est<-ggplot(data= Wentzell_potential, 
                  aes(y = width.potential,x = factor(forcing_species), color= factor(forcing_species ))) +
  geom_boxplot(alpha=0,position= position_dodge(1))+
  geom_point(pch = 16,alpha=0.15, position = position_jitterdodge(1))+
  scale_color_manual(values = color.pal)+
  theme_cowplot()+xlab("Species")+ggtitle("C")+
  theme(plot.title = element_text(size = 10,
                                  face = "bold"),text = element_text(size = 10 ),
        axis.title = element_text(face="bold"),legend.position = "none") + 
  labs(color="Species")+ylab("Potential width") )


lay_out(list(slope_est, 1, 1),
        list(dept_est,1,2),
        list(width_est,1,3))





