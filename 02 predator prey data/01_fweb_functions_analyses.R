

# function determining the time point at which species decline by 45% of its original abundance
extent.decline.foodweb<-function(data,species){
  data$Density<-cbind(data$Resource,data$Middle.sp, data$Consumer)
  n0<-mean(data$Density[1:200,species])
  percent.decline = 0.45
  n.final<-n0*(1-percent.decline)  
  
  time.point <- 200+which(round(data$Density[200:3000,species],2) < round(n.final,2))[1]
  
  
  return(time.point)
}

# function summarizing community level metrics .- commmunity ar1, community sd, community abundance, community trait
community.measures<-function(data, tipping.point, harvest.species){
  
  reps<-1
  community.trait.mean <- community.trait.sd <- community.ar1<- community.abundance <- community.abundance.CV <-community.eigen<-
    community.harvest.species<-community.abundance.sd<-array(dim=c(tipping.point,reps))
  for ( i in 1:reps){
    
    community.trait.mean[,i] <- apply(data[[i]]$trait[1:tipping.point,], 1, mean)
    community.ar1[,i]<-apply(data[[i]]$Density[1:tipping.point,],1, function(x) ar.ols(x, aic = FALSE, order.max = 1, dmean = FALSE, 
                                                                                       intercept = FALSE)$ar)
    
    community.trait.sd[,i] <- apply(data[[i]]$trait[1:tipping.point,],1, sd)
    community.abundance[,i] <- apply(data[[i]]$Density[1:tipping.point,],1, sum)
    community.abundance.CV[,i]<- apply(data[[i]]$Density[1:tipping.point,], 1, sd)
    community.abundance.sd[,i]<-apply(data[[i]]$Density[1:tipping.point,], 1, sd)
    community.eigen[,i] <- apply(data[[i]]$Density[1:tipping.point,], 1 , function(x) rda(x)$CA$u[1])
    community.harvest.species[,i] <- data[[i]]$Density[1:tipping.point,harvest.species]
    
  }
  output<-list(community.traitmean = community.trait.mean,
               community.traitSD= community.trait.sd, 
               total.abundance= community.abundance,
               CV.abundance= community.abundance.CV,
               sd.community=community.abundance.sd,
               eigen.values.community= community.eigen,
               AR.community = community.ar1,
               harvest.sp.density=community.harvest.species)
  
  return(output)
}


#variance-covariance function that estimates the covariance matrix and calculates the dominant eigenvalue

var.cov.mat<-function(windowsize, data, no.species, tipping.point){
  
  cov.mat<-array(dim=c(tipping.point, no.species, no.species))
  mw <- round(tipping.point * windowsize/100)
  dom.eigen<-numeric()
  omw <- tipping.point - mw + 1
  nMR<-array(dim=c(mw ,no.species, omw))
  for(k in 1:no.species){
    
    for (i in 1:omw) {
      Ytw <- data$Density[i:(i + mw - 1), k]
      nMR[ , k, i] <- Ytw
    }
  }
  
  
  for(t in 1:omw){
    
    for(i in 1:no.species){
      for(k in 1:no.species){
        
        cov.mat[t,i,k]<-  1/(windowsize-1)*sum(( nMR[, i, t] - mean(nMR[, i, t]))*(nMR[, k, t] -mean(nMR[, k, t])))
      }
      
    }
    dom.eigen[t]<- Re(eigen(cov.mat[t,,], only.values = T)$values[1])
  }
  
  return(dom.eigen)
}


# shift in dominant eigenvalue of covariance matrix estimated above
delta.shift.eigenvalue<-function(E){
  
  end.length<-1:length(E)
  
  eigen<- (E-mean(E))/(sd(E))
  
  smoothed.eigen<-gam_smoothing1(end.length,eigen,-1)[,6]
  #auc.eig<-auc(end.length,smoothed.eigen)
  
  time<-seq(1,length(eigen));
  tau.eigen<-cor.test(eigen, time, method="kendall", conf.level =0.95, alternative="two.sided")$estimate
  
  
  # if (auc.eig < 0){ auc.eig <- -log(-auc.eig)}else{
  #  auc.eig <- log(auc.eig)
  #}
  
  
  return(Tau.eigen=tau.eigen)
  
}


# modified early warning signal function that estimates species level EWS

genericEWS<-function(timeseries, winsize = 50, detrending = c("no", "gaussian", 
                                                              "loess", "linear", "first-diff"), bandwidth = NULL, span = NULL, 
                     degree = NULL, logtransform = FALSE, interpolate = FALSE, 
                     AR_n = FALSE, powerspectrum = FALSE) 
{
  timeseries <- data.matrix(timeseries)
  if (dim(timeseries)[2] == 1) {
    Y = timeseries
    timeindex = 1:dim(timeseries)[1]
  }
  else if (dim(timeseries)[2] == 2) {
    Y <- timeseries[, 2]
    timeindex <- timeseries[, 1]
  }
  else {
    warning("not right format of timeseries input")
  }
  if (interpolate) {
    YY <- approx(timeindex, Y, n = length(Y), method = "linear")
    Y <- YY$y
  }
  else {
    Y <- Y
  }
  if (logtransform) {
    Y <- log(Y + 1)
  }
  detrending <- match.arg(detrending)
  if (detrending == "gaussian") {
    if (is.null(bandwidth)) {
      bw <- round(bw.nrd0(timeindex))
    }
    else {
      bw <- round(length(Y) * bandwidth/100)
    }
    smYY <- ksmooth(timeindex, Y, kernel = "normal", bandwidth = bw, 
                    range.x = range(timeindex), x.points = timeindex)
    nsmY <- Y - smYY$y
    smY <- smYY$y
  }
  else if (detrending == "linear") {
    nsmY <- resid(lm(Y ~ timeindex))
    smY <- fitted(lm(Y ~ timeindex))
  }
  else if (detrending == "loess") {
    if (is.null(span)) {
      span <- 25/100
    }
    else {
      span <- span/100
    }
    if (is.null(degree)) {
      degree <- 2
    }
    else {
      degree <- degree
    }
    smYY <- loess(Y ~ timeindex, span = span, degree = degree, 
                  normalize = FALSE, family = "gaussian")
    smY <- predict(smYY, data.frame(x = timeindex), se = FALSE)
    nsmY <- Y - smY
  }
  else if (detrending == "first-diff") {
    nsmY <- diff(Y)
    timeindexdiff <- timeindex[1:(length(timeindex) - 1)]
  }
  else if (detrending == "no") {
    smY <- Y
    nsmY <- Y
  }
  mw <- round(length(Y) * winsize/100)
  omw <- length(nsmY) - mw + 1
  low <- 6
  high <- omw
  nMR <- matrix(data = NA, nrow = mw, ncol = omw)
  x1 <- 1:mw
  for (i in 1:omw) {
    Ytw <- nsmY[i:(i + mw - 1)]
    nMR[, i] <- Ytw
  }
  nARR <- numeric()
  nSD <- numeric()
  nSK <- numeric()
  rD<-numeric() #rD
  RoD<-numeric() # ROD
  nKURT <- numeric()
  nACF <- numeric()
  nDENSITYRATIO <- numeric()
  nSPECT <- matrix(0, nrow = omw, ncol = ncol(nMR))
  nCV <- numeric()
  smARall <- numeric()
  smARmaxeig <- numeric()
  detB <- numeric()
  timesss<-numeric()
  biomass<-ar1.true<-ar.squared<-l <- numeric()
  ARn <- numeric()
  nSD <- apply(nMR, 2, sd, na.rm = TRUE)
  for (i in 1:ncol(nMR)) {
    nYR <- ar.ols(nMR[, i], aic = FALSE, order.max = 1, dmean = FALSE, 
                  intercept = FALSE)
    
    l<-length(nMR[,i])
    data.level<-as.numeric(nMR[-1,i])
    data.lags<-as.numeric(nMR[-l,i])
    
    new.w<-lm(data.level ~ data.lags + I(data.lags^2/2))
    ar.squared[i]<-new.w$coefficients[3]
    ar1.true[i]<-new.w$coefficients[2]
    
    
    timesss[i]<-mean(nMR[,i], na.rm=TRUE)
    biomass[i]<- nMR[mw,i] 
    # RoD[i]<-(sqrt(sum((nMR[2:nrow(nMR),i]-nMR[1:nrow(nMR),i])^2))/(length(nMR[,i])-1))/nSD[i]
    nARR[i] <- nYR$ar
    nSK[i] <- abs(moments::skewness(nMR[, i], na.rm = TRUE))
    nKURT[i] <- moments::kurtosis(nMR[, i], na.rm = TRUE)
    nCV[i] <- nSD[i]/mean(nMR[, i])
    ACF <- acf(nMR[, i], lag.max = 1, type = c("correlation"), 
               plot = FALSE)
    nACF[i] <- ACF$acf[2]
    spectfft <- spec.ar(nMR[, i], n.freq = omw, plot = FALSE, 
                        order = 1)
    nSPECT[, i] <- spectfft$spec
    nDENSITYRATIO[i] <- spectfft$spec[low]/spectfft$spec[high]
    # if (AR_n) {
    #   ARall <- ar.ols(nMR[, i], aic = TRUE, order.max = 6, 
    #                   demean = F, intercept = F)
    #   smARall[i] <- ARall$ar[1]
    #   ARn[i] <- ARall$order
    #   roots <- Mod(polyroot(c(rev(-ARall$ar), 1)))
    #   smARmaxeig[i] <- max(roots)
    #   detB[i] <- (prod(roots))^(2/ARn[i])
    # }
  }
  timesss<-round(timesss,3)
  if(sd(timesss) == 0){
    z.stand.timess<- (timesss-mean(timesss))
  }else { z.stand.timess<- (timesss-mean(timesss))/sd(timesss)}
  
  nRETURNRATE = 1/nARR
  timevec <- seq(1, length(nARR))
  
  if(round(lm(timesss~timevec)$coefficients[2],3) < 0) 
  { z.stand.timess<- -1*z.stand.timess } else if( round(lm(z.stand.timess~timevec)$coefficients[2],1) >= 0 ){z.stand.timess <- z.stand.timess}
  
  Ktau.trait<-cor.test(timevec, z.stand.timess, alternative = c("two.sided"), 
                       method = c("kendall"), conf.level = 0.95)
  if(Ktau.trait[1] == 'NA')
  { Ktau.trait[1] <- 0}else {
    Ktau.trait[1] <-Ktau.trait[1]
  }
  KtAR <- cor.test(timevec, round(nARR,3), alternative = c("two.sided"), 
                   method = c("kendall"), conf.level = 0.95)
  KtAR.true<-cor.test(timevec, ar1.true, alternative = c("two.sided"), 
                      method = c("kendall"), conf.level = 0.95)
  KtAR2<-cor.test(timevec, ar.squared, alternative = c("two.sided"), 
                  method = c("kendall"), conf.level = 0.95)
  
  KtACF <- cor.test(timevec, round(nACF,3), alternative = c("two.sided"), 
                    method = c("kendall"), conf.level = 0.95)
  KtSD <- cor.test(timevec, round(nSD,3), alternative = c("two.sided"), 
                   method = c("kendall"), conf.level = 0.95)
  KtSK <- cor.test(timevec, nSK, alternative = c("two.sided"), 
                   method = c("kendall"), conf.level = 0.95)
  KtKU <- cor.test(timevec, nKURT, alternative = c("two.sided"), 
                   method = c("kendall"), conf.level = 0.95)
  KtDENSITYRATIO <- cor.test(timevec, nDENSITYRATIO, alternative = c("two.sided"), 
                             method = c("kendall"), conf.level = 0.95)
  KtRETURNRATE <- cor.test(timevec, nRETURNRATE, alternative = c("two.sided"), 
                           method = c("kendall"), conf.level = 0.95)
  KtCV <- cor.test(timevec,round(nCV,3), alternative = c("two.sided"), 
                   method = c("kendall"), conf.level = 0.95)
  out <- data.frame(timeindex[mw:length(nsmY)], nARR, nSD, ar.squared,ar1.true,
                    nSK, nKURT, nCV, nRETURNRATE, nDENSITYRATIO, 
                    nACF,timesss,biomass,
                    round(KtSD$estimate, digits = 3),  round(KtAR$estimate,digits = 3),
                    round(KtAR.true$estimate,digits = 3),round(KtAR2$estimate,digits = 3), round(KtACF$estimate,digits = 3),
                    round(Ktau.trait$estimate,digits = 3))
  
  
  colnames(out) <- c("timeindex", "ar1", "sd","artrue","ar2" , "sk", "kurt", 
                     "cv", "returnrate", "densratio", "acf1",
                     "ddata","adata","KTauSD", "KTauAR", "KTauAr.true","KTauAR.squared","KTauAcf", "Ktau.trait")
  return(out)
}



# effective ptential metrics such as slope, width and depth from the potential curve. 
# takes in time series data, and the species that is forced to collapse
epotential<-function(dat,forcing.species){
  

  if(forcing.species == 1 || forcing.species == 2 || forcing.species == 3){
  bi1<-numeric()
  for(i in 1:1){
    bi1[i]<-0.5*(erf((1-dat[[i]]$trait[200,forcing.species])/(sqrt(2*dat[[i]]$sa[1,forcing.species])))+
                       erf((1+dat[[i]]$trait[200,forcing.species])/(sqrt(2*dat[[i]]$sa[1,forcing.species]))))
  }
  
  Nj1<-matrix(NA,1,6)
  
  for(i in 1:1){
    Nj1[i,]<-c(dat[[i]]$Density[200,1],
               dat[[i]]$Density[200,2],
               dat[[i]]$Density[200,3],
               dat[[i]]$Density[200,4],
               dat[[i]]$Density[200,5],
               dat[[i]]$Density[200,6])
  }
  inter.comp.sp1<-intra.comp.sp1<-pred.comp.sp1<-numeric()
  V1<-list()
  slope.sp<-pot.depth.sp<-width.potential.sp<-numeric()
  N1<-seq(0,9,0.001)
  for(i in 1:1){
    inter.comp.sp1[i]<- sum(parameters$Pmatrix[forcing.species,]*Nj1[i,1:3]) - 
      parameters$Pmatrix[forcing.species,forcing.species]*Nj1[i,forcing.species]
    intra.comp.sp1[i]<-  parameters$Pmatrix[forcing.species,forcing.species]*Nj1[i,forcing.species]
    pred.comp.sp1[i]  <- sum( dat[[i]]$gammA[200,forcing.species,]*Nj1[i,4:5])
    
    
    
    V1[[i]] <-  list(  - ( bi1[i]*N1^2/2 - intra.comp.sp1[i]*N1^3/3 - inter.comp.sp1[i]*N1^2/2 +
                             pred.comp.sp1[i]*N1^2/2 ) )
    
    # metrics of stability :
    
    #1. depth of the potential
    pot.depth.sp[i]<- min(V1[[i]][[1]]) -  V1[[i]][[1]][1]
    
    
    #2. width of the potential 
    width.potential.sp[i]  <- (N1[which(V1[[i]][[1]] == min(V1[[i]][[1]]))] - N1[1])/max(N1)
    
    #3. slope of the potential
    slope.sp[i] <- pot.depth.sp[i]/width.potential.sp[i]
  }
  
  diffsp<- bi1+ pred.comp.sp1
  
  V.mean<- - ( mean(bi1)*N1^2/2 - mean(intra.comp.sp1)*N1^3/3 - mean(inter.comp.sp1)*N1^2/2 +
                  mean(pred.comp.sp1)*N1^2/2)

  } else if(forcing.species == 4 || forcing.species == 5){
    
    ## effective potential of species 4
    bi<- parameters$R[forcing.species]
    Nj4<-matrix(NA,1,6)
    # w/sqrt(2*sp2.dat.stoch[[1]]$var[2]+w^2)
    for(i in 1:1){
      Nj4[i,]<-c(dat[[i]]$Density[200,1],
                 dat[[i]]$Density[200,2],
                 dat[[i]]$Density[200,3],
                 dat[[i]]$Density[200,4],
                 dat[[i]]$Density[200,5],
                 dat[[i]]$Density[200,6])
    }
    inter.comp.sp4<-intra.comp.sp4<-pred.comp.sp4<-pred.cons.sp4<-numeric()
    slope.sp<-pot.depth.sp<-width.potential.sp<-numeric()
    
    N1<-seq(0,9,.001)
    V4<-list()
    basal.sp<-3
    top.pred<-1
    species<-forcing.species-basal.sp-top.pred+1
    for(i in 1:1){
      inter.comp.sp4[i] <- sum(parameters$A_matrix[species,]*Nj4[i,4:5]) - 
        parameters$A_matrix[species,species]*Nj4[i,4]
      intra.comp.sp4[i] <-  parameters$A_matrix[species,species]*Nj4[i,4]
      pred.comp.sp4[i]  <- -sum(dat[[i]]$gammA[200,,species]*Nj4[i,1:3])
      
      pred.cons.sp4[i]  <- sum(dat[[i]]$M.gammA[200,1]*Nj4[i,6])
      
      
      V4[[i]] <-  list(  - ( bi*N1^2/2 - intra.comp.sp4[i]*N1^3/3 - inter.comp.sp4[i]*N1^2/2 +
                               pred.comp.sp4[i]*N1^2/2 + pred.cons.sp4[i]*N1^2/2 ) )
      
      # metrics of stability :
      
      #1. depth of the potential
      pot.depth.sp[i]<- min(V4[[i]][[1]]) -  V4[[i]][[1]][1]
      
      
      #2. width of the potential 
      width.potential.sp[i]  <- (N1[which(V4[[i]][[1]] == min(V4[[i]][[1]]))] - N1[1])/max(N1)
      
      #3. slope of the potential
      slope.sp[i] <- pot.depth.sp[i]/width.potential.sp[i]
      
    }
    
    diffsp<- bi + pred.cons.sp4+pred.comp.sp4
    
    V.mean<- - ( bi*N1^2/2 - mean(intra.comp.sp4)*N1^3/3 - mean(inter.comp.sp4)*N1^2/2 +
                    mean(pred.comp.sp4)*N1^2/2 - mean(pred.cons.sp4)*N1^2/2 )
 #   lines(N1,V4.mean, typ='l', lwd=5,  col='darkblue')
    
    
  } else if(forcing.species == 6){
    
    bi=parameters$R[forcing.species]
    Nj6<-matrix(NA,1,6)
    # w/sqrt(2*sp2.dat.stoch[[1]]$var[2]+w^2)
    for(i in 1:1){
      Nj6[i,]<-c(dat[[i]]$Density[200,1],
                 dat[[i]]$Density[200,2],
                 dat[[i]]$Density[200,3],
                 dat[[i]]$Density[200,4],
                 dat[[i]]$Density[200,5],
                 dat[[i]]$Density[200,6])
    }
    inter.comp.sp6<-intra.comp.sp6<-pred.comp.sp6<-pred.cons.sp6<-numeric()
    slope.sp<-pot.depth.sp<-width.potential.sp<-numeric()
    V6<-list()
    
    N1<-seq(0,9,0.01)
    for(i in 1:1){
      pred.cons.sp6[i]  <- -sum(dat[[i]]$M.gammA[200,]*Nj6[i,4:5])
      intra.comp.sp6[i] <- 0.005*Nj6[i,6]
      V6[[i]] <-  list(  - ( bi*N1^2/2 - intra.comp.sp6[i]*N1^3/3+ pred.cons.sp6[i]*N1^2/2))
      
      # metrics of stability :
      
      #1. depth of the potential
      pot.depth.sp[i]<- min(V6[[i]][[1]]) -  V6[[i]][[1]][1]
      
      
      #2. width of the potential 
      width.potential.sp[i]  <- (N1[which(V6[[i]][[1]] == min(V6[[i]][[1]]))] - N1[1])/max(N1)
      
      #3. slope of the potential
      slope.sp[i] <- pot.depth.sp[i]/width.potential.sp[i]
      
    }
    diffsp<- bi+ pred.cons.sp6
    
    V.mean<- - ( bi*N1^2/2 - mean(intra.comp.sp6)*N1^3/3 + mean(pred.cons.sp6)*N1^2/2)
    
  }
  
  

  return(list(N1=N1, V.mean=V.mean,slope.sp=slope.sp,width.potential.sp=width.potential.sp,
              pot.depth.sp=pot.depth.sp,diffsp=diffsp))
}

is.there<-function(x=0, L.CI, U.CI){
  pos<- ifelse(x<U.CI, 1, -1) 
  negs<-ifelse(x>L.CI, -1, 1)
  return(pos+negs)}



#smoothing gam
gam_smoothing1<-function(years, timeseries,knots){
  if(length(which(timeseries<=0))==0){
    gam1<-gam(timeseries~s(as.vector(years), bs="cs", k=knots), family=gaussian(link="log"))}else{
      gam1<-gam(timeseries~s(as.vector(years), bs="cs", k=knots), family=gaussian)}
  time.series.fit<-predict(gam1, newdata=data.frame(years=years), type="response")
  
  X0<-predict(gam1, newdata=data.frame(years=years), type= "lpmatrix")
  X1<-predict(gam1, newdata=data.frame(years=years), type= "lpmatrix")
  Xi<-(X1-X0)
  df <- Xi%*%coef(gam1)              ## ith smooth derivative 
  df.sd <- rowSums(Xi%*%gam1$Vp*Xi)^.5 ## cheap diag(Xi%*%b$Vp%*%t(Xi))^.5
  #plot(years,df,type="l",ylim=range(c(df+2*df.sd,df-2*df.sd)))##plot 'em
  #lines(years,df+2*df.sd,lty=2);lines(years,df-2*df.sd,lty=2)
  splines<-data.frame(years=years,deriv=df,U.CI=df+2*df.sd, L.CI=df-2*df.sd)	
  splines$sign<-is.there(0, splines$L.CI, splines$U.CI)/2
  splines$fit<-time.series.fit
  return(splines)}




# Wentzell-Friedlin quasipotential function 
#returns the metrics of quasipotential such as slope, width and depth of the quasi potential
Congugate_gradient_method<-function(r, N, S, species_index){
  Vab<-matrix(NA, 2000, S)
  
  for( t in 1:2000){
    for(i in 1:S){
    Vab[t,i]<-sqrt(((N[t+5,i] - N[t,i]) - 0.5*(r[t+5,i]- r[t,i]))^2)^2 
    
  }
  }
  uab<-5*rowSums(Vab)
  
 
  
  Vp<-data.frame(state=   N[1:2000,species_index], Vab=uab)
  fit<- lm(Vab ~ state + I(state^2), data=Vp) 
  
 # effect_plot(fit, pred= state,plot.points = TRUE, interval = TRUE, x.label = expression(N[6], "state, N"), y.label = "V",
  #            point.alpha = 0.25, point.size = 2.5,line.thickness = 1.8 , colors = "lightsalmon3" )
  newdata <- data.frame('state'=seq(0,max(N[1:2000,species_index]),(max(N[1:2000,species_index])-0)/2000))
  p_newdat<-predict(fit, newdata,type="response")
  
  
  pot.depth.sp<- min(p_newdat)- p_newdat[1] #  min(V6[[i]][[1]]) -  V6[[i]][[1]][1]
  
  
  #2. width of the potential 
  width.potential.sp  <- abs(newdata$state[which(p_newdat == min(p_newdat)   )] - p_newdat[1])/max(newdata)  #(N1[which(V6[[i]][[1]] == min(V6[[i]][[1]]))] - N1[1])/max(N1)
  
  #3. slope of the potential
  slope.sp <- pot.depth.sp/width.potential.sp
  
  
  
  return(list(V=p_newdat,state=seq(0,max(N[1:2000,species_index]),(3-0)/2000) , pot.depth.sp=pot.depth.sp,width.potential.sp=width.potential.sp,slope.sp=slope.sp))
}



is.there<-function(x=0, L.CI, U.CI){
  pos<- ifelse(x<U.CI, 1, -1) 
  negs<-ifelse(x>L.CI, -1, 1)
  return(pos+negs)}


#stochastic environmental harvesting function
env.function<-function(name){
  
  dynamic.forcing<-numeric()
  x<-c(0.1,0.2,0.25,.4,0.8)
  variance.strength<-abs(rnorm(5000,0.2,0.1))
  sdd.strength <- seq(0, 10, by=((10-0)/(5000-1)))
  
  if(name == 'stochastic'){
    for(t in 1:5000){ 
      
      sdd<-sample(x, prob = c(0.3,0.3,0.3,0.1,0.1), size=1)
      
      
      dynamic.forcing[t]<-rnorm(1,sdd.strength[t],sd=variance.strength[t])
    }
  }
  
  #
  
  
  return(dynamic.forcing) 
}





# function for interspecific within trophic competition for prey-predator food.web module
matrix.competition<-function(middle.sp,basal.sp){
  
  
  Aspecies<- middle.sp
  Plantspecies<- basal.sp
  
  Amatrix<-matrix(runif(Aspecies^2, 0.0001, 0.001), nrow=Aspecies, ncol = Aspecies)
  diag(Amatrix)<-0.5
  
  Pmatrix<-matrix(runif(Plantspecies^2,  0.0000, 0.0001), nrow=Plantspecies, ncol = Plantspecies)
  diag(Pmatrix)<-0.01
  
  out<-return(list(Amatrix=Amatrix,Pmatrix=Pmatrix))
  
}

