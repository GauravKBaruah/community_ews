library(vegan)
library(MASS)
library(flux)
library(mgcv)


##round(length(timeseries)/4)	
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

# function for layout used for plotting
lay_out = function(...) {    
  x <- list(...)
  n <- max(sapply(x, function(x) max(x[[2]])))
  p <- max(sapply(x, function(x) max(x[[3]])))
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(n, p)))    
  
  for (i in seq_len(length(x))) {
    print(x[[i]][[1]], vp = grid::viewport(layout.pos.row = x[[i]][[2]], 
                                           layout.pos.col = x[[i]][[3]]))
  }
} 




# extract time series point where 45% decline of a species has happened
extent.decline.mutualism<-function(data,species){
  n0<-mean(data$Density[1:500,species])
  percent.decline = 0.45
  n.final<-n0*(1-percent.decline)  
  
  time.point <- 500+which(round(data$Density[500:2500,species],1) < round(n.final,1))[1]
  
  
  return(time.point)
}




#function that calculates and summarizes communiy measures such as community ar1, community trait and community abundance

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



#function for variance-covariance matrix and returns the dominant eigenvalue from this matrix
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



# calculates the shift in dominant eigenvalue from covariance matrix
delta.shift.eigenvalue<-function(E){
  
  end.length<-1:length(E)
  
  eigen<- (E-mean(E))/(sd(E))
  
  smoothed.eigen<-gam_smoothing1(end.length,eigen,-1)[,6]
#  auc.eig<-auc(end.length,smoothed.eigen)
  
  time<-seq(1,length(eigen));
  tau.eigen<-cor.test(eigen, time, method="kendall", conf.level =0.95, 
                      alternative="two.sided")$estimate
  
  
  # if (auc.eig < 0){ auc.eig <- -log(-auc.eig)}else{
  #  auc.eig <- log(auc.eig)
  #}
  
  
  return(Tau.eigen=tau.eigen)
  
}




# modified function from "early warning signals" package that estimates autocorrelation, variance etc,

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
  timesss<-round(timesss,2)
  if(sd(timesss) == 0){
    z.stand.timess<- (timesss-mean(timesss))
  }else { z.stand.timess<- (timesss-mean(timesss))/sd(timesss)}
  
  nRETURNRATE = 1/nARR
  timevec <- seq(1, length(nARR))
  if(round(lm(timesss~timevec)$coefficients[2],2) < 0) 
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






#effective potential curve from equation 8 of 1-Dimensional equation and all the metrics from the effective potential curve

grates<-function(dat,forcing.species){
  
  if(forcing.species==0){ 
    forcing.plant.no <- 0
    forcing.animal.no<-0
  } else if(forcing.species==5 ){ 
    forcing.animal.no = 6-forcing.species 
    forcing.plant.no = 0 
    bi1<-numeric()
    for(i in 1:1){
      bi1[i]<-mean(0.5*(erf((1-dat[[i]]$trait[450:500,forcing.species])/(sqrt(2*dat[[i]]$varA[forcing.animal.no])))+
                     erf((1+dat[[i]]$trait[450:500,forcing.species])/(sqrt(2*dat[[i]]$varA[forcing.animal.no])))))
    }
    
    Nj1<-matrix(NA,1,6)
    for(i in 1:1){
      Nj1[i,]<-c( mean(dat[[i]]$Density[450:500,1]),
                 mean(dat[[i]]$Density[450:500,2]),
                 mean(dat[[i]]$Density[450:500,3]),
                 mean(dat[[i]]$Density[450:500,4]),
                 mean(dat[[i]]$Density[450:500,5]),
                 mean(dat[[i]]$Density[450:500,6]))
    }
    inter.comp.sp1<-intra.comp.sp1<-mutualism.sp1<-numeric()
    #V.mean<-list()
    V1<-list()
    slope.sp<-pot.depth.sp<-width.potential.sp<-numeric()
    N1<-seq(0,5,0.01)
    
    for(i in 1:1){
      inter.comp.sp1[i]<- sum(parameters$Amatrix[forcing.animal.no,]*Nj1[i,5:6]) - 
        parameters$Amatrix[forcing.animal.no,forcing.animal.no]*Nj1[i,forcing.species]
      intra.comp.sp1[i]<-  parameters$Amatrix[forcing.animal.no,forcing.animal.no]*Nj1[i,forcing.species]
      mutualism.sp1[i]  <- sum( mean(dat[[i]]$Gamma[450:500,,forcing.animal.no])*Nj1[i,1:4])
      
      
      
      V1[[i]] <-  list(  - ( bi1[i]*N1^2/2 - intra.comp.sp1[i]*N1^3/3 - inter.comp.sp1[i]*N1^2/2 +
                               mutualism.sp1[i]*N1^2/2) )
      
      # metrics of stability :
      
      #1. depth of the potential
      pot.depth.sp[i]<- min(V1[[i]][[1]]) -  V1[[i]][[1]][1]
      
      
      #2. width of the potential 
      width.potential.sp[i]  <- (N1[which(V1[[i]][[1]] == min(V1[[i]][[1]]))] - N1[1])/max(N1)
      
      #3. slope of the potential
      slope.sp[i] <- pot.depth.sp[i]/width.potential.sp[i]
      
      
      
    }
    diffsp<-0*bi1+mutualism.sp1 #intra.comp.sp1-inter.comp.sp1
    
    V.mean<- - ( mean(bi1)*N1^2/2 - mean(intra.comp.sp1)*N1^3/3 - mean(inter.comp.sp1)*N1^2/2 +
                   mean(mutualism.sp1)*N1^2/2)
    
  } else if(forcing.species==6){
    forcing.animal.no <- 2
    forcing.plant.no<-0
    bi1<-numeric()
    for(i in 1:1){
      bi1[i]<-0.5*mean(erf((1-dat[[i]]$trait[450:500,forcing.species])/(sqrt(2*dat[[i]]$varA[forcing.animal.no])))+
                         erf((1+dat[[i]]$trait[450:500,forcing.species])/(sqrt(2*dat[[i]]$varA[forcing.animal.no]))))
    }
    
    Nj1<-matrix(NA,1,6)
    # w/sqrt(2*sp2.dat.stoch[[1]]$var[2]+w^2)
    
    for(i in 1:1){
      Nj1[i,]<-c(mean(dat[[i]]$Density[450:500,1]),
                 mean(dat[[i]]$Density[450:500,2]),
                 mean(dat[[i]]$Density[450:500,3]),
                 mean(dat[[i]]$Density[450:500,4]),
                 mean(dat[[i]]$Density[450:500,5]),
                 mean(dat[[i]]$Density[450:500,6]))
    }
    inter.comp.sp1<-intra.comp.sp1<-mutualism.sp1<-numeric()
    V1<-list()
    slope.sp<-pot.depth.sp<-width.potential.sp<-numeric()
    N1<-seq(0,5,0.01)
    
    for(i in 1:1){
      inter.comp.sp1[i]<- sum(parameters$Amatrix[forcing.animal.no,]*Nj1[i,5:6]) - 
        parameters$Amatrix[forcing.animal.no,forcing.animal.no]*Nj1[i,forcing.species]
      intra.comp.sp1[i]<-  parameters$Amatrix[forcing.animal.no,forcing.animal.no]*Nj1[i,forcing.species]
      mutualism.sp1[i]  <- sum( mean(dat[[i]]$Gamma[150:200,,forcing.animal.no])*Nj1[i,1:4])
      
      
      
      V1[[i]] <-  list(  - ( bi1*N1^2/2 - intra.comp.sp1[i]*N1^3/3 - inter.comp.sp1[i]*N1^2/2 +
                               mutualism.sp1[i]*N1^2/2) )
      
      # metrics of stability :
      
      #1. depth of the potential
      pot.depth.sp[i]<- min(V1[[i]][[1]]) -  V1[[i]][[1]][1]
      
      
      #2. width of the potential 
      width.potential.sp[i]  <- (N1[which(V1[[i]][[1]] == min(V1[[i]][[1]]))] - N1[1])/max(N1)
      
      #3. slope of the potential
      slope.sp[i] <- pot.depth.sp[i]/width.potential.sp[i]
      
      
      
    }
    
    
    diffsp<-0*bi1+mutualism.sp1 #+intra.comp.sp1-inter.comp.sp1
    
    V.mean<- - ( mean(bi1)*N1^2/2 - mean(intra.comp.sp1)*N1^3/3 - mean(inter.comp.sp1)*N1^2/2 +
                   mean(mutualism.sp1)*N1^2/2) 
    
  } else if(forcing.species == 1 ||forcing.species==3 || forcing.species==4){
    
    forcing.plant.no = forcing.species
    forcing.animal.no=0
    bi4<-0.5*mean((erf((1-dat[[1]]$trait[450:500,forcing.plant.no])/(sqrt(2*dat[[1]]$varP[forcing.plant.no])))+
                    erf((1+dat[[1]]$trait[450:500,forcing.plant.no])/(sqrt(2*dat[[1]]$varP[forcing.plant.no])))))
    
    
    
    Nj4<-matrix(NA,1,6)
    # w/sqrt(2*sp2.dat.stoch[[1]]$var[2]+w^2)
    for(i in 1:1){
      Nj4[i,]<-c(mean(dat[[i]]$Density[450:500,1]),
                 mean(dat[[i]]$Density[450:500,2]),
                 mean(dat[[i]]$Density[450:500,3]),
                 mean(dat[[i]]$Density[450:500,4]),
                 mean(dat[[i]]$Density[450:500,5]),
                 mean(dat[[i]]$Density[450:500,6]))
    }
    inter.comp.sp4<-intra.comp.sp4<-mutualism.sp4<-numeric()
    slope.sp<-pot.depth.sp<-width.potential.sp<-numeric()
    V4<-list()
    N1<-seq(0,5,0.01)
    
    for(i in 1:1){
      inter.comp.sp4[i] <- sum(parameters$Pmatrix[forcing.plant.no,]*Nj4[i,1:4]) - 
        parameters$Pmatrix[forcing.plant.no,forcing.plant.no]*Nj4[i,forcing.plant.no]
      intra.comp.sp4[i] <-  parameters$Pmatrix[forcing.plant.no,forcing.plant.no]*Nj4[i,forcing.plant.no]
      mutualism.sp4[i]  <- sum(mean(dat[[i]]$Gamma[450:500,forcing.plant.no,]*Nj4[i,5:6]))
      
      #pred.cons.sp4[i]  <- sum(sp4.dat.stoch[[i]]$M.gammA[200,1]*Nj4[i,6])
      
      
      V4[[i]] <-  list(  -  ( bi4*N1^2/2 - intra.comp.sp4[i]*N1^3/3 - inter.comp.sp4[i]*N1^2/2 +
                                mutualism.sp4[i]*N1^2/2) )
      
      # metrics of stability :
      
      #1. depth of the potential
      pot.depth.sp[i]<- min(V4[[i]][[1]]) -  V4[[i]][[1]][1]
      
      
      #2. width of the potential 
      width.potential.sp[i]  <- (N1[which(V4[[i]][[1]] == min(V4[[i]][[1]]))] - N1[1])/max(N1)
      
      #3. slope of the potential
      slope.sp[i] <- pot.depth.sp[i]/width.potential.sp[i]
      
    }
    
    diffsp<- 0*bi4+ mutualism.sp4 #- intra.comp.sp4-inter.comp.sp4
    
    V.mean<- - ( mean(bi4)*N1^2/2 - mean(intra.comp.sp4)*N1^3/3 - mean(inter.comp.sp4)*N1^2/2 +
                   mean(mutualism.sp4)*N1^2/2)
    
  } else if(forcing.species ==2 ){
    bi4<- numeric()
    forcing.plant.no = forcing.species
    forcing.animal.no=0
    bi4<-0.5*mean(erf((1-dat[[1]]$trait[450:500,forcing.plant.no])/(sqrt(2*dat[[1]]$varP[forcing.plant.no])))+
                erf((1+dat[[1]]$trait[450:500,forcing.plant.no])/(sqrt(2*dat[[1]]$varP[forcing.plant.no]))))
    
    
    
    Nj4<-matrix(NA,1,6)
    # w/sqrt(2*sp2.dat.stoch[[1]]$var[2]+w^2)
    for(i in 1:1){
      Nj4[i,]<-c(mean(dat[[i]]$Density[450:500,1]),
                 mean(dat[[i]]$Density[450:500,2]),
                 mean(dat[[i]]$Density[450:500,3]),
                 mean(dat[[i]]$Density[450:500,4]),
                 mean(dat[[i]]$Density[450:500,5]),
                 mean(dat[[i]]$Density[450:500,6]))
    }
    inter.comp.sp4<-intra.comp.sp4<-mutualism.sp4<-numeric()
    slope.sp<-pot.depth.sp<-width.potential.sp<-numeric()
    V4<-list()
    N1<-seq(0,5,0.01)
    
    for(i in 1:1){
      inter.comp.sp4[i] <- sum(parameters$Pmatrix[forcing.plant.no,]*Nj4[i,1:4]) - 
        parameters$Pmatrix[forcing.plant.no,forcing.plant.no]*Nj4[i,forcing.plant.no]
      intra.comp.sp4[i] <-  parameters$Pmatrix[forcing.plant.no,forcing.plant.no]*Nj4[i,forcing.plant.no]
      mutualism.sp4[i]  <- sum(mean(dat[[i]]$Gamma[450:500,forcing.plant.no,]*Nj4[i,5:6]))
      
      #pred.cons.sp4[i]  <- sum(sp4.dat.stoch[[i]]$M.gammA[200,1]*Nj4[i,6])
      
      
      V4[[i]] <-  list(  -  ( bi4*N1^2/2 - intra.comp.sp4[i]*N1^3/3 - inter.comp.sp4[i]*N1^2/2 +
                                mutualism.sp4[i]*N1^2/2) )
      
      # metrics of stability :
      
      #1. depth of the potential
      pot.depth.sp[i]<- min(V4[[i]][[1]]) -  V4[[i]][[1]][1]
      
      
      #2. width of the potential 
      width.potential.sp[i]  <- (N1[which(V4[[i]][[1]] == min(V4[[i]][[1]]))] - N1[1])/max(N1)
      
      #3. slope of the potential
      slope.sp[i] <- pot.depth.sp[i]/width.potential.sp[i]
      
    }
    
    diffsp<-  mutualism.sp4 #- intra.comp.sp4-inter.comp.sp4
    
    V.mean<- - ( mean(bi4)*N1^2/2 - mean(intra.comp.sp4)*N1^3/3 - mean(inter.comp.sp4)*N1^2/2 +
                   mean(mutualism.sp4)*N1^2/2)
    
  }
  
  
  
  return(list(N1=N1, V.mean=V.mean,slope.sp=slope.sp,width.potential.sp=width.potential.sp,
              pot.depth.sp=pot.depth.sp,diffsp=diffsp))
}




# Wentzell potential function from stochastic deviation theory. returns the stability-landscape metrics

Congugate_gradient_method<-function(r, N, S, species_index){
  Vab<-matrix(NA, 4000, S)
  
  for( t in 1:4000){
    for(i in 1:S){
      Vab[t,i]<-sqrt(((N[t+1,i] - N[t,i])/5 - 0.5*(r[t+1,i]- r[t,i]))^2)^2 
      
    }
  }
  uab<-5*rowSums(Vab)
  
  
  
  Vp<-data.frame(state=   N[1:4000,species_index], Vab=uab)
  fit<- lm(Vab ~ state + I(state^2), data=Vp) 
  newdata <- data.frame('state'=seq(0,max(N[1:4000,species_index]),(5-0)/4000))
  p_newdat<-predict(fit, newdata,type="response")
  
  
  pot.depth.sp<- min(p_newdat)- p_newdat[1] #  min(V6[[i]][[1]]) -  V6[[i]][[1]][1]
  
  
  #2. width of the potential 
  width.potential.sp  <- abs(newdata$state[which(p_newdat == min(p_newdat)   )] - p_newdat[1])/max(newdata)  #(N1[which(V6[[i]][[1]] == min(V6[[i]][[1]]))] - N1[1])/max(N1)
  
  #3. slope of the potential
  slope.sp <- pot.depth.sp/width.potential.sp
  
  
  
  return(list(V=p_newdat, pot.depth.sp=pot.depth.sp,width.potential.sp=width.potential.sp,slope.sp=slope.sp))
}




#stochastic fold harvesting function
# the dynamic strength of forcing is not constant but stochastically varies over time
env.function<-function(name){
  #tt<-seq(0,9000,2)
  dynamic.forcing<-numeric()
  x<-c(1,0.8,0.5,0.2,0.1,0.05)
  variance.strength<-abs(rnorm(5000,0.2,0.1))
  sdd.strength <- seq(0, 10, by=((10-0)/(5000-1)))
  
  if(name == 'stochastic'){
    for(t in 1:5000){ 
      
      sdd<-sample(x, prob = c(0.1,0.1,0.1,0.2,0.2,0.3), size=1)
      
      
      dynamic.forcing[t]<-rnorm(1,sdd.strength[t],sd=variance.strength[t])
    }
  }
  
  else if (name=='linear'){
    for(t in 1:3000){ 
      
      dynamic.forcing[t]<-rnorm(1,sdd.strength[t],sd=0.1)  
    }
  }
  
  
  return(dynamic.forcing) 
}


expand.matrix <- function(A){
  m <- nrow(A)
  n <- ncol(A)
  B <- matrix(0,nrow = m, ncol = m)
  C <- matrix(0,nrow = n, ncol = n)
  cbind(rbind(B,t(A)),rbind(A,C))
}


#error function
erf<-function(x) {return(2 * pnorm(x * sqrt(2)) - 1)}



#function that gives the competition interaction matrix 
#returns competition matrix for two guilds of species for the mutualistic community
matrix.competition<-function(matrix){
  
  
  Aspecies<- dim(matrix)[2]
  Plantspecies<- dim(matrix)[1]
  
  Amatrix<-matrix(runif(Aspecies^2, 0.0001, 0.001), nrow=Aspecies, ncol = Aspecies)
  diag(Amatrix)<-0.5
  
  Pmatrix<-matrix(runif(Plantspecies^2, 0.0001, 0.001), nrow=Plantspecies, ncol = Plantspecies)
  diag(Pmatrix)<-0.5
  
  out<-return(list(Amatrix=Amatrix,Pmatrix=Pmatrix))
  
}




