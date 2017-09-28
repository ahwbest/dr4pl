
library(drc)

### When drc itself fails
setwd("C:\\Users\\Hyowon\\Presentation\\Dose_response_modelling")
setEPS()
postscript("PL_drc_error.eps",width=8,height=6, horizontal = FALSE)

par(mar=c(5,5,3,2)+0.1,mfrow=c(2,2))

setwd("C:\\Users\\Hyowon\\Research\\Dose_response_modelling\\Data")
df.drm <- read.csv(paste("drc_error1.csv",sep=""))

plot(x=log10(df.drm$Dose),
     y=df.drm$Response,
     xlab="",
     ylab="Response",
     main="(a)",
     lwd=2,
     pch=16,
     type="p",
     bty="l",
     xaxt="n")

mtext(text="Dose",side=1,line=4, cex = 0.7)

at.x=unique(log10(df.drm$Dose))
axis(1,at=at.x,labels=10^at.x,las=2)

df.drm <- read.csv(paste("drc_error2.csv",sep=""))

plot(x=log10(df.drm$Dose),
     y=df.drm$Response,
     xlab="",
     ylab="Response",
     main="(b)",
     lwd=2,
     pch=16,
     type="p",
     bty='l',
     xaxt="n")
mtext(text="Dose",side=1,line=4, cex = 0.7)

at.x=unique(log10(df.drm$Dose))
axis(1,at=at.x,labels=10^at.x,las=2)

df.drm <- read.csv(paste("drc_error3.csv",sep=""))

plot(x=log10(df.drm$Dose),
     y=df.drm$Response,
     xlab="",
     ylab="Response",
     main="(c)",
     lwd=2,
     pch=16,
     type="p",
     bty='l',
     xaxt="n")
mtext(text="Dose",side=1,line=4, cex = 0.7)

at.x=unique(log10(df.drm$Dose))
axis(1,at=at.x,labels=10^at.x,las=2)

dev.off()       


### When constrained optimization fails
setwd("C:\\Users\\Hyowon\\Presentation\\Dose_response_modelling")
setEPS()
postscript("4PL_constraint.eps",width=8,height=4)

par(mar=c(5,5,4,2)+0.1,mfrow=c(1,2))
for(i in 1:2) {
  setwd("C:\\Users\\Hyowon\\Research\\Dose_response_modelling\\Data")
  
  df.drm <- read.csv(paste("drc_constraint",i,".csv",sep=""))
  
  drm.plot <- drm(Response~Dose,
                  data=df.drc,
                  fct=LL.4(names=c("Slope","Lower Limit","Upper Limit","IC50")))
  
   drm.plot <- drm(Response~Dose,
                   data=df.drm,
                   fct=LL.4(names=c("Slope","Lower Limit","Upper Limit","IC50")),
                   lowerl=c(-Inf,0,-Inf,-Inf),
                   upperl=c(Inf,Inf,Inf,Inf))
  
  if(i==1) {
    char.main <- "(a)"
  }else {
    char.main <- "(b)"
  }
  plot(drm.plot,
       xlab="Dose",
       ylab="",
       main=char.main,
       xlim=c(0,max(df.drm$Dose)*10),
       ylim=c(0,drm.plot$coefficients[3]*1.3),
       lwd=2,
       legend=FALSE,
       pch=16,
       type="all",
       bty='l',
       yaxs='i')
  mtext(text="Response",side=2,line=4)
  legend(x="topright",
         legend=paste("Lower Limit=",round(drm.plot$coefficients[2],2),sep=""),
         bty="n")
}

dev.off()       

