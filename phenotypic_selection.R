library("gsg")
library("mgcv")
library("data.table")

addalpha <- function(colors, alpha=1.0) {
  r <- col2rgb(colors, alpha=T)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}

#Read BLUPs
bf<-fread("/Users/pelle/Dropbox/SwAsp/GWAS/traits/BLUPs/SwAsp94_phenology.txt",head=T)
blups<-fread("/Users/pelle/Dropbox/SwAsp/GWAS/traits/BLUPs/SwAsp94_all_traits.txt",head=T)
blups$pop<-trunc((as.numeric(substr(blups$id,6,9))-1)/10)+1
blups<-cbind(blups,bf[,c(13,14)])

##Calculate mean fitness
Wbar<-function(x){
  return((x-min(x))/(max(x)-min(x)))
}

#Estimate fitness landscapes
pdf("/Users/pelle/Dropbox/SwAsp/projects/climate_adaptation/figures and tables/pdf/wbar.pdf",width=8,height=8)

par(mfrow=c(2,2))
cols<-rep("grey75",94)
cols[blups$pop<7]<-addalpha(c("indianred1"),0.9)
cols[blups$pop>8]<-addalpha(c("lightsteelblue1"),0.9)

sim.data<-NULL
for(i in seq(-1.5,2,0.1)){
  for(j in seq(-3,2.5,0.1)){
    sim.data<-rbind(sim.data,c(i,j))
  }
}
sim.data<-as.data.frame(sim.data)
names(sim.data)<-c("z1","z2")

dfE1<-NULL
dfS1$W<-Wbar(blups$heightE09)
dfE$z<-scale(blups$bsE06)
dfE1$pop<-blups$pop
dfE1$anc<-as.numeric((blups$pop>6))
dfE1<-as.data.frame(dfE1)
dfE1$sites<-1

mE1<-gam(W~s(z1),data=dfE1)
fE1<-fitness.landscape(mod=mE1,phenotype=c("z1"),points=sim.data,PI.method="boot.para",n.boot=1000)
plot(fE1$points[,1],fE1$Wbar,type="l",ylim=c(0,1),las=1,xlab="Bud set",ylab="Fitness",main="South common garden (Ekebo)")
points(W~z,dfE1,bg=cols,pch=21,cex=1.2)
lines(fE1$points[,1],fE1$WbarPI[2,],lty=3)
lines(fE1$points[,1],fE1$WbarPI[1,],lty=3)
gradE1<-gam.gradients(mod=mE1,phenotype=c("z1"),se.method="boot.para",n.boot=1000,standardized=T)

dfS1<-NULL
dfS1$W<-Wbar(blups$heightS08)
dfS1$z2<-scale(blups$bfS08)
dfS1$pop<-blups$pop
dfS1$anc<-as.numeric((blups$pop>6))
dfS1<-as.data.frame(dfS1)
dfS1$sites<-0

mS1<-gam(W~s(z1)+s(z2)+s1*z2,data=dfS1)
fS1<-fitness.landscape(mod=mS1,phenotype=c("z"),points=as.data.frame(seq(-1.5,1.5,0.05),col.names="z"),PI.method="boot",n.boot=1000)
plot(fS1$points[,1],fS1$Wbar,type="l",ylim=c(0,1),las=1,xlab="Bud set",ylab="Fitness",main="North common garden (Savar)")
points(W~z,dfS1,bg=cols,pch=21,cex=1.2)
lines(fS1$points[,1],fS1$WbarPI[2,],lty=3)
lines(fS1$points[,1],fS1$WbarPI[1,],lty=3)
gradS1<-gam.gradients(mod=mS1,phenotype=c("z1"),se.method="boot.para",n.boot=1000,standardized=T)

dfE2<-NULL
dfE2$W<-Wbar(blups$heightE09)
dfE2$z<-scale(blups$bfE05)
dfE2$pop<-blups$pop
dfE2$anc<-as.numeric((blups$pop>6))
dfE2<-as.data.frame(dfE2)
dfE2$sites<-1

mE2<-gam(W~s(z),data=dfE2)
fE2<-fitness.landscape(mod=mE2,phenotype=c("z"),points=as.data.frame(seq(-3,2.5,0.1),col.names="z"),PI.method="boot",n.boot=1000)
plot(fE2$points[,1],fE2$Wbar,type="l",ylim=c(0,1),las=1,xlab="Bud flush",ylab="Fitness",main="")
points(W~z,dfE2,bg=cols,pch=21,cex=1.2)
lines(fE2$points[,1],fE2$WbarPI[2,],lty=3)
lines(fE2$points[,1],fE2$WbarPI[1,],lty=3)
gradE2<-gam.gradients(mod=mE2,phenotype=c("z"),se.method="boot.para",n.boot=250,standardized=T)

dfS2<-NULL
dfS2$W<-Wbar(blups$heightS08)
dfS2$z<-scale(blups$bfS08)
dfS2$pop<-blups$pop
dfS2$anc<-as.numeric((blups$pop>6))
dfS2<-as.data.frame(dfS2)
dfS2$sites<-0

mS2<-gam(W~s(z),data=dfS2)
fS2<-fitness.landscape(mod=mS2,phenotype=c("z"),points=as.data.frame(seq(-2,2.5,0.05),col.names="z"),PI.method="boot",n.boot=1000)
plot(fS2$points[,1],fS2$Wbar,type="l",ylim=c(0,1),las=1,xlab="Bud flush",ylab="Fitness",main="")
points(W~z,dfS2,bg=cols,pch=21,cex=1.2)
lines(fS2$points[,1],fS2$WbarPI[2,],lty=3)
lines(fS2$points[,1],fS2$WbarPI[1,],lty=3)
gradS2<-gam.gradients(mod=mS2,phenotype=c("z"),se.method="boot.para",n.boot=250,standardized=T)

dev.off()
