library("devtools")
library("data.table")
library("hierfstat")
library("gdm")
library("rgdal")
library("rgeos")
library("rworldmap")
library("RColorBrewer")
library("scales")
library("maptools")
library("fossil")
library("vegan")
library("envirem")
library("geosphere")
load("~/Dropbox/SwAsp/projects/climate_adaptation/analyses/disimilarity_modelling/disimilarity_modelling_new.RData")

##Read WorldClim data
midH<-stack(list.files(path="~/Dropbox/RTools/envirem/Europe_holo_ccsm4_2.5arcmin",pattern=".tif",full.names=T))
current<-stack(list.files(path="~/Dropbox/RTools/envirem/Europe_current_2.5arcmin",pattern=".tif",full.names=T))
r70.45<-stack(list.files(path="~/Dropbox/RTools/envirem/Scandinavia_ccsm4_rcp45_2070_2.5arcmin",pattern=".tif",full.names=T))
r70.85<-stack(list.files(path="~/Dropbox/RTools/envirem/Scandinavia_ccsm4_rcp85_2070_2.5arcmin",pattern=".tif",full.names=T))

##Read SwAsp population data
env.pop<-read.table("~/Dropbox/SwAsp/SwAsp_Resequencing/asp201_94/environmental_variation/raw_data/SwAsp_Clone_Origin_data_bypop.txt",head=T)

##Populations
SwAsp.coords<-env.pop[,c(3,2)]
names(SwAsp.coords)<-c("x","y")

##Common gardens
CG<-data.frame(site=c("Ekebo","Savar"),longitude=c(13.11,20.55),latitude=c(55.95,63.90))
CG.coords<-CG[,c(2,3)]
names(CG.coords)<-c("x","y")
points<-SpatialPoints(CG.coords,proj4string = current@crs)
CG.env<-extract(current,points,df=T)

##Current climate
points<-SpatialPoints(SwAsp.coords,proj4string = current@crs)
SwAsp.env<-extract(current,points,df=T)
SwAsp.env<-cbind(SwAsp.coords,SwAsp.env)
SwAsp.env<-SwAsp.env[,c(3,1,2,4:19)]
names(SwAsp.env)[c(1,2,3)]<-c("Pop","Longitude","Latitude")
names(SwAsp.env)[4:19]<-substr(names(SwAsp.env)[4:19],19,50)

##Common garden
CG.coords<-CG[,c(2,3)]
names(CG.coords)<-c("x","y")
points<-SpatialPoints(CG.coords,proj4string = r70.85@crs)
CG.env<-extract(r70.85,points,df=T)
CG.env$ID<-c("Ekebo","Savar")
                
##Create raster data for Sweden
#Using shape file
#Sweden shapefile
e<-extent(c(5.5,30,54,71))
data(wrld_simpl)
SWE<-subset(wrld_simpl, NAME %in% c("Sweden"))

##P. tremula distribution
ptremula <- readOGR("~/Dropbox/RTools/Populus_tremula/Populus_tremula_EUFORGEN.shp")
ptremula<-crop(ptremula,extent(e))

##Current climate data
r.SwAsp<-crop(current,extent(SWE))
r.SwAsp<-mask(r.SwAsp,SWE)
r.SwAsp<-crop(r.SwAsp,extent(ptremula))
r.SwAsp<-mask(r.SwAsp,ptremula)
names(r.SwAsp)<-substr(names(r.SwAsp),19,50)
r.SwAsp.table<-as.data.frame(rasterToPoints(r.SwAsp))

##RCP4.5
r.SwAsp70.45<-crop(r70.45,extent(SWE))
r.SwAsp70.45<-mask(r.SwAsp70.45,SWE)
r.SwAsp70.45<-crop(r.SwAsp70.45,extent(ptremula))
r.SwAsp70.45<-mask(r.SwAsp70.45,ptremula)
names(r.SwAsp70.45)<-substr(names(r.SwAsp70.45),28,50)

##RCP8.5
r.SwAsp70.85<-crop(r70.85,extent(SWE))
r.SwAsp70.85<-mask(r.SwAsp70.85,SWE)
r.SwAsp70.85<-crop(r.SwAsp70.85,extent(ptremula))
r.SwAsp70.85<-mask(r.SwAsp70.85,ptremula)
names(r.SwAsp70.85)[1:3]<-substr(names(r.SwAsp70.85)[1:3],29,50)
names(r.SwAsp70.85)[4:16]<-substr(names(r.SwAsp70.85)[4:16],28,50)

##Mid-Holocene
#r.SwAspMidH<-crop(midH,extent(SWE))
#r.SwAspMidH<-mask(r.SwAspMidH,SWE)
#r.SwAspMidH<-crop(r.SwAspMidH,extent(ptremula))
#r.SwAspMidH<-mask(r.SwAspMidH,ptremula)
#names(r.SwAspMidH)<-names(r.SwAsp)

##Plot current climate data
col1<-addalpha("indianred3",alpha=0.9)
col2<-addalpha("indianred1",alpha=0.9)
col3<-addalpha(rgb(1,0.85,0.85),alpha=0.9)
col4<-addalpha(rgb(0.85,0.95,1.0),alpha=0.9)
col5<-addalpha("lightsteelblue2",alpha=0.9)
col6<-addalpha("lightsteelblue4",alpha=0.9)
cols<-rev((colorRampPalette(c(col1,col2,col3,col4,col5,col6),alpha=TRUE)(100)))

cols=terrain.colors(100,alpha=1)

pdf(file="~/Dropbox/SwAsp/projects/climate_adaptation/figures and tables/pdf/current_climate.pdf",width=12,height=8)
par(mfrow=c(1,2))
plot(ptremula,col="grey95",border="grey60",main="Continentality")
plot(r.SwAsp[[9]],add=T,col=cols,zlim=c(16,30))
points(CG[,c(2,3)],pch=21,bg="grey85",cex=1.5)
text(CG[,c(2,3)],labels=c("Savar", "Ekebo"),cex=0.85,adj=c(0.5,1.5))

plot(ptremula,col="grey99",border="grey60",main="Climatic moisture index")
plot(r.SwAsp[[8]],add=T)
text(env.pop[,c(3,2)],labels=as.character(env.pop$Pop),cex=1.2)
dev.off()

###Read significant SNPs 
p.snps<-scan("~/Dropbox/SwAsp/projects/climate_adaptation/analyses/disimilarity_modelling/climate_significant_clumps.txt",what="character")
##read genotype data and kinship matrix
geno<-fread("gunzip -c ~/Dropbox/SwAsp/GWAS/gemma/data/SwAsp_94samples.filter.maf005.pos.recode.gt.beagle.traw.gz",head=T)
##LD pruned data
geno<-fread("gunzip -c ~/Dropbox/SwAsp/GWAS/gemma/data/SwAsp_94_pruned.20k.2k.0.9.traw.gz",head=T)

##Read SwAsp environmental variables
env.var<-read.table("~/Dropbox/SwAsp/SwAsp_Resequencing/asp201_94/environmental_variation/raw_data/SwAsp_Clone_Origin_data_94_gwas_samples.txt",head=T)
env.var$Pop<- trunc((env.var$Clone-1)/10)+1

##Read SwAsp population data
env.pop<-read.table("~/Dropbox/SwAsp/SwAsp_Resequencing/asp201_94/environmental_variation/raw_data/SwAsp_Clone_Origin_data_bypop.txt",head=T)

##Make genotype data compatible with hierfstat
fstat.data<-t(geno[geno$SNP %in% p.snps,-c(1:6)])
fstat.data[fstat.data==0]<-11
fstat.data[fstat.data==1]<-12
fstat.data[fstat.data==2]<-11
fstat.data<-as.data.frame(cbind(env.var$Pop,fstat.data))

##Read pairwise Fst for all SNPs
fst.all<-read.table("~/Dropbox/SwAsp/SwAsp_Resequencing/asp201_94/Fst/pairwise_fst.txt",head=T)

##Caclulate pairwise Fst for selected markers
fst.selected<-pairwise.WCfst(fstat.data)

##Prepare Fst distance matrices 
D.ref<-matrix(rep(NA,144),ncol=12,nrow=12)
D.ref[lower.tri(D.ref)]<-fst.all$fst.wc

D.selected<-fst.selected

##Test for isolation by distance
D.spatial<-earth.dist(SwAsp.env[,c("Longitude","Latitude")],dist=T)
D.env<-vegdist(SwAsp.env[,c(11,12)],method="euclidian")

#Transform Fst to Fst/(1-Fst)
transform.fst<-function(x){
  return(x/(1-x))
}
Fst.ref<-D.ref
Fst.ref[lower.tri(Fst.ref)]<-transform.fst(Fst.ref[lower.tri(Fst.ref)])
Fst.ref<-as.dist(Fst.ref)

Fst.selected<-D.selected
Fst.selected[lower.tri(Fst.selected)]<-transform.fst(Fst.selected[lower.tri(Fst.selected)])
Fst.selected<-as.dist(Fst.selected)

mantel.ref<-mantel(D.ref,D.spatial)
mantel.selected<-mantel(D.selected,D.spatial)

##Plot genetic vs spatial/environmental distance
pdf("~/Dropbox/SwAsp/projects/climate_adaptation/figures and tables/pdf/ibd.pdf",width=12,height=8)
par(mfcol=c(1,2))
par(las=1)
par(mar=c(5,5,4,1))
col1<-addalpha(c("lightsteelblue1"),0.9)
plot(Fst.ref~D.spatial,bg=col1,pch=21,cex=1.3,xlab="Geographic distance (km)",ylab="",main="Reference SNPs")
title(ylab="Genetic distance (Fst/(1-Fst)", mgp = c(3.5, 1, 0))
abline(lm(Fst.ref~D.spatial),col="grey75")

col2<-addalpha(c("lightsteelblue3"),0.5)
plot(Fst.ref~D.env,bg=col2,pch=21,cex=1.3,0.5,xlab="Environmental distance",ylab="")
title(ylab="Genetic distance (Fst/(1-Fst)", mgp = c(3.5, 1, 0))
abline(lm(Fst.ref~D.env),col="grey75")

col1<-addalpha(c("indianred1"),0.9)
plot(Fst.selected~D.spatial,bg=col1,pch=21,cex=1.3,xlab="Geographic distance (km)",ylab="",main="Associated SNPs")
title(ylab="Genetic distance (Fst/(1-Fst)", mgp = c(3.5, 1, 0))
abline(lm(Fst.selected~D.spatial),col="grey75")

col2<-addalpha(c("indianred3"),0.5)
plot(Fst.selected~D.env,bg=col2,pch=21,cex=1.3,xlab="Environmental distance",ylab="")
title(ylab="Genetic distance (Fst/(1-Fst)", mgp = c(3.5, 1, 0))
abline(lm(Fst.selected~D.env),col="grey75")
dev.off()

##Scale Fst to 0-1 and format data 
scale.dist<-function(x){
  return((x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)))
}
D.ref<-scale.dist(Fst.ref)
D.selected<-scale.dist(Fst.selected)
Fstmax<-max(Fst.selected)
Fstmin<-min(Fst.selected)

#Prepare data for GDM
##Environmental variables to keep
keep<-which(names(SwAsp.env) %in% c("climaticMoistureIndex","continentality","growingDegDays5"))

refdist<-data.frame(cbind(SwAsp.env$Pop,as.matrix(D.ref)))
names(refdist)[1]<-"Pop"
gdm.ref<-formatsitepair(refdist,bioFormat=3,XColumn="Longitude",YColumn="Latitude",predData=SwAsp.env[,c(1:3,keep)],siteColumn="Pop")

seldist<-data.frame(cbind(SwAsp.env$Pop,as.matrix(D.selected)))
names(seldist)[1]<-"Pop"
gdm.selected<-formatsitepair(seldist,bioFormat=3,XColumn="Longitude",YColumn="Latitude",predData=SwAsp.env[,c(1:3,keep)],siteColumn="Pop")

##Fit GDM models
pdf("~/Dropbox/SwAsp/projects/climate_adaptation/figures and tables/pdf/gdm_ref.pdf",width=6,height=20)
gdm.ref.fit<-gdm(gdm.ref,geo=T)
gdm.ref.fit$explained
plot(gdm.ref.fit,plot.layout=c(5,1))
dev.off()

pdf("~/Dropbox/SwAsp/projects/climate_adaptation/figures and tables/pdf/gdm_selected.pdf",width=6,height=20)
gdm.selected.fit<-gdm(gdm.selected,geo=T)
gdm.selected.fit$explained
plot(gdm.selected.fit,plot.layout=c(5,1))
dev.off()

##Test environmental variables one by one
expl2=vector(mode = "integer", length = ncol(SwAsp.env)-3)
for (i in c(4:19)){
  print(c(i-3,names(SwAsp.env[i])))
  gdm.test=formatsitepair(seldist, bioFormat=3, XColumn="Longitude", YColumn="Latitude",predData=SwAsp.env[c(1:3,i)], siteColumn="Pop")
  gdm.test.fit=gdm(gdm.test, geo=T)
  expl2[i-3]=gdm.test.fit$explained
}
cbind(c(1:16),expl2)

##Transform GDM
keep2<-which(names(r.SwAsp.table) %in% c("climaticMoistureIndex","continentality","growingDegDays5"))
gdm.ref.trans<-gdm.transform(gdm.ref.fit,r.SwAsp.table[,c(1,2,keep2)])
gdm.selected.trans<-gdm.transform(gdm.selected.fit,r.SwAsp.table[,c(1,2,keep2)])

##Modified pcaToRaster function from Fitzpatrick & Keller (2015)
pcaToRaster <- function(snpTrans, x, y){
  pca <- prcomp(snpTrans, center=TRUE, scale.=FALSE)

  ##assigns to colors, edit as needed to maximize color contrast, etc.
  a1 <- pca$x[,1]; a2 <- pca$x[,2]; a3 <- pca$x[,3]
  r <- a1+a2; g <- -a2; b <- a3+a2-a1
  ##scales colors
  scalR <- (r-min(r))/(max(r)-min(r))*255
  scalG <- (g-min(g))/(max(g)-min(g))*225
  scalB <- (b-min(b))/(max(b)-min(b))*255

  ##assigns color to raster
  ##Taken from https://stackoverflow.com/questions/19627344/how-to-create-a-raster-from-a-data-frame-in-r
  outrast<-as.data.frame(cbind(x,y,scalR,scalG,scalB))
  return(rasterFromXYZ(outrast))
}

##Plot spatial variation in genetic composition
ref.pca<-pcaToRaster(gdm.ref.trans, r.SwAsp.table$x, r.SwAsp.table$y)
plotRGB(ref.pca,r=3,g=2,b=1)

selected.pca<-pcaToRaster(gdm.selected.trans, r.SwAsp.table$x, r.SwAsp.table$y)
plotRGB(selected.pca,r=1,g=2,b=3)

##RGBdiffMap from Fitzpatrick & Keller (2015)
RGBdiffMap <- function(predMap1, predMap2, x, y){
  require(vegan)
  PCA1 <- prcomp(predMap1, center=TRUE, scale.=F)
  PCA2 <- prcomp(predMap2, center=TRUE, scale.=F)
  diffProcrust <- procrustes(PCA1, PCA2, scale=F, symmetrical=FALSE)
  residMap <- residuals(diffProcrust)
  return(rasterFromXYZ(as.data.frame(cbind(x,y,residMap))))
}

##Add alpha channel to colors, from https://github.com/mylesmharrison/colorRampPaletteAlpha
addalpha <- function(colors, alpha=1.0) {
  r <- col2rgb(colors, alpha=T)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}

##Set up color palette
cols<-rev(colorRampPalette(brewer.pal(11,"Spectral"),alpha=T)(100))

col1<-addalpha("indianred4",alpha=0.8)
col2<-addalpha("indianred1",alpha=0.8)
col3<-addalpha("ivory1",alpha=0.8)
col4<-addalpha("indianred1",alpha=0.4)
col5<-addalpha("lightsteelblue2",alpha=0.8)
col6<-addalpha("lightsteelblue4",alpha=0.8)
col7<-addalpha("snow1",alpha=0.8)
col8<-addalpha("lightsteelblue2",alpha=0.4)

cols2<-rev((colorRampPalette(c(col7,col3,col8,col5,col6),alpha=TRUE)(100)))
cols1<-(colorRampPalette(c(col7,col3,col4,col2,col1),alpha=T)(100))

##Plot difference in genetic turnover
genetic.turnover<-RGBdiffMap(gdm.ref.trans,gdm.selected.trans, r.SwAsp.table$x, r.SwAsp.table$y)
plot(ptremula,col="grey95",border="grey60")
plot(genetic.turnover,col=cols,axes=F,box=F,add=T)

#Predict genetic changes to future climate
gdm.selected.pred.midH<-predict.gdm(gdm.selected.fit,r.SwAsp[[keep]],time=T,predRasts=r.SwAspMidH[[keep]])-0.0626508
plot(ptremula,col="grey95",border="grey60")
plot(gdm.selected.pred.midH,col=cols,axes=F,box=F,add=T)

keep3<-which(names(r.SwAsp) %in% c("climaticMoistureIndex","continentality","growingDegDays5"))
keep4<-which(names(r.SwAsp70.85) %in% c("climaticMoistureIndex","continentality","growingDegDays5"))

gdm.selected.pred.45<-predict.gdm(gdm.selected.fit,r.SwAsp[[keep3]],time=T,predRasts=r.SwAsp70.45[[keep3]])
gdm.selected.pred.85<-predict.gdm(gdm.selected.fit,r.SwAsp[[keep3]],time=T,predRasts=r.SwAsp70.85[[keep4]])

##Take average across all 2070 scenarios
gdm.selected.pred<-as.data.frame(rasterToPoints(gdm.selected.pred.45))
gdm.selected.pred[[3]]<-apply(cbind(as.data.frame(rasterToPoints(gdm.selected.pred.45))[[3]],as.data.frame(rasterToPoints(gdm.selected.pred.85))[[3]]),1,mean)
gdm.selected.pred[[3]]<-gdm.selected.pred[[3]]#*(Fstmax-Fstmin)
gdm.selected.pred<-rasterFromXYZ(gdm.selected.pred)

par(mfrow=c(2,2))
plot(ptremula,col="grey99",border="grey60",main="RCP4.5")
plot(gdm.selected.pred.45,col=cols1,axes=F,box=F,add=T,zlim=c(0,0.04))
plot(ptremula,col="grey99",border="grey60",main="RCP8.5")
plot(gdm.selected.pred.85,col=cols1,axes=F,box=F,add=T,zlim=c(0,0.04))

##Mahalanobis distance
mdist<-function(x1,x2,y1,y2,sy1,sy2){
  d<-sqrt(((x1-y1)/sy1)^2+((x2-y2)/sy2)^2)

  return(d)
}
sy1<-sd(r.SwAsp$climaticMoistureIndex@data@values,na.rm=T)
sy2<-sd(r.SwAsp$growingDegDays5@data@values,na.rm=T)

tmp1<-mdist(CG.env45[2,]$ccsm4_rcp45_2070_2.5arcmin_climaticMoistureIndex,CG.env45[2,]$ccsm4_rcp45_2070_2.5arcmin_growingDegDays5,r.SwAsp$climaticMoistureIndex,r.SwAsp$growingDegDays5,sy1,sy2)
tmp2<-mdist(CG.env85[2,]$ccsm4_rcp85_2070_2.5arcmin_climaticMoistureIndex,CG.env85[2,]$ccsm4_rcp85_2070_2.5arcmin_growingDegDays5,r.SwAsp$climaticMoistureIndex,r.SwAsp$growingDegDays5,sy1,sy2)

##Calculate distance from future climate at Sävar to current day climate
plot(ptremula,col="grey99",border="grey60",main="")
plot(tmp1,col=cols2,add=T,zlim=c(0,8))
points(CG[2,c(2,3)],pch=21,bg="grey60")
plot(ptremula,col="grey99",border="grey60",main="")
plot(tmp2,col=cols2,add=T,zlim=c(0,8))
points(CG[2,c(2,3)],pch=21,bg="grey60")

save.image(file="~/Dropbox/SwAsp/projects/climate_adaptation/analyses/disimilarity_modelling/disimilarity_modelling_new.RData")

