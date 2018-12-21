library("data.table")
source("~/Dropbox/RScripts/manhattan_plot.R")

addalpha <- function(colors, alpha=1.0) {
  r <- col2rgb(colors, alpha=T)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}

##hapFLK results
hfiles<-list.files("~/Dropbox/SwAsp/hapFLK/results",pattern="*_fit15.hapflk",full.names = TRUE)
hapflk<-NULL
for(i in hfiles){
  hf<-fread(i,head=T)
  hapflk<-rbind(hapflk,hf)
}
##Sort chromosome order
hapflk<-setorder(hapflk,chr,pos)

##Read LFMM results
lfmm<-fread("~/Dropbox/SwAsp/projects/climate_adaptation/analyses/lfmm/lfmm_analyses.txt")

##FLK results
hfiles<-list.files("~/Dropbox/SwAsp/hapFLK/results",pattern="*_fit15.flk",full.names = TRUE)
flk<-NULL
for(i in hfiles){
  hf<-fread(i,head=T)
  flk<-rbind(flk,hf)
}
rm(hf)
##Sort chromosome order
flk<-setorder(flk,chr,pos)

pdf("~/Dropbox/SwAsp/projects/climate_adaptation/figures and tables/pdf/hapFLK.pdf",height=6,width=14)
col1<-addalpha(c("lightsteelblue3"),0.7)
col2<-addalpha(c("lightsteelblue1"),0.7)
idx<-(hapflk$hapflk>1)
#idx2<-hapflk$hapflk>4.414
idx2<-hapflk$hapflk>3
manhattan(hapflk[idx,],chr="chr",bp="pos",p="hapflk",logp=F,ymin=1,ymax=12,col=c(col1,col2),cex=0.75,highlight = idx2[idx],suggestivecolor = "red")
dev.off()

##Plot cumulative hapFLK for all and selected loci
ehapflk<-ecdf(hapflk$hapflk)
idx1<-(lfmm$p.cmi<1e-6)
idx2<-(lfmm$p.cont<1e-6)
idx3<-(lfmm$p.cont<1e-6&lfmm$CHR!=10)
ehapflk1<-ecdf(hapflk$hapflk[idx1])
ehapflk2<-ecdf(hapflk$hapflk[idx2])
ehapflk3<-ecdf(hapflk$hapflk[idx3])
pdf("~/Dropbox/SwAsp/projects/climate_adaptation/figures and tables/pdf/cumulative_hapflk.pdf",height=7,width=7)
plot(ehapflk,main="",xlab="hapFLK",ylab="Cumulative proportion", las=1,col="lightsteelblue3",lwd=2,las=1)
lines(ehapflk1,col="indianred4",lwd=2,do.points=F,verticals=T)
lines(ehapflk2,col="indianred2",lwd=2,do.points=F,verticals=T)
lines(ehapflk3,col="indianred2",lwd=2,lty=6,do.points=F,verticals=T)
dev.off()
