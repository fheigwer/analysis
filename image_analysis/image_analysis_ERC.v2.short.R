#!/usr/bin/R
library(EBImage)
library(FNN)
library(geometry)
args=commandArgs(trailingOnly = TRUE)
options("scipen"=100, "digits"=4,warn=-1)
dir=args[length(args)];


for(i in 1:(length(args)-1)){
nuc_image_name=args[i];
identifier=sub(".*/(.+)_\\w+.tif","\\1",args[i],perl=TRUE);

TEST_IMAGE_nuclei=TEST_IMAGE_nuclei_raw=readImage(nuc_image_name)#[500:1500,500:1500]
TEST_IMAGE_nuclei=TEST_IMAGE_nuclei/max(TEST_IMAGE_nuclei)
TEST_IMAGE_nuclei=gblur(1-(log(TEST_IMAGE_nuclei)/min(log(TEST_IMAGE_nuclei))),sigma=1)

  nuclei_binary=fillHull(dilate(erode(thresh(TEST_IMAGE_nuclei,w=21,h=21,offset=0.1),makeBrush(5))))
  nuclei_objects = bwlabel(nuclei_binary)  
  rm(nuclei_binary)
  res = rgbImage(blue=TEST_IMAGE_nuclei)
  res = paintObjects(nuclei_objects, res, opac=c(0.3),thick=0.3, col='white')
 


rm(TEST_IMAGE_nuclei)



writeImage(res, file=paste(dir,"/",identifier,"_segmented",".tif",sep=""), type="tiff", quality = 100, 8 )

if(max(nuclei_objects)>10){
DNA_features=cbind(
  (computeFeatures.basic(nuclei_objects,basic.quantiles=c(0),refnames="nuclei",ref=TEST_IMAGE_nuclei_raw)),
  (computeFeatures.shape(nuclei_objects,ref=TEST_IMAGE_nuclei_raw)),
  (computeFeatures.moment(nuclei_objects,ref=TEST_IMAGE_nuclei_raw))
)

colnames(DNA_features)=paste("DNA",colnames(DNA_features),sep=".")

mean.feat=apply(DNA_features,2,mean)
names(mean.feat)=paste(names(mean.feat),"mean",sep=".")

cell_features=c("cells"=max(nuclei_objects),mean.feat)

write.table(t(cell_features),file=paste(dir,"/",identifier,".qtab",sep=""),sep="\t",quote=FALSE,row.names=FALSE)

}else{
    cell_features=rep(0,times=232)
    names(cell_features)=c("cells"
                           ,"DNA.b.mean.mean"
                           ,"DNA.b.sd.mean"
                           ,"DNA.b.mad.mean"
                           ,"DNA.b.q0.mean"
                           ,"DNA.s.area.mean"
                           ,"DNA.s.perimeter.mean"
                           ,"DNA.s.radius.mean.mean"
                           ,"DNA.s.radius.sd.mean"
                           ,"DNA.s.radius.min.mean"
                           ,"DNA.s.radius.max.mean"
                           ,"DNA.m.cx.mean"
                           ,"DNA.m.cy.mean"
                           ,"DNA.m.majoraxis.mean"
                           ,"DNA.m.eccentricity.mean"
                           ,"DNA.m.theta.mean"
    )
    write.table(t(cell_features),file=paste(dir,"/",identifier,".qtab",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
}
}
