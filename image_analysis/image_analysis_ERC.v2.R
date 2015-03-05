#!/usr/bin/R
library(EBImage)
library(FNN)
library(geometry)
args=commandArgs(trailingOnly = TRUE)
options("scipen"=100, "digits"=4,warn=-1)


dir=args[length(args)];


for(i in 1:(length(args)-1)){
imgages=unlist(strsplit(x=args[i],split="::"))

nuc_image_name=imgages[1] 
body_image_name=imgages[2] 
mito_image_name=imgages[3]

identifier=sub(".*/(.+)_\\w+.tif","\\1",nuc_image_name,perl=TRUE);

TEST_IMAGE_nuclei=TEST_IMAGE_nuclei_raw=readImage(nuc_image_name)#[500:1500,500:1500]
TEST_IMAGE_nuclei=TEST_IMAGE_nuclei/max(TEST_IMAGE_nuclei)
TEST_IMAGE_nuclei=gblur(1-(log(TEST_IMAGE_nuclei)/min(log(TEST_IMAGE_nuclei))),sigma=1)

TEST_IMAGE_actin=TEST_IMAGE_actin_raw=readImage(body_image_name)#[500:1500,500:1500]
TEST_IMAGE_actin=TEST_IMAGE_actin/max(TEST_IMAGE_actin)
TEST_IMAGE_actin=gblur(1-(log(TEST_IMAGE_actin)/min(log(TEST_IMAGE_actin))),sigma=1)

TEST_IMAGE_tubulin=TEST_IMAGE_tubulin_raw=readImage(mito_image_name)#[500:1500,500:1500]
TEST_IMAGE_tubulin=TEST_IMAGE_tubulin/max(TEST_IMAGE_tubulin)
TEST_IMAGE_tubulin=gblur(1-(log(TEST_IMAGE_tubulin)/min(log(TEST_IMAGE_tubulin))),sigma=1)


#display(TEST_IMAGE_nuclei)

if(max(bwlabel(thresh(TEST_IMAGE_nuclei,w=50,h=50,offset=0.05)))<1600){
  
  nuclei_binary=fillHull(dilate(erode(thresh(TEST_IMAGE_nuclei,w=31,h=31,offset=0.065),makeBrush(5))))
  actin_binary= TEST_IMAGE_actin 
  actin_binary [actin_binary > 0.2]  =1 #thresh(TEST_IMAGE_actin,w=80,h=80,offset=0.03)
  actin_binary [actin_binary < 1]  =0 
  tubulin_binary= TEST_IMAGE_tubulin
  tubulin_binary [tubulin_binary > 0.3]  =1 #thresh(TEST_IMAGE_fak,w=80,h=80,offset=0.03)
  tubulin_binary [tubulin_binary < 1]  =0 #tubulin_binary= thresh(TEST_IMAGE_tubulin,w=80,h=80,offset=0.07)
  nuclei_objects = bwlabel(nuclei_binary)  
  body_binary=dilate(erode(actin_binary+tubulin_binary))
  cell_bodies_objects= propagate(TEST_IMAGE_actin, seeds= nuclei_objects, mask=body_binary,lambda=3e-4)
  
  img = rgbImage(green=TEST_IMAGE_actin , blue=TEST_IMAGE_nuclei, red=TEST_IMAGE_tubulin)
  res = paintObjects(cell_bodies_objects, img, opac=c(0.3),thick=0.3,col='green')
  res = paintObjects(nuclei_objects, res, opac=c(0.3),thick=0.3, col='white')
 rm(body_binary,tubulin_binary,actin_binary,nuclei_binary)
  
}else{
  
  nuclei_binary=fillHull(dilate(erode(thresh(TEST_IMAGE_nuclei,w=21,h=21,offset=0.1),makeBrush(5))))
  actin_binary= TEST_IMAGE_actin 
  actin_binary [actin_binary > 0.3]  =1 #thresh(TEST_IMAGE_actin,w=80,h=80,offset=0.03)
  actin_binary [actin_binary < 1]  =0 
  tubulin_binary= TEST_IMAGE_tubulin
  tubulin_binary [tubulin_binary > 0.25]  =1 #thresh(TEST_IMAGE_fak,w=80,h=80,offset=0.03)
  tubulin_binary [tubulin_binary < 1]  =0 #tubulin_binary= thresh(TEST_IMAGE_tubulin,w=80,h=80,offset=0.07)
  nuclei_objects = bwlabel(nuclei_binary)  
  body_binary=dilate(erode(actin_binary+tubulin_binary))
  cell_bodies_objects= propagate(TEST_IMAGE_actin, seeds= nuclei_objects, mask=body_binary,lambda=3e-4)
  rm(body_binary,tubulin_binary,actin_binary,nuclei_binary)
  img = rgbImage(green=TEST_IMAGE_actin , blue=TEST_IMAGE_nuclei, red=TEST_IMAGE_tubulin)
  res = paintObjects(cell_bodies_objects, img, opac=c(0.3),thick=0.3,col='green')
  res = paintObjects(nuclei_objects, res, opac=c(0.3),thick=0.3, col='white')
 
}

rm(TEST_IMAGE_actin,TEST_IMAGE_nuclei,TEST_IMAGE_tubulin)



writeImage(res, file=paste(dir,"/",identifier,"_segmented",".tif",sep=""), type="tiff", quality = 100, 8 )

if(max(cell_bodies_objects)>10){
actin_features=cbind(
  (computeFeatures.basic(cell_bodies_objects,basic.quantiles=c(0),refnames="actin",ref=TEST_IMAGE_actin_raw)),
  (computeFeatures.shape(cell_bodies_objects,ref=TEST_IMAGE_actin_raw)),
  (computeFeatures.moment(cell_bodies_objects,ref=TEST_IMAGE_actin_raw)),
  (computeFeatures.haralick(cell_bodies_objects,ref=TEST_IMAGE_actin_raw,haralick.scales=c(1)))
)

DNA_features=cbind(
  (computeFeatures.basic(nuclei_objects,basic.quantiles=c(0),refnames="nuclei",ref=TEST_IMAGE_nuclei_raw)),
  (computeFeatures.shape(nuclei_objects,ref=TEST_IMAGE_nuclei_raw)),
  (computeFeatures.moment(nuclei_objects,ref=TEST_IMAGE_nuclei_raw)),
  (computeFeatures.haralick(nuclei_objects,ref=TEST_IMAGE_nuclei_raw,haralick.scales=c(1)))
)
tubulin_features=cbind(
  (computeFeatures.basic(cell_bodies_objects,basic.quantiles=c(0),refnames="tubulin",ref=TEST_IMAGE_tubulin_raw)),
  (computeFeatures.haralick(cell_bodies_objects,ref=TEST_IMAGE_tubulin_raw,haralick.scales=c(1)))
)

colnames(tubulin_features)=paste("tubulin",colnames(tubulin_features),sep=".")
colnames(actin_features)=paste("actin",colnames(actin_features),sep=".")
colnames(DNA_features)=paste("DNA",colnames(DNA_features),sep=".")
nearest.neighbours=get.knn(cbind(actin_features[,"actin.m.cx"],actin_features[,"actin.m.cy"]),k=10)[["nn.dist"]][,c(1,2,10)]
colnames(nearest.neighbours)=c("dist.1st.nn","dist.2n.nn","dist.10st.nn")

displacement<-function(x){
  return(sqrt((x[1]-x[3])^2+(x[2]-x[4])^2))
}
nucleus.displacement=apply(cbind(actin_features[,"actin.m.cx"],actin_features[,"actin.m.cy"],DNA_features[,"DNA.m.cx"],DNA_features[,"DNA.m.cy"]),1,displacement)
names(nucleus.displacement)=c("nuclear.displacement")

cell_features=cbind(tubulin_features,actin_features,DNA_features,nearest.neighbours,nucleus.displacement)



save(cell_features,file=paste(dir,"/",identifier,"_single_cell.RData",sep=""))

var.feat=apply(cell_features,2,var)
mean.feat=apply(cell_features,2,mean)
median.feat=apply(cell_features,2,median)
names(var.feat)=paste(names(var.feat),"var",sep=".")
names(mean.feat)=paste(names(mean.feat),"mean",sep=".")
names(median.feat)=paste(names(median.feat),"median",sep=".")

cell_features=c("cells"=max(cell_bodies_objects),mean.feat,median.feat,var.feat)

write.table(t(cell_features),file=paste(dir,"/",identifier,".tab",sep=""),sep="\t",quote=FALSE,row.names=FALSE)

}else{
    cell_features=rep(0,times=232)
    names(cell_features)=c("cells","tubulin.b.mean.var"
                           ,"tubulin.b.sd.var"
                           ,"tubulin.b.mad.var"
                           ,"tubulin.b.q0.var"
                           ,"tubulin.h.asm.s1.var"
                           ,"tubulin.h.con.s1.var"
                           ,"tubulin.h.cor.s1.var"
                           ,"tubulin.h.var.s1.var"
                           ,"tubulin.h.idm.s1.var"
                           ,"tubulin.h.sav.s1.var"
                           ,"tubulin.h.sva.s1.var"
                           ,"tubulin.h.sen.s1.var"
                           ,"tubulin.h.ent.s1.var"
                           ,"tubulin.h.dva.s1.var"
                           ,"tubulin.h.den.s1.var"
                           ,"tubulin.h.f12.s1.var"
                           ,"tubulin.h.f13.s1.var"
                           ,"actin.b.mean.var"
                           ,"actin.b.sd.var"
                           ,"actin.b.mad.var"
                           ,"actin.b.q0.var"
                           ,"actin.s.area.var"
                           ,"actin.s.perimeter.var"
                           ,"actin.s.radius.mean.var"
                           ,"actin.s.radius.sd.var"
                           ,"actin.s.radius.min.var"
                           ,"actin.s.radius.max.var"
                           ,"actin.m.cx.var"
                           ,"actin.m.cy.var"
                           ,"actin.m.majoraxis.var"
                           ,"actin.m.eccentricity.var"
                           ,"actin.m.theta.var"
                           ,"actin.h.asm.s1.var"
                           ,"actin.h.con.s1.var"
                           ,"actin.h.cor.s1.var"
                           ,"actin.h.var.s1.var"
                           ,"actin.h.idm.s1.var"
                           ,"actin.h.sav.s1.var"
                           ,"actin.h.sva.s1.var"
                           ,"actin.h.sen.s1.var"
                           ,"actin.h.ent.s1.var"
                           ,"actin.h.dva.s1.var"
                           ,"actin.h.den.s1.var"
                           ,"actin.h.f12.s1.var"
                           ,"actin.h.f13.s1.var"
                           ,"DNA.b.mean.var"
                           ,"DNA.b.sd.var"
                           ,"DNA.b.mad.var"
                           ,"DNA.b.q0.var"
                           ,"DNA.s.area.var"
                           ,"DNA.s.perimeter.var"
                           ,"DNA.s.radius.mean.var"
                           ,"DNA.s.radius.sd.var"
                           ,"DNA.s.radius.min.var"
                           ,"DNA.s.radius.max.var"
                           ,"DNA.m.cx.var"
                           ,"DNA.m.cy.var"
                           ,"DNA.m.majoraxis.var"
                           ,"DNA.m.eccentricity.var"
                           ,"DNA.m.theta.var"
                           ,"DNA.h.asm.s1.var"
                           ,"DNA.h.con.s1.var"
                           ,"DNA.h.cor.s1.var"
                           ,"DNA.h.var.s1.var"
                           ,"DNA.h.idm.s1.var"
                           ,"DNA.h.sav.s1.var"
                           ,"DNA.h.sva.s1.var"
                           ,"DNA.h.sen.s1.var"
                           ,"DNA.h.ent.s1.var"
                           ,"DNA.h.dva.s1.var"
                           ,"DNA.h.den.s1.var"
                           ,"DNA.h.f12.s1.var"
                           ,"DNA.h.f13.s1.var"
                           ,"dist.1st.nn.var"
                           ,"dist.2n.nn.var"
                           ,"dist.10st.nn.var"
                           ,"nucleus.displacement.var"
                           ,"tubulin.b.mean.mean"
                           ,"tubulin.b.sd.mean"
                           ,"tubulin.b.mad.mean"
                           ,"tubulin.b.q0.mean"
                           ,"tubulin.h.asm.s1.mean"
                           ,"tubulin.h.con.s1.mean"
                           ,"tubulin.h.cor.s1.mean"
                           ,"tubulin.h.var.s1.mean"
                           ,"tubulin.h.idm.s1.mean"
                           ,"tubulin.h.sav.s1.mean"
                           ,"tubulin.h.sva.s1.mean"
                           ,"tubulin.h.sen.s1.mean"
                           ,"tubulin.h.ent.s1.mean"
                           ,"tubulin.h.dva.s1.mean"
                           ,"tubulin.h.den.s1.mean"
                           ,"tubulin.h.f12.s1.mean"
                           ,"tubulin.h.f13.s1.mean"
                           ,"actin.b.mean.mean"
                           ,"actin.b.sd.mean"
                           ,"actin.b.mad.mean"
                           ,"actin.b.q0.mean"
                           ,"actin.s.area.mean"
                           ,"actin.s.perimeter.mean"
                           ,"actin.s.radius.mean.mean"
                           ,"actin.s.radius.sd.mean"
                           ,"actin.s.radius.min.mean"
                           ,"actin.s.radius.max.mean"
                           ,"actin.m.cx.mean"
                           ,"actin.m.cy.mean"
                           ,"actin.m.majoraxis.mean"
                           ,"actin.m.eccentricity.mean"
                           ,"actin.m.theta.mean"
                           ,"actin.h.asm.s1.mean"
                           ,"actin.h.con.s1.mean"
                           ,"actin.h.cor.s1.mean"
                           ,"actin.h.var.s1.mean"
                           ,"actin.h.idm.s1.mean"
                           ,"actin.h.sav.s1.mean"
                           ,"actin.h.sva.s1.mean"
                           ,"actin.h.sen.s1.mean"
                           ,"actin.h.ent.s1.mean"
                           ,"actin.h.dva.s1.mean"
                           ,"actin.h.den.s1.mean"
                           ,"actin.h.f12.s1.mean"
                           ,"actin.h.f13.s1.mean"
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
                           ,"DNA.h.asm.s1.mean"
                           ,"DNA.h.con.s1.mean"
                           ,"DNA.h.cor.s1.mean"
                           ,"DNA.h.var.s1.mean"
                           ,"DNA.h.idm.s1.mean"
                           ,"DNA.h.sav.s1.mean"
                           ,"DNA.h.sva.s1.mean"
                           ,"DNA.h.sen.s1.mean"
                           ,"DNA.h.ent.s1.mean"
                           ,"DNA.h.dva.s1.mean"
                           ,"DNA.h.den.s1.mean"
                           ,"DNA.h.f12.s1.mean"
                           ,"DNA.h.f13.s1.mean"
                           ,"dist.1st.nn.mean"
                           ,"dist.2n.nn.mean"
                           ,"dist.10st.nn.mean"
                           ,"nucleus.displacement.mean"
                           ,"tubulin.b.mean.median"
                           ,"tubulin.b.sd.median"
                           ,"tubulin.b.mad.median"
                           ,"tubulin.b.q0.median"
                           ,"tubulin.h.asm.s1.median"
                           ,"tubulin.h.con.s1.median"
                           ,"tubulin.h.cor.s1.median"
                           ,"tubulin.h.var.s1.median"
                           ,"tubulin.h.idm.s1.median"
                           ,"tubulin.h.sav.s1.median"
                           ,"tubulin.h.sva.s1.median"
                           ,"tubulin.h.sen.s1.median"
                           ,"tubulin.h.ent.s1.median"
                           ,"tubulin.h.dva.s1.median"
                           ,"tubulin.h.den.s1.median"
                           ,"tubulin.h.f12.s1.median"
                           ,"tubulin.h.f13.s1.median"
                           ,"actin.b.mean.median"
                           ,"actin.b.sd.median"
                           ,"actin.b.mad.median"
                           ,"actin.b.q0.median"
                           ,"actin.s.area.median"
                           ,"actin.s.perimeter.median"
                           ,"actin.s.radius.mean.median"
                           ,"actin.s.radius.sd.median"
                           ,"actin.s.radius.min.median"
                           ,"actin.s.radius.max.median"
                           ,"actin.m.cx.median"
                           ,"actin.m.cy.median"
                           ,"actin.m.majoraxis.median"
                           ,"actin.m.eccentricity.median"
                           ,"actin.m.theta.median"
                           ,"actin.h.asm.s1.median"
                           ,"actin.h.con.s1.median"
                           ,"actin.h.cor.s1.median"
                           ,"actin.h.var.s1.median"
                           ,"actin.h.idm.s1.median"
                           ,"actin.h.sav.s1.median"
                           ,"actin.h.sva.s1.median"
                           ,"actin.h.sen.s1.median"
                           ,"actin.h.ent.s1.median"
                           ,"actin.h.dva.s1.median"
                           ,"actin.h.den.s1.median"
                           ,"actin.h.f12.s1.median"
                           ,"actin.h.f13.s1.median"
                           ,"DNA.b.mean.median"
                           ,"DNA.b.sd.median"
                           ,"DNA.b.mad.median"
                           ,"DNA.b.q0.median"
                           ,"DNA.s.area.median"
                           ,"DNA.s.perimeter.median"
                           ,"DNA.s.radius.mean.median"
                           ,"DNA.s.radius.sd.median"
                           ,"DNA.s.radius.min.median"
                           ,"DNA.s.radius.max.median"
                           ,"DNA.m.cx.median"
                           ,"DNA.m.cy.median"
                           ,"DNA.m.majoraxis.median"
                           ,"DNA.m.eccentricity.median"
                           ,"DNA.m.theta.median"
                           ,"DNA.h.asm.s1.median"
                           ,"DNA.h.con.s1.median"
                           ,"DNA.h.cor.s1.median"
                           ,"DNA.h.var.s1.median"
                           ,"DNA.h.idm.s1.median"
                           ,"DNA.h.sav.s1.median"
                           ,"DNA.h.sva.s1.median"
                           ,"DNA.h.sen.s1.median"
                           ,"DNA.h.ent.s1.median"
                           ,"DNA.h.dva.s1.median"
                           ,"DNA.h.den.s1.median"
                           ,"DNA.h.f12.s1.median"
                           ,"DNA.h.f13.s1.median"
                           ,"dist.1st.nn.median"
                           ,"dist.2n.nn.median"
                           ,"dist.10st.nn.median"
                           ,"nucleus.displacement.median"
    )
    write.table(t(cell_features),file=paste(dir,"/",identifier,".tab",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
}
}