#!/usr/bin/R
library(EBImage)
args=commandArgs(trailingOnly = TRUE)
options(warn=-1)
nuc_image_name=args[1] 
body_image_name=args[2] 
mito_image_name=args[3]
dir=args[4]; 
identifier=args[5]; 

TEST_IMAGE_nuclei=TEST_IMAGE_nuclei_raw=readImage(nuc_image_name)[500:1500,500:1500]
TEST_IMAGE_nuclei=TEST_IMAGE_nuclei/max(TEST_IMAGE_nuclei)
TEST_IMAGE_nuclei=gblur(1-(log(TEST_IMAGE_nuclei)/min(log(TEST_IMAGE_nuclei))),sigma=1)

TEST_IMAGE_actin=TEST_IMAGE_actin_raw=readImage(body_image_name)[500:1500,500:1500]
TEST_IMAGE_actin=TEST_IMAGE_actin/max(TEST_IMAGE_actin)
TEST_IMAGE_actin=gblur(1-(log(TEST_IMAGE_actin)/min(log(TEST_IMAGE_actin))),sigma=1)

TEST_IMAGE_tubulin=TEST_IMAGE_tubulin_raw=readImage(mito_image_name)[500:1500,500:1500]
TEST_IMAGE_tubulin=TEST_IMAGE_tubulin/max(TEST_IMAGE_tubulin)
TEST_IMAGE_tubulin=gblur(1-(log(TEST_IMAGE_tubulin)/min(log(TEST_IMAGE_tubulin))),sigma=1)


#display(TEST_IMAGE_nuclei)

if(max(bwlabel(thresh(TEST_IMAGE_nuclei,w=50,h=50,offset=0.05)))<400){
  
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
    
  img = rgbImage(green=TEST_IMAGE_actin , blue=TEST_IMAGE_nuclei, red=TEST_IMAGE_tubulin)
  res = paintObjects(cell_bodies_objects, img, opac=c(0.3),thick=0.3,col='green')
  res = paintObjects(nuclei_objects, res, opac=c(0.3),thick=0.3, col='white')
 
}





writeImage(res, file=paste(dir,"/",identifier,"_segmented",".tif",sep=""), type="tiff", quality = 100, 8 )


actin_features=cbind(
  (computeFeatures.basic(cell_bodies_objects,basic.quantiles=c(0),refnames="actin",ref=TEST_IMAGE_actin_raw)),
  (computeFeatures.shape(cell_bodies_objects,ref=TEST_IMAGE_actin_raw)),
  (computeFeatures.moment(cell_bodies_objects,ref=TEST_IMAGE_actin_raw))
)

DNA_features=cbind(
  (computeFeatures.basic(nuclei_objects,basic.quantiles=c(0),refnames="nuclei",ref=TEST_IMAGE_nuclei_raw)),
  (computeFeatures.shape(nuclei_objects,ref=TEST_IMAGE_nuclei_raw)),
  (computeFeatures.moment(nuclei_objects,ref=TEST_IMAGE_nuclei_raw))
)
tubulin_features=cbind(
  (computeFeatures.basic(cell_bodies_objects,basic.quantiles=c(0),refnames="tubulin",ref=TEST_IMAGE_tubulin_raw)),
  (computeFeatures.haralick(cell_bodies_objects,ref=TEST_IMAGE_tubulin_raw))
)

colnames(tubulin_features)=paste("tubulin",colnames(tubulin_features),sep=".")
colnames(actin_features)=paste("actin",colnames(actin_features),sep=".")
colnames(DNA_features)=paste("DNA",colnames(DNA_features),sep=".")
cell_features=cbind(tubulin_features,actin_features,DNA_features)

save(cell_features,file=paste(dir,"/",identifier,"_single_cell.RData",sep=""))

cell_features=c("cells"=max(cell_bodies_objects),colMeans(cell_features))
write.table(t(cell_features),file=paste(dir,"/",identifier,".tab",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
