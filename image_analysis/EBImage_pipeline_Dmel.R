#!/usr/bin/R
library(EBImage)
args=commandArgs(trailingOnly = TRUE)
options(warn=-1)

nuc_image_name=args[1] 
body_image_name=args[2] 
dir=args[3]; 
identifier=args[4]; 

TEST_IMAGE_nuclei=readImage(nuc_image_name)[500:1500,500:1500]
TEST_IMAGE_nuclei=TEST_IMAGE_nuclei/max(TEST_IMAGE_nuclei)
TEST_IMAGE_nuclei=gblur(1-(log(TEST_IMAGE_nuclei)/min(log(TEST_IMAGE_nuclei))),sigma=1)

TEST_IMAGE_body=readImage(body_image_name)[500:1500,500:1500]
TEST_IMAGE_body=TEST_IMAGE_body/max(TEST_IMAGE_body)
TEST_IMAGE_body=gblur(1-(log(TEST_IMAGE_body)/min(log(TEST_IMAGE_body))),sigma=1)

nuclei_binary=thresh(TEST_IMAGE_nuclei,w=5,h=5,offset=0.1)
cell_bodies_binary= dilate(thresh(TEST_IMAGE_body,w=5,h=5,offset=0.05))
nuclei_binary=nuclei_binary*cell_bodies_binary

nuclei_objects = bwlabel(nuclei_binary)
cell_bodies_objects= propagate(TEST_IMAGE_body, seeds= nuclei_objects, mask=cell_bodies_binary)


cell_features=colMeans(computeFeatures.basic(cell_bodies_objects,basic.quantiles=c(0),refnames="actin",ref=TEST_IMAGE_body))
names(cell_features)=paste("alpha_tubulin",names(cell_features))
cell_features=c("#cells"=max(cell_bodies_objects),cell_features,colMeans(computeFeatures.shape(cell_bodies_objects,ref=TEST_IMAGE_body)),colMeans(computeFeatures.moment(cell_bodies_objects,ref=TEST_IMAGE_body)))
write.table(t(cell_features),file=paste(dir,"/",identifier,".tab",sep=""),sep="\t",quote=FALSE, row.names=FALSE)

img = rgbImage(green=TEST_IMAGE_body, blue=TEST_IMAGE_nuclei)
res = paintObjects(cell_bodies_objects, img, opac=c(0.3),thick=0.3,col='green')
res = paintObjects(nuclei_objects, res, opac=c(0.3),thick=0.3, col='blue')
writeImage(res, file=paste(dir,"/",identifier,"_segmented",".tif",sep=""), type="tiff", quality = 100, 8 )




