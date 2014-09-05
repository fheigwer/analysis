#!/usr/bin/R
library(EBImage)
args=commandArgs(trailingOnly = TRUE)
options(warn=-1)

nuc_image_name="/Users/b110-mm06/Desktop/IC-50-images/A01fld01wvDAPIDAPI.tif"		#args[1] 
body_image_name="/Users/b110-mm06/Desktop/IC-50-images/A01fld01wvCy3Cy3.tif"		#args[2]
tubulin= "/Users/b110-mm06/Desktop/IC-50-images/A01fld01wvFITCFITC.tif"				#args[3]
fak=	"/Users/b110-mm06/Desktop/IC-50-images/A01fld01wvCy5Cy5.tif"				#args[4]





dir="/Users/b110-mm06/Desktop/"	#args[4]; 
identifier="A01fld01"	#args[5]; 

TEST_IMAGE_nuclei=TEST_IMAGE_nuclei_raw=readImage(nuc_image_name)[500:1500,500:1500]
TEST_IMAGE_nuclei=TEST_IMAGE_nuclei/max(TEST_IMAGE_nuclei)
TEST_IMAGE_nuclei=gblur(1-(log(TEST_IMAGE_nuclei)/min(log(TEST_IMAGE_nuclei))),sigma=1)

TEST_IMAGE_actin=TEST_IMAGE_actin_raw=readImage(body_image_name)[500:1500,500:1500]
TEST_IMAGE_actin=TEST_IMAGE_actin/max(TEST_IMAGE_actin)
TEST_IMAGE_actin=gblur(1-(log(TEST_IMAGE_actin)/min(log(TEST_IMAGE_actin))),sigma=1)

TEST_IMAGE_tubulin=TEST_IMAGE_tubulin_raw=readImage(tubulin)[500:1500,500:1500]
TEST_IMAGE_tubulin=TEST_IMAGE_tubulin/max(TEST_IMAGE_tubulin)
TEST_IMAGE_tubulin=gblur(1-(log(TEST_IMAGE_tubulin)/min(log(TEST_IMAGE_tubulin))),sigma=1)

TEST_IMAGE_fak=TEST_IMAGE_fak_raw=readImage(fak)[500:1500,500:1500]
TEST_IMAGE_fak=TEST_IMAGE_fak/max(TEST_IMAGE_fak)
TEST_IMAGE_fak=gblur(1-(log(TEST_IMAGE_fak)/min(log(TEST_IMAGE_fak))),sigma=1)

nuclei_binary=thresh(TEST_IMAGE_nuclei,w=50,h=50,offset=0.05)

actin_binary= dilate(thresh(TEST_IMAGE_actin,w=50,h=50,offset=0.02))
tubulin_binary= dilate(thresh(TEST_IMAGE_tubulin,w=50,h=50,offset=0.02))

body_binary=actin_binary+tubulin_binary

nuclei_binary=nuclei_binary*body_binary

nuclei_objects = bwlabel(nuclei_binary)

cell_bodies_objects= propagate(TEST_IMAGE_actin+TEST_IMAGE_tubulin, seeds= nuclei_objects, mask=body_binary)

img = rgbImage(green=TEST_IMAGE_actin,red=TEST_IMAGE_tubulin, blue=TEST_IMAGE_nuclei)
res = paintObjects(cell_bodies_objects, img, opac=c(0.3),thick=0.3,col='green')
res = paintObjects(nuclei_objects, res, opac=c(0.3),thick=0.3, col='blue')
display(res)


writeImage(res, file=paste(dir,"/",identifier,"_segmented",".tif",sep=""), type="tiff", quality = 100, 8 )

actin_features=c(
		colMeans(computeFeatures.basic(cell_bodies_objects,basic.quantiles=c(0),refnames="actin",ref=TEST_IMAGE_actin_raw)),
		colMeans(computeFeatures.shape(cell_bodies_objects,ref=TEST_IMAGE_actin_raw)),
		colMeans(computeFeatures.moment(cell_bodies_objects,ref=TEST_IMAGE_actin_raw))
		)

DNA_features=c(
		colMeans(computeFeatures.basic(nuclei_objects,basic.quantiles=c(0),refnames="nuclei",ref=TEST_IMAGE_nuclei_raw)),
		colMeans(computeFeatures.shape(nuclei_objects,ref=TEST_IMAGE_nuclei_raw)),
		colMeans(computeFeatures.moment(nuclei_objects,ref=TEST_IMAGE_nuclei_raw))
		)
tubulin_features=c(
		colMeans(computeFeatures.basic(cell_bodies_objects,basic.quantiles=c(0),refnames="tubulin",ref=TEST_IMAGE_tubulin_raw)),
		colMeans(computeFeatures.shape(cell_bodies_objects,ref=TEST_IMAGE_tubulin_raw)),
		colMeans(computeFeatures.moment(cell_bodies_objects,ref=TEST_IMAGE_tubulin_raw))
		)	
FAK_features=c(
		colMeans(computeFeatures.basic(cell_bodies_objects,basic.quantiles=c(0),refnames="FAK",ref=TEST_IMAGE_fak_raw)),
		colMeans(computeFeatures.shape(cell_bodies_objects,ref=TEST_IMAGE_fak_raw)),
		colMeans(computeFeatures.moment(cell_bodies_objects,ref=TEST_IMAGE_fak_raw)),
		colMeans(computeFeatures.haralick(cell_bodies_objects,ref=TEST_IMAGE_fak_raw))
				)				
names(tubulin_features)=paste("tubulin",names(tubulin_features),sep=".")
names(actin_features)=paste("actin",names(actin_features),sep=".")
names(FAK_features)=paste("FAK",names(FAK_features),sep=".")
names(DNA_features)=paste("DNA",names(DNA_features),sep=".")


cell_features=c("#cells"=max(cell_bodies_objects),tubulin_features,actin_features,FAK_features,DNA_features)
write.table(t(cell_features),file=paste(dir,"/",identifier,".tab",sep=""),sep="\t",quote=FALSE, row.names=FALSE)