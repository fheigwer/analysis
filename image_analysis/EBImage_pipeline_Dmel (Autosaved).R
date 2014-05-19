#!/usr/bin/R
library(EBImage)
args=commandArgs(trailingOnly = TRUE)
options(warn=-1)
dir= "/Users/b110-mm06/Desktop/Projects/it_infrastructure/InCell_Test/4/20131223_HumanDmel_test_X038_BD19_2014.01.09.01.11.48" #args[1]
for(i in 1:24){
	for(j in c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P")){
		print(nuc_image_name)
nuc_image_name=paste(dir,"/",j,i,"fld1wvDAPIDAPI.tif",sep="");
body_image_name=paste(dir,"/",j,i,"fld1wvFITCFITC.tif",sep="");

identifier=paste(j,i,sep="");

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
#cell_bodies_objects[which(cell_bodies_objects>0)]=1 

cell_features=colMeans(computeFeatures.basic(cell_bodies_objects,basic.quantiles=c(0),refnames="actin",ref=TEST_IMAGE_body))
names(cell_features)=paste("alpha_tubulin",names(cell_features))
cell_features=c("#cells"=max(cell_bodies_objects),cell_features,colMeans(computeFeatures.shape(cell_bodies_objects,ref=TEST_IMAGE_body)),colMeans(computeFeatures.moment(cell_bodies_objects,ref=TEST_IMAGE_body)))
save(cell_features,file=paste(dir,"/",identifier,".RData",sep=""))

img = rgbImage(green=TEST_IMAGE_body, blue=TEST_IMAGE_nuclei)
res = paintObjects(cell_bodies_objects, img, opac=c(0.3),thick=0.3,col='green')
res = paintObjects(nuclei_objects, res, opac=c(0.3),thick=0.3, col='blue')
writeImage(res, file=paste(dir,"/",identifier,"_segmented",".tif",sep=""), type="tiff", quality = 100, 8 )
	}
}



