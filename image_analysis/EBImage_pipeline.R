#!/usr/bin/R
library(EBImage)
library(CRImage)
args=commandArgs(trailingOnly = TRUE)
options(warn=-1)
nuc_image_name=args[1]
body_image_name=args[2]
mito_image_name=args[3]
  
identifier=body_image_name#gsub("(..).+$","\\1",body_image_name)

TEST_IMAGE_nuclei=readImage(nuc_image_name)[500:1500, 500:1500]
TEST_IMAGE_nuclei=TEST_IMAGE_nuclei/max(TEST_IMAGE_nuclei)
TEST_IMAGE_nuclei=1-(log(TEST_IMAGE_nuclei)/min(log(TEST_IMAGE_nuclei)))


TEST_IMAGE_mito=readImage(mito_image_name)[500:1500, 500:1500]
TEST_IMAGE_mito=TEST_IMAGE_mito/max(TEST_IMAGE_mito)
TEST_IMAGE_mito=1-(log(TEST_IMAGE_mito)/min(log(TEST_IMAGE_mito)))

#filtermask= makeBrush(3, shape='gaussian')
nuclei_binary=thresh(TEST_IMAGE_nuclei,w=10,h=10,offset=0.09)
nuclei_binary=fillHull(nuclei_binary)
nuclei_objects = bwlabel(nuclei_binary)
nucleic_features=computeFeatures.shape(nuclei_objects)
cleaned_nuclei=rmObjects(nuclei_objects,which(log(nucleic_features[,1])<calculateOtsu(log(nucleic_features[,1]))))

TEST_IMAGE_body=readImage(body_image_name)[500:1500, 500:1500]
TEST_IMAGE_body=TEST_IMAGE_body/max(TEST_IMAGE_body)
TEST_IMAGE_body=1-(log(TEST_IMAGE_body)/min(log(TEST_IMAGE_body)))
cell_bodies_binary= fillHull(thresh(TEST_IMAGE_body,w=30,h=30,offset=0.02))
cell_bodies_objects= propagate(TEST_IMAGE_body, seeds=cleaned_nuclei, mask=cell_bodies_binary)

cell_features=colMeans(computeFeatures.basic(cell_bodies_objects,basic.quantiles=c(0),refnames="actin",ref=TEST_IMAGE_body))
names(cell_features)=paste("actin",names(cell_features))
mito_features=colMeans(computeFeatures.basic(cell_bodies_objects,basic.quantiles=c(0),refnames="tubulin",ref=TEST_IMAGE_mito))
names(mito_features)=paste("tubulin",names(mito_features))

cell_features=c("#cells"=max(cell_bodies_objects),mito_features,cell_features,colMeans(computeFeatures.shape(cell_bodies_objects,ref=TEST_IMAGE_body)),colMeans(computeFeatures.moment(cell_bodies_objects,ref=TEST_IMAGE_body)))
save(cell_features,file=paste(identifier,".RData",sep=""))

img = rgbImage(green=TEST_IMAGE_body, blue=TEST_IMAGE_nuclei,red=TEST_IMAGE_mito)
res = paintObjects(cell_bodies_objects, img, opac=c(0.5),col='green')
res = paintObjects(cleaned_nuclei, res, opac=c(0.5), col='blue')
save(res,file=paste(identifier,"_segmented",".RData",sep=""))
