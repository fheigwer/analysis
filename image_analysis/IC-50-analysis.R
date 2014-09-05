#!/usr/bin/R
library(EBImage)
args=commandArgs(trailingOnly = TRUE)
options(warn=-1)

for(i in c("I","J","K","L","M","N","O","P")){
  for(j in 1:24){
    for(k in 1:4){
      if(j<10){
        nuc_image_name=  paste("/Users/b110-mm06/Desktop/IC-50-images/plate1/",i,"0",j,"fld0",k,"wvDAPIDAPI.tif"   ,sep=""  )  #args[1] 
        body_image_name=  paste("/Users/b110-mm06/Desktop/IC-50-images/plate1/",i,"0",j,"fld0",k,"wvCy3Cy3.tif" ,sep=""  )		#args[2]
        tubulin=   paste("/Users/b110-mm06/Desktop/IC-50-images/plate1/",i,"0",j,"fld0",k,"wvFITCFITC.tif"	 ,sep=""  )			#args[3]
        fak=	paste("/Users/b110-mm06/Desktop/IC-50-images/plate1/",i,"0",j,"fld0",k,"wvCy5Cy5.tif" ,sep=""  )				#args[4]
        dir="/Users/b110-mm06/Desktop/IC-50-results/plate_1/"	#args[4]; 
        identifier=  paste("",i,"0",j,"fld0",k,sep="")  #args[5]; 
      }else{
        nuc_image_name=  paste("/Users/b110-mm06/Desktop/IC-50-images/plate1/",i,j,"fld0",k,"wvDAPIDAPI.tif" ,sep=""  ) #args[1] 
        body_image_name=  paste("/Users/b110-mm06/Desktop/IC-50-images/plate1/",i,j,"fld0",k,"wvCy3Cy3.tif",sep="" ) 	#args[2]
        tubulin=   paste("/Users/b110-mm06/Desktop/IC-50-images/plate1/",i,j,"fld0",k,"wvFITCFITC.tif"		,sep="")		#args[3]
        fak=	paste("/Users/b110-mm06/Desktop/IC-50-images/plate1/",i,j,"fld0",k,"wvCy5Cy5.tif",	sep="")			#args[4]
        dir="/Users/b110-mm06/Desktop/IC-50-results/plate_1/"	#args[4]; 
        identifier=  paste("",i,j,"fld0",k,sep="")	#args[5]; 
      }


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

#display(TEST_IMAGE_nuclei)

if(max(bwlabel(thresh(TEST_IMAGE_nuclei,w=50,h=50,offset=0.05)))<400){
  
  nuclei_binary=fillHull(dilate(erode(thresh(TEST_IMAGE_nuclei,w=31,h=31,offset=0.065),makeBrush(5))))
  actin_binary= TEST_IMAGE_actin 
  actin_binary [actin_binary > 0.2]  =1 #thresh(TEST_IMAGE_actin,w=80,h=80,offset=0.03)
  actin_binary [actin_binary < 1]  =0 
  #fak_binary= TEST_IMAGE_fak 
  #fak_binary [fak_binary > 0.2]  =1 #thresh(TEST_IMAGE_fak,w=80,h=80,offset=0.03)
  #fak_binary [fak_binary < 1]  =0 #tubulin_binary= thresh(TEST_IMAGE_tubulin,w=80,h=80,offset=0.07)
  fak_binary= fillHull(dilate(erode(thresh(TEST_IMAGE_fak,w=60,h=60,offset=0),makeBrush(5))))
  tubulin_binary= TEST_IMAGE_tubulin
  tubulin_binary [tubulin_binary > 0.3]  =1 #thresh(TEST_IMAGE_fak,w=80,h=80,offset=0.03)
  tubulin_binary [tubulin_binary < 1]  =0 #tubulin_binary= thresh(TEST_IMAGE_tubulin,w=80,h=80,offset=0.07)
  nuclei_objects = bwlabel(nuclei_binary)  
  body_binary=dilate(erode(actin_binary+fak_binary+tubulin_binary))
  cell_bodies_objects= propagate(TEST_IMAGE_actin, seeds= nuclei_objects, mask=body_binary,lambda=3e-4)
  
  img = rgbImage(green=TEST_IMAGE_actin , blue=TEST_IMAGE_nuclei)
  res = paintObjects(cell_bodies_objects, img, opac=c(0.3),thick=0.3,col='green')
  res = paintObjects(nuclei_objects, img, opac=c(0.3),thick=0.3, col='white')
  res = paintObjects(cell_bodies_objects, res, opac=c(0.3),thick=0.3,col='cyan')
 
  
}else{
  
  nuclei_binary=fillHull(dilate(erode(thresh(TEST_IMAGE_nuclei,w=21,h=21,offset=0.1),makeBrush(5))))
  actin_binary= TEST_IMAGE_actin 
  actin_binary [actin_binary > 0.3]  =1 #thresh(TEST_IMAGE_actin,w=80,h=80,offset=0.03)
  actin_binary [actin_binary < 1]  =0 
  #fak_binary= TEST_IMAGE_fak 
  #fak_binary [fak_binary > 0.25]  =1 #thresh(TEST_IMAGE_fak,w=80,h=80,offset=0.03)
  #fak_binary [fak_binary < 1]  =0 #tubulin_binary= thresh(TEST_IMAGE_tubulin,w=80,h=80,offset=0.07)
  fak_binary= fillHull(dilate(erode(thresh(TEST_IMAGE_fak,w=60,h=60,offset=0),makeBrush(5))))
  tubulin_binary= TEST_IMAGE_tubulin
  tubulin_binary [tubulin_binary > 0.25]  =1 #thresh(TEST_IMAGE_fak,w=80,h=80,offset=0.03)
  tubulin_binary [tubulin_binary < 1]  =0 #tubulin_binary= thresh(TEST_IMAGE_tubulin,w=80,h=80,offset=0.07)
  nuclei_objects = bwlabel(nuclei_binary)  
  body_binary=dilate(erode(actin_binary+fak_binary+tubulin_binary))
  cell_bodies_objects= propagate(TEST_IMAGE_actin, seeds= nuclei_objects, mask=body_binary,lambda=3e-4)
  
  
  img = rgbImage(green=TEST_IMAGE_actin , blue=TEST_IMAGE_nuclei)
  res = paintObjects(cell_bodies_objects, img, opac=c(0.3),thick=0.3,col='green')
  res = paintObjects(nuclei_objects, img, opac=c(0.3),thick=0.3, col='white')
  res = paintObjects(cell_bodies_objects, res, opac=c(0.3),thick=0.3,col='cyan')
 
}





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
  colMeans(computeFeatures.haralick(cell_bodies_objects,ref=TEST_IMAGE_tubulin_raw))
)	
FAK_features=c(
  colMeans(computeFeatures.basic(cell_bodies_objects,basic.quantiles=c(0),refnames="FAK",ref=TEST_IMAGE_fak_raw)),
  colMeans(computeFeatures.haralick(cell_bodies_objects,ref=TEST_IMAGE_fak_raw))
)				
names(tubulin_features)=paste("tubulin",names(tubulin_features),sep=".")
names(actin_features)=paste("actin",names(actin_features),sep=".")
names(FAK_features)=paste("FAK",names(FAK_features),sep=".")
names(DNA_features)=paste("DNA",names(DNA_features),sep=".")


cell_features=c("cells"=max(cell_bodies_objects),tubulin_features,actin_features,FAK_features,DNA_features)
write.table(t(cell_features),file=paste(dir,"/",identifier,".tab",sep=""),sep="\t",quote=FALSE, row.names=identifier)
    }
  }
}
