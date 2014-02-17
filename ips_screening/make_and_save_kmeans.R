library(stats)
library(gplots)
library(ggplot2)
library(sm)
library(class)

raw_data=read.table("/Users/b110-mm06/Desktop/Projects/ips_kinome_screen/data/All_Nuclei_ctrl.tab",sep="\t",dec=".",header=TRUE)
raw_data=raw_data[(which(!is.na(rowMeans(raw_data)))),]

names=as.vector(read.table("/Users/b110-mm06/Desktop/Projects/ips_kinome_screen/column_names.txt",sep="\t",dec=".",header=TRUE))


drop_v <- c(1,2,4,5,58,59,67,62,66,16,17,18,25,26,27,8,10)
names=colnames(names[,-c(drop_v)])
raw_data=raw_data[,-c(drop_v)]

colnames(raw_data)=names
correlation_matrix=cor(raw_data,method="spearman")

keep <- c(1,4,20,19,25,26,47,48)
raw_data=raw_data[,c(keep)]


#to_drop=list()
#for(i in 1:dim(correlation_matrix)[1]){
#  temp=correlation_matrix[correlation_matrix[,i]>0.90 ,i]
#  temp=temp[temp<1]
#  for(tempele in names(temp)){
#    to_drop[[tempele]]<-1
#  }
#}

#for(i in as.vector(names(to_drop))){
# if(!is.na(i)){
#   raw_data <- raw_data[,c(which(colnames(raw_data)!=i))]
#  }
#}

 print(raw_clust <-kmeans(raw_data,6, nstart = 30))
 
 pc=princomp(raw_data)
 
 
 
 panel.smooth <- function(x,y, ...)
{    smoothScatter(x,y,add=T)
}

save(raw_clust,file="kmeans_result.RData")
pairs(raw_data,panel=panel.smooth,col=raw_clust$cluster)

tiff(file='raw_data.tif')  
 pairs(raw_data,panel=panel.smooth,col=raw_clust$cluster)
  
# Write the file  
dev.off()  

svg(file='raw_data.svg', height=20, width=20, onefile=TRUE)  
  
# Plot your graph  
 pairs(raw_data,panel=panel.smooth,col=raw_clust$cluster)
  
# Write the file  
dev.off()  


svg(file='pca.svg', height=20, width=20, onefile=TRUE)  
  
# Plot your graph  
 pairs(pc$scores,panel=panel.smooth,col=raw_clust$cluster)
  
# Write the file  
dev.off()  


tiff(file='pca.tif')  

 pairs(pc$scores,panel=panel.smooth,col=raw_clust$cluster)
 
  
# Write the file  
dev.off()  
