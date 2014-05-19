library(stats)
library(gplots)
library(ggplot2)
library(sm)
library(class)

raw_data=read.table("/Users/b110-mm06/Desktop/Projects/ips_kinome_screen/data/Christian_IPSC/All_ctrl_non_rand.tab",sep="\t",dec=".",header=TRUE)
raw_data=raw_data[(which(!is.na(rowMeans(raw_data)))),]
names=as.vector(read.table("/Users/b110-mm06/Desktop/Projects/ips_kinome_screen/column_names.txt",sep="\t",dec=".",header=TRUE))
drop_v <- c(1,2,4,5,58,59,67,62,66,16,17,18,25,26,27,8,10)
names=colnames(names[,-c(drop_v)])
raw_data=raw_data[,-c(drop_v)]
colnames(raw_data)=names
keep <- c(1,4,20,19,25,26,47,48)
raw_data=raw_data[,c(keep)]

for(y in 1:dim(raw_data)[2]){raw_data[,y]=(raw_data[,y]-min(raw_data[,y]))/(max(raw_data[,y])-min(raw_data[,y]))}
raw_clust <-kmeans(raw_data,6, nstart = 30)
save(raw_data,file="raw_controls_non_rand.RData")
save(raw_clust,file="kmeans_result_non_rand.RData")



raw_data=read.table("/Users/b110-mm06/Desktop/Projects/ips_kinome_screen/data/Christian_IPSC/All_ctrl.tab",sep="\t",dec=".",header=TRUE)
raw_data=raw_data[(which(!is.na(rowMeans(raw_data)))),]
names=as.vector(read.table("/Users/b110-mm06/Desktop/Projects/ips_kinome_screen/column_names.txt",sep="\t",dec=".",header=TRUE))
drop_v <- c(1,2,4,5,58,59,67,62,66,16,17,18,25,26,27,8,10)
names=colnames(names[,-c(drop_v)])
raw_data=raw_data[,-c(drop_v)]
colnames(raw_data)=names
keep <- c(1,4,20,19,25,26,47,48)
raw_data=raw_data[,c(keep)]
for(y in 1:dim(raw_data)[2]){raw_data[,y]=(raw_data[,y]-min(raw_data[,y]))/(max(raw_data[,y])-min(raw_data[,y]))}
raw_clust <-kmeans(raw_data,6, nstart = 30)
save(raw_data,file="raw_controls.RData")
save(raw_clust,file="kmeans_result.RData")

