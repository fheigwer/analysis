library(class)
args=commandArgs(trailingOnly = TRUE)
options(warn=-1)
i=args[1]

sample=paste("/home/mount/fheigwer/ips_kinome_screen/data/Christian_IPSC/",i,"/DefaultOUT_Nuclei_ctrl.tab",sep="")

load(file="/home/mount/fheigwer/ips_kinome_screen/kmeans_result_non_rand.RData")
load(file="/home/mount/fheigwer/ips_kinome_screen/raw_controls_non_rand.RData")
sample_data=read.table(sample,sep="\t",dec=".",header=TRUE)
sample_data=sample_data[(which(!is.na(rowMeans(sample_data)))),]
	
origin=sample_data[,1]
names=as.vector(read.table("/home/mount/fheigwer/ips_kinome_screen/column_names.txt",sep="\t",dec=".",header=TRUE))
names=colnames(names)
colnames(sample_data)=names
sample_data=sample_data[,colnames(raw_clust$centers)]


for(y in 1:dim(raw_data)[2]){raw_data[,y]=(raw_data[,y]-min(raw_data[,y]))/(max(raw_data[,y])-min(raw_data[,y]))}
for(y in 1:dim(sample_data)[2]){sample_data[,y]=(sample_data[,y]-min(sample_data[,y]))/(max(sample_data[,y])-min(sample_data[,y]))}
save(sample_data,file=paste(i,"sample_data_ctrl_ctrltrain_non_rand.RData"))
raw_levels=knn(raw_data, sample_data, raw_clust$cluster, k = 100, prob=F)
clustering=cbind(origin,raw_levels)
save(clustering,file=paste(i,"clustering_ctrl_ctrltrain_non_rand.RData"))

