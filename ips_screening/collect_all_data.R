
sample_data_complete=c()
clustering_complete=c()
names_complete=c()

for(i in c("10301","10302","10303","10401","10402","10403")){
  
  load(paste('~/Desktop/Projects/ips_kinome_screen/',i,' sample_data_ctrltrain.RData',sep=""))
  assign(paste("raw_data_sample_", i, sep=""),sample_data)
  sample_data_complete=rbind(sample_data_complete,sample_data)
  save(sample_data,file=paste("raw_data_sample_", i,".RData", sep=""))
  
  load(paste('~/Desktop/Projects/ips_kinome_screen/',i,' clustering_ctrltrain.RData',sep=""))
  assign(paste("clustering_sample_", i, sep=""),clustering)
  save(clustering,file=paste("clustering_sample_", i,".RData", sep=""))
  clustering_complete=rbind(clustering_complete,clustering)
  
  a=as.data.frame(read.table(file=paste("~/Desktop/Projects/ips_kinome_screen/data/Christian_IPSC/",i,"/wellnames.txt",sep="")))
  b=as.data.frame(read.table(file=paste("~/Desktop/Projects/ips_kinome_screen/data/Christian_IPSC/",i,"/wellnames_rna.txt",sep="")))
  c=paste(i,as.matrix(b[,3]),sep="_")
  c[grep("pos",a[,3])]<-"pos"
  c[grep("neg",a[,3])]<-"neg"
  b[,3]<-c
  save(b,file=paste("wellnames",i,".RData",sep=""))
  names_complete=rbind(names_complete,b)
  
  sample=paste("~/Desktop/Projects/ips_kinome_screen/data/Christian_IPSC/",i,"/DefaultOUT_Nuclei_ctrl.tab",sep="")
  sample_data=read.table(sample,sep="\t",dec=".",header=TRUE)
  origin=sample_data[(which(!is.na(rowMeans(sample_data)))),1]
  save(origin,file=paste("origin",i,".RData",sep=""))  
}

save(clustering_complete,file="all_samples_classifications.RData")
save(sample_data_complete,file="all_samples_feature_data.RData")
save(names_complete,file="all_samples_annotation.RData")
