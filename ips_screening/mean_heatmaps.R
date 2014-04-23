library(scales)
library(scatterplot3d)
library(rgl)
library(pheatmap)

load('~/Desktop/Projects/ips_kinome_screen/10301 sample_data_ctrl_ctrltrain.RData')
load('~/Desktop/Projects/ips_kinome_screen/10301 clustering_ctrl_ctrltrain.RData')

a=read.table(file="~/Desktop/Projects/ips_kinome_screen/data/Christian_IPSC/10301/wellnames.txt")
plate_well_name=c()
for(i in 1:nrow(a)){
   plate_well_name=c(plate_well_name,paste(b[i,3])) #paste("10303",b[i,2],b[i,3],sep="_"))
}
pos=a[,1][which(a[,3]=="pos")]
neg=a[,1][which(a[,3]=="neg")]
colnames(sample_data)=c("Nuclear Area", "Nuclear Extent", "max. Oct4 Intensity", "max. Hoechst Intensity", "ø Hoechst Intensity", "ø Oct4 Intensity", "Distance from first neighbor", "Number of Neighbors")
		
poswell=c()
posmatrix=c(0,0,0,0,0,0,0,0)#c(0,0,0,0,0,0)
for(i in pos){
	poswell=c(poswell,which(clustering[,1]==i))
	postats=c()
	for(j in 1:6){
		postats=c(postats,length(which(clustering[which(clustering[,1]==i),2]==j)))
	}
	#postats=postats/sum(postats)*100
	postats=colMeans(sample_data[poswell,])
	posmatrix=rbind(posmatrix,postats)
}
colnames(posmatrix)=c("Nuclear Area", "Nuclear Extent", "max. Oct4 Intensity", "max. Hoechst Intensity", "ø Hoechst Intensity", "ø Oct4 Intensity", "Distance from first neighbor", "Number of Neighbors")
#colnames(posmatrix)=c("1","2","3","4","5","6")
posmatrix1=posmatrix[-1,]
rownames(posmatrix1)=plate_well_name[which(a[,3]=="pos")] #rep("pos",times=dim(posmatrix1)[1])

negwell=c()
negmatrix=c(0,0,0,0,0,0,0,0)#c(0,0,0,0,0,0)
for(i in neg){
	negwell=c(negwell,which(clustering[,1]==i))
	negtats=c()
	for(j in 1:6){
		negtats=c(negtats,length(which(clustering[which(clustering[,1]==i),2]==j)))
	}
	#negtats=negtats/sum(negtats)*100
	negtats=colMeans(sample_data[negwell,])
	negmatrix=rbind(negmatrix,negtats)
}
}
colnames(negmatrix)=c("Nuclear Area", "Nuclear Extent", "max. Oct4 Intensity", "max. Hoechst Intensity", "ø Hoechst Intensity", "ø Oct4 Intensity", "Distance from first neighbor", "Number of Neighbors")
#colnames(negmatrix)=c("1","2","3","4","5","6")
negmatrix1=negmatrix[-1,]
rownames(negmatrix1)=plate_well_name[which(a[,3]=="neg")]#rep("neg",times=dim(negmatrix1)[1])

load('~/Desktop/Projects/ips_kinome_screen/10302 sample_data_ctrl_ctrltrain.RData')
load('~/Desktop/Projects/ips_kinome_screen/10302 clustering_ctrl_ctrltrain.RData')
a=read.table(file="~/Desktop/Projects/ips_kinome_screen/data/Christian_IPSC/10302/wellnames.txt")
plate_well_name=c()
for(i in 1:nrow(a)){
   plate_well_name=c(plate_well_name,paste(b[i,3])) #paste("10303",b[i,2],b[i,3],sep="_"))
}
pos=a[,1][which(a[,3]=="pos")]
neg=a[,1][which(a[,3]=="neg")]
colnames(sample_data)=c("Nuclear Area", "Nuclear Extent", "max. Oct4 Intensity", "max. Hoechst Intensity", "ø Hoechst Intensity", "ø Oct4 Intensity", "Distance from first neighbor", "Number of Neighbors")
		
poswell=c()
posmatrix=c(0,0,0,0,0,0,0,0)#c(0,0,0,0,0,0)
for(i in pos){
	poswell=c(poswell,which(clustering[,1]==i))
	postats=c()
	for(j in 1:6){
		postats=c(postats,length(which(clustering[which(clustering[,1]==i),2]==j)))
	}
	#postats=postats/sum(postats)*100
	postats=colMeans(sample_data[poswell,])
	posmatrix=rbind(posmatrix,postats)
}
colnames(posmatrix)=c("Nuclear Area", "Nuclear Extent", "max. Oct4 Intensity", "max. Hoechst Intensity", "ø Hoechst Intensity", "ø Oct4 Intensity", "Distance from first neighbor", "Number of Neighbors")
#colnames(posmatrix)=c("1","2","3","4","5","6")
posmatrix2=posmatrix[-1,]
rownames(posmatrix2)=plate_well_name[which(a[,3]=="pos")]#rep("pos",times=dim(posmatrix2)[1])

negwell=c()
negmatrix=c(0,0,0,0,0,0,0,0)#c(0,0,0,0,0,0)
for(i in neg){
	negwell=c(negwell,which(clustering[,1]==i))
	negtats=c()
	for(j in 1:6){
		negtats=c(negtats,length(which(clustering[which(clustering[,1]==i),2]==j)))
	}
	#negtats=negtats/sum(negtats)*100
	negtats=colMeans(sample_data[negwell,])
	negmatrix=rbind(negmatrix,negtats)
}
colnames(negmatrix)=c("Nuclear Area", "Nuclear Extent", "max. Oct4 Intensity", "max. Hoechst Intensity", "ø Hoechst Intensity", "ø Oct4 Intensity", "Distance from first neighbor", "Number of Neighbors")
#colnames(negmatrix)=c("1","2","3","4","5","6")
negmatrix2=negmatrix[-1,]
rownames(negmatrix2)=plate_well_name[which(a[,3]=="neg")]#rep("neg",times=dim(negmatrix2)[1])

load('~/Desktop/Projects/ips_kinome_screen/10303 sample_data_ctrl_ctrltrain.RData')
load('~/Desktop/Projects/ips_kinome_screen/10303 clustering_ctrl_ctrltrain.RData')
a=read.table(file="~/Desktop/Projects/ips_kinome_screen/data/Christian_IPSC/10303/wellnames.txt")
plate_well_name=c()
for(i in 1:nrow(a)){
  plate_well_name=c(plate_well_name,paste(b[i,3])) #paste("10303",b[i,2],b[i,3],sep="_"))
}
pos=a[,1][which(a[,3]=="pos")]
neg=a[,1][which(a[,3]=="neg")]
colnames(sample_data)=c("Nuclear Area", "Nuclear Extent", "max. Oct4 Intensity", "max. Hoechst Intensity", "ø Hoechst Intensity", "ø Oct4 Intensity", "Distance from first neighbor", "Number of Neighbors")
		
poswell=c()
posmatrix=c(0,0,0,0,0,0,0,0)#c(0,0,0,0,0,0)
for(i in pos){
	poswell=c(poswell,which(clustering[,1]==i))
	postats=c()
	for(j in 1:6){
		postats=c(postats,length(which(clustering[which(clustering[,1]==i),2]==j)))
	}
	#postats=postats/sum(postats)*100
	postats=colMeans(sample_data[poswell,])
	posmatrix=rbind(posmatrix,postats)
}
colnames(posmatrix)=c("Nuclear Area", "Nuclear Extent", "max. Oct4 Intensity", "max. Hoechst Intensity", "ø Hoechst Intensity", "ø Oct4 Intensity", "Distance from first neighbor", "Number of Neighbors")
#colnames(posmatrix)=c("1","2","3","4","5","6")
posmatrix3=posmatrix[-1,]
rownames(posmatrix3)=plate_well_name[which(a[,3]=="pos")]#rep("pos",times=dim(posmatrix3)[1])

negwell=c()
negmatrix=c(0,0,0,0,0,0,0,0)#c(0,0,0,0,0,0)
for(i in neg){
	negwell=c(negwell,which(clustering[,1]==i))
	negtats=c()
	for(j in 1:6){
		negtats=c(negtats,length(which(clustering[which(clustering[,1]==i),2]==j)))
	}
	#negtats=negtats/sum(negtats)*100
	negtats=colMeans(sample_data[negwell,])
	negmatrix=rbind(negmatrix,negtats)
}
colnames(negmatrix)=c("Nuclear Area", "Nuclear Extent", "max. Oct4 Intensity", "max. Hoechst Intensity", "ø Hoechst Intensity", "ø Oct4 Intensity", "Distance from first neighbor", "Number of Neighbors")
#colnames(negmatrix)=c("1","2","3","4","5","6")
negmatrix3=negmatrix[-1,]
rownames(negmatrix3)=plate_well_name[which(a[,3]=="neg")]#rep("neg",times=dim(negmatrix3)[1])


load('~/Desktop/Projects/ips_kinome_screen/10301 sample_data_ctrltrain.RData')
load('~/Desktop/Projects/ips_kinome_screen/10301 clustering_ctrltrain.RData')
b=read.table(file="~/Desktop/Projects/ips_kinome_screen/data/Christian_IPSC/10301/wellnames_rna.txt")
plate_well_name=c()
for(i in 1:nrow(b)){
   plate_well_name=c(plate_well_name,paste(b[i,3])) #paste("10303",b[i,2],b[i,3],sep="_"))
}
colnames(sample_data)=c("Nuclear Area", "Nuclear Extent", "max. Oct4 Intensity", "max. Hoechst Intensity", "ø Hoechst Intensity", "ø Oct4 Intensity", "Distance from first neighbor", "Number of Neighbors")
sample=b[,1][which(b[,3]!="empty")]
samplewell=c()
samplematrix=c(0,0,0,0,0,0,0,0)#c(0,0,0,0,0,0)
for(i in sample){
	samplewell=c(samplewell,which(clustering[,1]==i))
	samplestats=c()
	for(j in 1:6){
		samplestats=c(samplestats,length(which(clustering[which(clustering[,1]==i),2]==j)))
	}
	#samplestats=samplestats/sum(samplestats)*100
	samplestats=colMeans(sample_data[samplewell,])
	samplematrix=rbind(samplematrix,samplestats)
}
colnames(samplematrix)=c("Nuclear Area", "Nuclear Extent", "max. Oct4 Intensity", "max. Hoechst Intensity", "ø Hoechst Intensity", "ø Oct4 Intensity", "Distance from first neighbor", "Number of Neighbors")
#colnames(samplematrix)=c("1","2","3","4","5","6")
samplematrix1=samplematrix[-1,]
rownames(samplematrix1)=plate_well_name[which(b[,3]!="empty")]

load('~/Desktop/Projects/ips_kinome_screen/10302 sample_data_ctrltrain.RData')
load('~/Desktop/Projects/ips_kinome_screen/10302 clustering_ctrltrain.RData')
b=read.table(file="~/Desktop/Projects/ips_kinome_screen/data/Christian_IPSC/10302/wellnames_rna.txt")
plate_well_name=c()
for(i in 1:nrow(b)){
    plate_well_name=c(plate_well_name,paste(b[i,3])) #paste("10303",b[i,2],b[i,3],sep="_"))
}
colnames(sample_data)=c("Nuclear Area", "Nuclear Extent", "max. Oct4 Intensity", "max. Hoechst Intensity", "ø Hoechst Intensity", "ø Oct4 Intensity", "Distance from first neighbor", "Number of Neighbors")
sample=b[,1][which(b[,3]!="empty")]
samplewell=c()
samplematrix=c(0,0,0,0,0,0,0,0)#c(0,0,0,0,0,0)
for(i in sample){
	samplewell=c(samplewell,which(clustering[,1]==i))
	samplestats=c()
	for(j in 1:6){
		samplestats=c(samplestats,length(which(clustering[which(clustering[,1]==i),2]==j)))
	}
	#samplestats=samplestats/sum(samplestats)*100
	samplestats=colMeans(sample_data[samplewell,])
	samplematrix=rbind(samplematrix,samplestats)
}
colnames(samplematrix)=c("Nuclear Area", "Nuclear Extent", "max. Oct4 Intensity", "max. Hoechst Intensity", "ø Hoechst Intensity", "ø Oct4 Intensity", "Distance from first neighbor", "Number of Neighbors")
#colnames(samplematrix)=c("1","2","3","4","5","6")
samplematrix2=samplematrix[-1,]
rownames(samplematrix2)=plate_well_name[which(b[,3]!="empty")]

load('~/Desktop/Projects/ips_kinome_screen/10303 sample_data_ctrltrain.RData')
load('~/Desktop/Projects/ips_kinome_screen/10303 clustering_ctrltrain.RData')
b=read.table(file="~/Desktop/Projects/ips_kinome_screen/data/Christian_IPSC/10303/wellnames_rna.txt")
plate_well_name=c()
for(i in 1:nrow(b)){
  plate_well_name=c(plate_well_name,paste(b[i,3])) #paste("10303",b[i,2],b[i,3],sep="_"))
}
colnames(sample_data)=c("Nuclear Area", "Nuclear Extent", "max. Oct4 Intensity", "max. Hoechst Intensity", "ø Hoechst Intensity", "ø Oct4 Intensity", "Distance from first neighbor", "Number of Neighbors")
sample=b[,1][which(b[,3]!="empty")]
samplewell=c()
samplematrix=c(0,0,0,0,0,0,0,0)#c(0,0,0,0,0,0)
for(i in sample){
	samplewell=c(samplewell,which(clustering[,1]==i))
	samplestats=c()
	for(j in 1:6){
		samplestats=c(samplestats,length(which(clustering[which(clustering[,1]==i),2]==j)))
	}
	#samplestats=samplestats/sum(samplestats)*100
	samplestats=colMeans(sample_data[samplewell,])
	samplematrix=rbind(samplematrix,samplestats)
}
colnames(samplematrix)=c("Nuclear Area", "Nuclear Extent", "max. Oct4 Intensity", "max. Hoechst Intensity", "ø Hoechst Intensity", "ø Oct4 Intensity", "Distance from first neighbor", "Number of Neighbors")
#colnames(samplematrix)=c("1","2","3","4","5","6")
samplematrix3=samplematrix[-1,]
rownames(samplematrix3)=plate_well_name[which(b[,3]!="empty")]


posmatrix=rbind(posmatrix1,posmatrix2,posmatrix3)
negmatrix=rbind(negmatrix1,negmatrix2,negmatrix3)
samplematrix=rbind(samplematrix1,samplematrix2,samplematrix3)

posMean=colMeans(posmatrix)
negMean=colMeans(negmatrix)
total=rbind(samplematrix,posmatrix,negmatrix)

dist_to_pos=c()
for(i in 1:dim(samplematrix)[1]){
	dist_to_pos=c(dist_to_pos,(sum((posMean-samplematrix[i,])^2))^(1/2))
}
names(dist_to_pos)=rownames(samplematrix)
dist_to_neg=c()
for(i in 1:dim(samplematrix)[1]){
	dist_to_neg=c(dist_to_neg,(sum((negMean-samplematrix[i,])^2))^(1/2))
}
names(dist_to_neg)=rownames(samplematrix)

plot_matrix=rbind(colMeans(posmatrix),colMeans(negmatrix),samplematrix[plot_names,])
row.names(plot_matrix)[1:2]=c("Oct4","non treated")
#svg(file="Heatmap_mean.svg")
#pheatmap(plot_matrix,cellheight=10)
#dev.off()
#dev.new()
norm_0_1=function(x){
	return((x-min(x))/(max(x)-min(x)))
	}
svg(file="Distplot_mean_values_103.svg")
plot(dist_to_pos, dist_to_neg,pch="*", xlab="lg ( euclidean distance to positive controls)", ylab="lg ( euclidean distance to negative controls)")
lines(c(log(median(dist_to_pos)), log(median(dist_to_pos))),c(0,max(dist_to_neg)),col="red")
lines(c(0, max(dist_to_pos)),c(log(median(dist_to_neg)), log(median(dist_to_neg))),col="red")
#points(log(dist_to_pos)[c(names(sort(dist_to_pos/dist_to_neg))[1:100])],log(dist_to_neg)[c(names(sort(dist_to_pos/dist_to_neg))[1:100])],col="green")
#points(log(dist_to_pos)[c("TAF1","SMG1")],log(dist_to_neg)[c("TAF1","SMG1")],col="red",pch="x")
text(c(2.25,2.25,4.25,4.25),c(0,3.75,0,4),c("Oct4-/+ like","Oct4- like","Oct4+ like","Indep. phenotype"))
dev.off()


realtive_hitlist_103=c(names(sort(dist_to_pos/dist_to_neg))[1:100])




#identify(log(dist_to_pos), log(dist_to_neg),labels=names(dist_to_pos),n=3)
