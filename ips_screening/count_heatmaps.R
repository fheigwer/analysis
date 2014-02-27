library(scales)
library(scatterplot3d)
library(rgl)
library(pheatmap)

load('~/Desktop/Projects/ips_kinome_screen/10401 sample_data_ctrl_ctrltrain.RData')
load('~/Desktop/Projects/ips_kinome_screen/10401 clustering_ctrl_ctrltrain.RData')
a=read.table(file="~/Desktop/Projects/ips_kinome_screen/data/Christian_IPSC/10401/wellnames.txt")
pos=a[,1][which(a[,3]=="pos")]
neg=a[,1][which(a[,3]=="neg")]
colnames(sample_data)=c("Nuclear Area", "Nuclear Extent", "max. Oct4 Intensity", "max. Hoechst Intensity", "ø Hoechst Intensity", "ø Oct4 Intensity", "Distance from first neighbor", "Number of Neighbors")
		
poswell=c()
posmatrix=c(0,0,0,0,0,0)
for(i in pos){
	poswell=c(poswell,which(clustering[,1]==i))
	postats=c()
	for(j in 1:6){
		postats=c(postats,length(which(clustering[which(clustering[,1]==i),2]==j)))
	}
	postats=postats/sum(postats)*100
	#postats=colMeans(sample_data[poswell,])
	posmatrix=rbind(posmatrix,postats)
}
#colnames(posmatrix)=c("1","2","3","4","5","6")
posmatrix1=posmatrix[-1,]
rownames(posmatrix1)=rep("pos",times=dim(posmatrix1)[1])

negwell=c()
negmatrix=c(0,0,0,0,0,0)
for(i in neg){
	negwell=c(negwell,which(clustering[,1]==i))
	negtats=c()
	for(j in 1:6){
		negtats=c(negtats,length(which(clustering[which(clustering[,1]==i),2]==j)))
	}
	negtats=negtats/sum(negtats)*100
	#negtats=colMeans(sample_data[negwell,])
	negmatrix=rbind(negmatrix,negtats)
}
colnames(negmatrix)=c("1","2","3","4","5","6")
negmatrix1=negmatrix[-1,]
rownames(negmatrix1)=rep("neg",times=dim(negmatrix1)[1])

load('~/Desktop/Projects/ips_kinome_screen/10402 sample_data_ctrl_ctrltrain.RData')
load('~/Desktop/Projects/ips_kinome_screen/10402 clustering_ctrl_ctrltrain.RData')
a=read.table(file="~/Desktop/Projects/ips_kinome_screen/data/Christian_IPSC/10402/wellnames.txt")
pos=a[,1][which(a[,3]=="pos")]
neg=a[,1][which(a[,3]=="neg")]
colnames(sample_data)=c("Nuclear Area", "Nuclear Extent", "max. Oct4 Intensity", "max. Hoechst Intensity", "ø Hoechst Intensity", "ø Oct4 Intensity", "Distance from first neighbor", "Number of Neighbors")
		
poswell=c()
posmatrix=c(0,0,0,0,0,0)
for(i in pos){
	poswell=c(poswell,which(clustering[,1]==i))
	postats=c()
	for(j in 1:6){
		postats=c(postats,length(which(clustering[which(clustering[,1]==i),2]==j)))
	}
	postats=postats/sum(postats)*100
	#postats=colMeans(sample_data[poswell,])
	posmatrix=rbind(posmatrix,postats)
}
#colnames(posmatrix)=c("1","2","3","4","5","6")
posmatrix2=posmatrix[-1,]
rownames(posmatrix2)=rep("pos",times=dim(posmatrix2)[1])

negwell=c()
negmatrix=c(0,0,0,0,0,0)
for(i in neg){
	negwell=c(negwell,which(clustering[,1]==i))
	negtats=c()
	for(j in 1:6){
		negtats=c(negtats,length(which(clustering[which(clustering[,1]==i),2]==j)))
	}
	negtats=negtats/sum(negtats)*100
	#negtats=colMeans(sample_data[negwell,])
	negmatrix=rbind(negmatrix,negtats)
}
colnames(negmatrix)=c("1","2","3","4","5","6")
negmatrix2=negmatrix[-1,]
rownames(negmatrix2)=rep("neg",times=dim(negmatrix2)[1])

load('~/Desktop/Projects/ips_kinome_screen/10403 sample_data_ctrl_ctrltrain.RData')
load('~/Desktop/Projects/ips_kinome_screen/10403 clustering_ctrl_ctrltrain.RData')
a=read.table(file="~/Desktop/Projects/ips_kinome_screen/data/Christian_IPSC/10403/wellnames.txt")
pos=a[,1][which(a[,3]=="pos")]
neg=a[,1][which(a[,3]=="neg")]
colnames(sample_data)=c("Nuclear Area", "Nuclear Extent", "max. Oct4 Intensity", "max. Hoechst Intensity", "ø Hoechst Intensity", "ø Oct4 Intensity", "Distance from first neighbor", "Number of Neighbors")
		
poswell=c()
posmatrix=c(0,0,0,0,0,0)
for(i in pos){
	poswell=c(poswell,which(clustering[,1]==i))
	postats=c()
	for(j in 1:6){
		postats=c(postats,length(which(clustering[which(clustering[,1]==i),2]==j)))
	}
	postats=postats/sum(postats)*100
	#postats=colMeans(sample_data[poswell,])
	posmatrix=rbind(posmatrix,postats)
}
#colnames(posmatrix)=c("1","2","3","4","5","6")
posmatrix3=posmatrix[-1,]
rownames(posmatrix3)=rep("pos",times=dim(posmatrix3)[1])

negwell=c()
negmatrix=c(0,0,0,0,0,0)
for(i in neg){
	negwell=c(negwell,which(clustering[,1]==i))
	negtats=c()
	for(j in 1:6){
		negtats=c(negtats,length(which(clustering[which(clustering[,1]==i),2]==j)))
	}
	negtats=negtats/sum(negtats)*100
	#negtats=colMeans(sample_data[negwell,])
	negmatrix=rbind(negmatrix,negtats)
}
colnames(negmatrix)=c("1","2","3","4","5","6")
negmatrix3=negmatrix[-1,]
rownames(negmatrix3)=rep("neg",times=dim(negmatrix3)[1])


load('~/Desktop/Projects/ips_kinome_screen/10401 sample_data_ctrltrain.RData')
load('~/Desktop/Projects/ips_kinome_screen/10401 clustering_ctrltrain.RData')
b=read.table(file="~/Desktop/Projects/ips_kinome_screen/data/Christian_IPSC/10401/wellnames_rna.txt")
colnames(sample_data)=c("Nuclear Area", "Nuclear Extent", "max. Oct4 Intensity", "max. Hoechst Intensity", "ø Hoechst Intensity", "ø Oct4 Intensity", "Distance from first neighbor", "Number of Neighbors")
sample=b[,1][which(b[,3]!="empty")]
samplewell=c()
samplematrix=c(0,0,0,0,0,0)
for(i in sample){
	samplewell=c(samplewell,which(clustering[,1]==i))
	samplestats=c()
	for(j in 1:6){
		samplestats=c(samplestats,length(which(clustering[which(clustering[,1]==i),2]==j)))
	}
	samplestats=samplestats/sum(samplestats)*100
	#postats=colMeans(sample_data[poswell,])
	samplematrix=rbind(samplematrix,samplestats)
}
colnames(samplematrix)=c("1","2","3","4","5","6")
samplematrix1=samplematrix[-1,]
rownames(samplematrix1)=b[,3][which(b[,3]!="empty")]

load('~/Desktop/Projects/ips_kinome_screen/10402 sample_data_ctrltrain.RData')
load('~/Desktop/Projects/ips_kinome_screen/10402 clustering_ctrltrain.RData')
b=read.table(file="~/Desktop/Projects/ips_kinome_screen/data/Christian_IPSC/10402/wellnames_rna.txt")
colnames(sample_data)=c("Nuclear Area", "Nuclear Extent", "max. Oct4 Intensity", "max. Hoechst Intensity", "ø Hoechst Intensity", "ø Oct4 Intensity", "Distance from first neighbor", "Number of Neighbors")
sample=b[,1][which(b[,3]!="empty")]
samplewell=c()
samplematrix=c(0,0,0,0,0,0)
for(i in sample){
	samplewell=c(samplewell,which(clustering[,1]==i))
	samplestats=c()
	for(j in 1:6){
		samplestats=c(samplestats,length(which(clustering[which(clustering[,1]==i),2]==j)))
	}
	samplestats=samplestats/sum(samplestats)*100
	#postats=colMeans(sample_data[poswell,])
	samplematrix=rbind(samplematrix,samplestats)
}
colnames(samplematrix)=c("1","2","3","4","5","6")
samplematrix2=samplematrix[-1,]
rownames(samplematrix2)=b[,3][which(b[,3]!="empty")]

load('~/Desktop/Projects/ips_kinome_screen/10403 sample_data_ctrltrain.RData')
load('~/Desktop/Projects/ips_kinome_screen/10403 clustering_ctrltrain.RData')
b=read.table(file="~/Desktop/Projects/ips_kinome_screen/data/Christian_IPSC/10403/wellnames_rna.txt")
colnames(sample_data)=c("Nuclear Area", "Nuclear Extent", "max. Oct4 Intensity", "max. Hoechst Intensity", "ø Hoechst Intensity", "ø Oct4 Intensity", "Distance from first neighbor", "Number of Neighbors")
sample=b[,1][which(b[,3]!="empty")]
samplewell=c()
samplematrix=c(0,0,0,0,0,0)
for(i in sample){
	samplewell=c(samplewell,which(clustering[,1]==i))
	samplestats=c()
	for(j in 1:6){
		samplestats=c(samplestats,length(which(clustering[which(clustering[,1]==i),2]==j)))
	}
	samplestats=samplestats/sum(samplestats)*100
	#postats=colMeans(sample_data[poswell,])
	samplematrix=rbind(samplematrix,samplestats)
}
colnames(samplematrix)=c("1","2","3","4","5","6")
samplematrix3=samplematrix[-1,]
rownames(samplematrix3)=b[,3][which(b[,3]!="empty")]


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


svg(file="Heatmap_relative_cell_number_104.svg",height=30)
pheatmap(rbind(posmatrix,negmatrix,samplematrix[c(names(sort(dist_to_pos/dist_to_neg))[1:100]),]),cellheight=10)
dev.off()

norm_0_1=function(x){
	return((x-min(x))/(max(x)-min(x)))
	}
tiff(file="Distplot_relative_cell_number_104.tif")
plot(log(dist_to_pos), log(dist_to_neg),pch="*",xlim=c(2,4.5),ylim=c(0,4), xlab="lg ( euclidean distance to positive controls)", ylab="lg ( euclidean distance to negative controls)")
lines(c(log(median(dist_to_pos)), log(median(dist_to_pos))),c(0,max(dist_to_neg)),col="red")
lines(c(0, max(dist_to_pos)),c(log(median(dist_to_neg)), log(median(dist_to_neg))),col="red")
points(log(dist_to_pos)[c(names(sort(dist_to_pos/dist_to_neg))[1:100])],log(dist_to_neg)[c(names(sort(dist_to_pos/dist_to_neg))[1:100])],col="green")
points(log(dist_to_pos)[c("TAF1","SMG1")],log(dist_to_neg)[c("TAF1","SMG1")],col="red",pch="x")
text(c(2.25,2.25,4.25,4.25),c(0,3.75,0,4),c("Oct4-/+ like","Oct4- like","Oct4+ like","Indep. phenotype"))
dev.off()


realtive_hitlist_104=c(names(sort(dist_to_pos/dist_to_neg))[1:100])




#identify(log(dist_to_pos), log(dist_to_neg),labels=names(dist_to_pos),n=3)
