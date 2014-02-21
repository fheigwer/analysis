load('~/Desktop/Projects/ips_kinome_screen/10301 sample_data_ctrl.RData')
load('~/Desktop/Projects/ips_kinome_screen/10301 clustering_ctrl.RData')

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
	#postats=colMeans(sample_data[poswell,])
	posmatrix=rbind(posmatrix,postats)
}
#colnames(posmatrix)=c("1","2","3","4","5","6")
posmatrix=posmatrix[-1,]
rownames(posmatrix)=rep("pos",times=dim(posmatrix)[1])

negwell=c()
negmatrix=c(0,0,0,0,0,0)
for(i in neg){
	negwell=c(negwell,which(clustering[,1]==i))
	negtats=c()
	for(j in 1:6){
		negtats=c(negtats,length(which(clustering[which(clustering[,1]==i),2]==j)))
	}
	#negtats=colMeans(sample_data[negwell,])
	negmatrix=rbind(negmatrix,negtats)
}
colnames(negmatrix)=c("1","2","3","4","5","6")
negmatrix=negmatrix[-1,]
rownames(negmatrix)=rep("neg",times=dim(negmatrix)[1])

#tiff(file="population_heatmap.tif",height=2024,width=2024)
#pheatmap(rbind(posmatrix,negmatrix),mar=c(10,10,10,10),cex=2)
#dev.off()

dev.new()
scatterplot3d(sample_data[poswell,1],sample_data[poswell,2],sample_data[poswell,6],mar=c(10,10,10,10),color=alpha("black",0.3),box=FALSE,pch=".",xlab=colnames(sample_data)[1],ylab=colnames(sample_data)[2],zlab=colnames(sample_data)[6],main="positive controls")
dev.new()
scatterplot3d(sample_data[negwell,1],sample_data[negwell,2],sample_data[negwell,6],mar=c(10,10,10,10),color=alpha("black",0.3),box=FALSE,pch=".",xlab=colnames(sample_data)[1],ylab=colnames(sample_data)[2],zlab=colnames(sample_data)[6],main="negative controls")


library(scales)
library(scatterplot3d)
library(rgl)

x=cbind(rbind(matrix(rnorm(100,mean=1,sd=0.8),nrow=100),matrix(rnorm(100,mean=3,sd=0.8),nrow=100),matrix(rnorm(100,mean=10,sd=0.8),nrow=100))
		,rbind(matrix(rnorm(100,mean=1,sd=0.8),nrow=100),matrix(rnorm(100,mean=3,sd=0.8),nrow=100),matrix(rnorm(100,mean=10,sd=0.8),nrow=100))
		,rbind(matrix(rnorm(100,mean=1,sd=0.8),nrow=100),matrix(rnorm(100,mean=3,sd=0.8),nrow=100),matrix(rnorm(100,mean=10,sd=0.8),nrow=100))
		)
		
y=cbind(	rbind(matrix(rnorm(700,mean=1,sd=0.8),nrow=7),matrix(rnorm(500,mean=3,sd=0.8),nrow=5),matrix(rnorm(300,mean=10,sd=0.8),nrow=3))
		,	rbind(matrix(rnorm(700,mean=1,sd=0.8),nrow=7),matrix(rnorm(500,mean=3,sd=0.8),nrow=5),matrix(rnorm(300,mean=10,sd=0.8),nrow=3))
		,	rbind(matrix(rnorm(700,mean=1,sd=0.8),nrow=7),matrix(rnorm(500,mean=3,sd=0.8),nrow=5),matrix(rnorm(300,mean=10,sd=0.8),nrow=3))
		)
		
z=cbind(	rbind(matrix(rnorm(300,mean=1,sd=0.8),nrow=3),matrix(rnorm(500,mean=3,sd=0.8),nrow=5),matrix(rnorm(700,mean=10,sd=0.8),nrow=7))
		,	rbind(matrix(rnorm(300,mean=1,sd=0.8),nrow=3),matrix(rnorm(500,mean=3,sd=0.8),nrow=5),matrix(rnorm(700,mean=10,sd=0.8),nrow=7))
		,	rbind(matrix(rnorm(300,mean=1,sd=0.8),nrow=3),matrix(rnorm(500,mean=3,sd=0.8),nrow=5),matrix(rnorm(700,mean=10,sd=0.8),nrow=7))
		)
tiff(file="scatter_black.tif")
	scatterplot3d(x,xlab="Oct 4 intensity",ylab="Cell size",zlab="Cell roundness")
dev.off()
 
train_clust=kmeans(x,3,nstart=10)
colnames(train_clust$centers)=c("Oct 4 intensity","Cell size","Cell roundness")
tiff(file="scatter_colored.tif")
 plot3d(x,xlab="Oct 4 intensity",ylab="Cell size",zlab="Cell roundness")#,col=train_clust$cluster)
points3d(train_clust$centers,xlab="Oct 4 intensity",ylab="Cell size",zlab="Cell roundness",col=1:3,size=10)
# points3d( y[1:15,1:3],xlab="Oct 4 intensity",ylab="Cell size",zlab="Cell roundness",col="blue",size=5)
# points3d(z[1:15,1:3],xlab="Oct 4 intensity",ylab="Cell size",zlab="Cell roundness",col="cyan",size=5)
dev.off()

png(file="Starplot_test.png",width = 4024, height = 2024)
	stars(train_clust$centers[,c(1,2,3)],cex=3,bg="grey",scale=TRUE,full=TRUE,draw.segments=TRUE,key.loc=c(-2,4),mar=c(50,50,50,50),col.segments=c("#d53e4f","#e6f598","#fdae61"))
dev.off()

levels_y=knn(x, y[1:15,1:3], train_clust$cluster, k = 100, prob=F)
levels_z=knn(x, z[1:15,1:3], train_clust$cluster, k = 100, prob=F)





poswell=c()
posmatrix=c(3,5,7)
negmatrix=c(7,5,3)
names(posmatrix)=c("black","red","green")
totalmatrix=rbind(posmatrix,negmatrix)
rownames(totalmatrix)=c("sample 1","sample 2")
pheatmap(totalmatrix,co)


for(i in 1:15){
	poswell=c(poswell,which(clustering[,1]==i))
	postats=c()
	for(j in 1:3){
		postats=c(postats,length(which(levels_y==i)))
	}
	#postats=colMeans(sample_data[poswell,])
	posmatrix=rbind(posmatrix,postats)
}

posmatrix=posmatrix[-1,]
rownames(posmatrix)=rep("pos",times=dim(posmatrix)[1])

negwell=c()
negmatrix=c(0,0,0)
for(i in 1:15){
	negwell=c(negwell,which(clustering[,1]==i))
	negtats=c()
	for(j in 1:3){
		negtats=c(negtats,length(which(levels_z==i)))
	}
	#negtats=colMeans(sample_data[negwell,])
	negmatrix=rbind(negmatrix,negtats)
}
colnames(negmatrix)=c("black","red","green")
negmatrix=negmatrix[-1,]
rownames(negmatrix)=rep("neg",times=dim(negmatrix)[1])

tiff(file="population_heatmap.tif",height=2024,width=2024)
pheatmap(rbind(posmatrix,negmatrix),mar=c(10,10,10,10),cex=2)
dev.off()









