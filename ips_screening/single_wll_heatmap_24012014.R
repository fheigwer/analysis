library(stats)
library(gplots)
library(ggplot2)
library(sm)
library(class)

# load data  G3
POSITIVE <- read.table("/Users/b110-mm06/Desktop/Projects/ips_kinome_screen/data/singlecell/sample_well_data_131204/10301_G03/DefaultOUT_Nuclei.csv",sep="\t",dec=".",header=TRUE)
POSITIVE<-POSITIVE[(which(!is.na(rowMeans(POSITIVE)))),]
NEGATIVE <- read.table("/Users/b110-mm06/Desktop/Projects/ips_kinome_screen/data/singlecell/sample_well_data_131204/10301_G04/DefaultOUT_Nuclei.csv",sep="\t",dec=".",header=TRUE)
NEGATIVE<-NEGATIVE[(which(!is.na(rowMeans(NEGATIVE)))),]
SAMPLE <- read.table("/Users/b110-mm06/Desktop/Projects/ips_kinome_screen/data/singlecell/sample_well_data_131204/10303_A07/DefaultOUT_Nuclei.csv",sep="\t",dec=".",header=TRUE)
SAMPLE<-SAMPLE[(which(!is.na(rowMeans(SAMPLE)))),]

TRAINDATA=rbind(POSITIVE,NEGATIVE)

drop_v <- as.vector(c(
			grep(".*AreaShape_EulerNumber.*", colnames(TRAINDATA), perl=T),
            grep(".*Children_Green_Count.*", colnames(TRAINDATA), perl=T),
            grep(".*Children_Red_Count.*", colnames(TRAINDATA), perl=T),
            grep(".*Number_Object_Number.*", colnames(TRAINDATA), perl=T),
            grep(".*id.*", colnames(TRAINDATA), perl=T),
            grep(".*ImageNumber.*", colnames(TRAINDATA), perl=T),
            grep(".*ObjectNumber.*", colnames(TRAINDATA), perl=T),
            grep(".*Location_Center_X.*", colnames(TRAINDATA), perl=T),
            grep(".*Location_Center_Y.*", colnames(TRAINDATA), perl=T),
            grep(".*Parent_Nuclei.*", colnames(TRAINDATA), perl=T),
            grep("AreaShape_Center_X", colnames(TRAINDATA), perl=T),
            grep("AreaShape_Center_Y", colnames(TRAINDATA), perl=T),
            grep(".*AreaShape_Orientation.*", colnames(TRAINDATA),perl=T),
            grep(".*AreaShape_Compactness.*", colnames(TRAINDATA),perl=T)
            ))
TRAINDATA <- TRAINDATA [,-c(drop_v)]

shape_v=as.vector(c(grep(".*AreaShape.*", colnames(TRAINDATA),perl=T)))
fluo_v=as.vector(c(grep(".*Intensity.*", colnames(TRAINDATA),perl=T)))
neighbors_v=as.vector(c(grep(".*Neighbors.*", colnames(TRAINDATA),perl=T)))
textures_v=as.vector(c(grep(".*Texture.*", colnames(TRAINDATA),perl=T)))

shape_data <- as.matrix(TRAINDATA)[,sort(colnames(as.matrix(TRAINDATA)[,shape_v]))]
fluo_data <- as.matrix(TRAINDATA)[,sort(colnames(as.matrix(TRAINDATA)[,fluo_v]))]
neighbors_data <- as.matrix(TRAINDATA)[,sort(colnames(as.matrix(TRAINDATA)[,neighbors_v]))]
textures_data <- as.matrix(TRAINDATA)[,sort(colnames(as.matrix(TRAINDATA)[,textures_v]))]

TRAINDATA=cbind(shape_data,fluo_data,neighbors_data)#,textures_data)

correlation_matrix=cor(TRAINDATA,method="spearman")

to_drop=list()
for(i in 1:dim(correlation_matrix)[1]){
	temp=correlation_matrix[correlation_matrix[,i]>0.90 ,i]
	temp=temp[temp<1]
	for(tempele in names(temp)){
		to_drop[[tempele]]<-1
	}
}

for(i in as.vector(names(to_drop))){
	if(!is.na(i)){
		TRAINDATA <- TRAINDATA[,c(which(colnames(TRAINDATA)!=i))]
	}
}


SAMPLE <- SAMPLE [,-c(drop_v)]
shape_data <- as.matrix(SAMPLE)[,sort(colnames(as.matrix(SAMPLE)[,shape_v]))]
fluo_data <- as.matrix(SAMPLE)[,sort(colnames(as.matrix(SAMPLE)[,fluo_v]))]
neighbors_data <- as.matrix(SAMPLE)[,sort(colnames(as.matrix(SAMPLE)[,neighbors_v]))]
textures_data <- as.matrix(SAMPLE)[,sort(colnames(as.matrix(SAMPLE)[,textures_v]))]
SAMPLE=cbind(shape_data,fluo_data,neighbors_data)#,textures_data)

for(i in as.vector(names(to_drop))){
	if(!is.na(i)){
		SAMPLE <- SAMPLE[,c(which(colnames(SAMPLE)!=i))]
	}
}





TRAINED_CLUST <-kmeans(TRAINDATA,6, nstart = 30)

SAMPLE.levels=knn(TRAINDATA, SAMPLE, TRAINED_CLUST$cluster, k = 100, prob=F) 


plot(density(as.numeric(TRAINED_CLUST$cluster)),ylim=c(0,0.5),col="green")
lines(density(as.numeric(SAMPLE.levels)),col="red")




#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#pairs(combined_data[random_row,],col=combined_clust$cluster[random_row])
#g3_4 <- combined_clust$centers
#for( i in 1:dim(combined_clust$centers)[1]){
#	for(j in 1:dim(combined_clust$centers)[2]){
#		g3_4[i,j]=(g3_4[i,j]-min(combined_clust$centers[,j]))/(max(combined_clust$centers[,j])-min(combined_clust$centers[,j]))
#	}
#}
#dev.new()
#####################################################################################################################################
pc=princomp(TRAINDATA)
par(mfrow=c(2,7))
for(i in 1:14){
	sm.density.compare(TRAINDATA[,i], TRAINED_CLUST$cluster)
}

plot3d(pc$scores[,c(1,2,3)], col=TRAINED_CLUST$cluster)
#####################################################################################################################################
plot(g3_4[1,1:30],type="l",ylim=c(0,1.5))
lines(g3_4[2,1:30],type="l",col=rainbow(7)[1])
lines(g3_4[3,1:30],type="l",col=rainbow(7)[2])
lines(g3_4[4,1:30],type="l",col=rainbow(7)[3])
lines(g3_4[5,1:30],type="l",col=rainbow(7)[4])
lines(g3_4[6,1:30],type="l",col=rainbow(7)[5])
lines(g3_4[7,1:30],type="l",col=rainbow(7)[6])
                 
train <- G3.data	
A7.data <- combined_data
cl <- G3.clust$cluster
A7.levels=knn(train, A7.data, cl, k = 30, prob=TRUE)                 
pc=princomp(combined_data)
plot3d(pc$scores[,c(1,2,3)], col=G4.levels)

plot(density(G3.clust$cluster),ylim=c(0,1))
lines(density(as.numeric(G4.levels)),col="green")
lines(density(as.numeric(A7.levels)),col="red")

                 
                 
            
           