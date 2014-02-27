library(misc3d)
library(scales)
library(scatterplot3d)
load('~/Desktop/Projects/ips_kinome_screen/raw_controls.RData')
load('~/Desktop/Projects/ips_kinome_screen/kmeans_result.RData')

colnames(raw_clust$centers)=c("Nuclear Area", "Nuclear Extent", "max. Oct4 Intensity", "max. Hoechst Intensity", "ø Hoechst Intensity", "ø Oct4 Intensity", "Distance from first neighbor", "Number of Neighbors")
colnames(raw_data)=c("Nuclear Area", "Nuclear Extent", "max. Oct4 Intensity", "max. Hoechst Intensity", "ø Hoechst Intensity", "ø Oct4 Intensity", "Distance from first neighbor", "Number of Neighbors")


tiff(file="Area_MaxInt_Dist_train.tif",width = 2024, height = 2024)
	scatterplot3d(raw_data[,1],raw_data[,3],raw_data[,7],mar=c(10,10,10,10),color=alpha(raw_clust$cluster,0.5),box=FALSE,pch=".",xlab=colnames(raw_data)[1],ylab=colnames(raw_data)[3],zlab=colnames(raw_data)[7],cex.symbols=3, cex.axis=3,cex.lab=2.5)
dev.off()


tiff(file="Areaext_Meanint_Neighbour_train.tif",width = 2024, height = 2024)
	scatterplot3d(raw_data[,2],raw_data[,8],raw_data[,6],mar=c(10,10,10,10),color=alpha(raw_clust$cluster,0.5),box=FALSE,pch=".",xlab=colnames(raw_data)[2],ylab=colnames(raw_data)[8],zlab=colnames(raw_data)[6],cex.symbols=3, cex.axis=3,cex.lab=2.5)
dev.off()

tiff(file="MeanintG_MeanintB_Max_train.tif",width = 2024, height = 2024)
	scatterplot3d(raw_data[,5],raw_data[,6],raw_data[,3],mar=c(10,10,10,10),color=alpha(raw_clust$cluster,0.5),box=FALSE,pch=".",xlab=colnames(raw_data)[5],ylab=colnames(raw_data)[6],zlab=colnames(raw_data)[3],cex.symbols=3, cex.axis=3,cex.lab=2.5)
dev.off()

png(file="Starplot.png",width = 4024, height = 2024)
	stars(raw_clust$centers[,c(1,2,3,6,4,5,7,8)],cex=3,bg="grey",scale=TRUE,full=TRUE,draw.segments=TRUE,key.loc=c(-1.6,5),mar=c(1,30,1,0),col.segments=c("#d53e4f","#f46d43","#abdda4","#e6f598","#3288bd","#66c2a5","#fee08b","#fdae61"))
dev.off()

################################################################################################################################################################
################################################################################################################################################################
for(i in 3:4){
	for(j in 1:3){
#		load(paste('~/Desktop/Projects/ips_kinome_screen/10',i,'0',j,' sample_data_ctrl_ctrltrain.RData',sep=""))
#		load(paste('~/Desktop/Projects/ips_kinome_screen/10',i,'0',j,' clustering_ctrl_ctrltrain.RData',sep=""))
#		colnames(sample_data)=c("Nuclear Area", "Nuclear Extent", "max. Oct4 Intensity", "max. Hoechst Intensity", "ø Hoechst Intensity", "ø Oct4 Intensity", "Distance from first neighbor", "Number of Neighbors")
#		sample=sample(1:(length(clustering[,2])),size=50000)
#		tiff(file=paste("paired_plot_",i,"_",j,"_ctrl.tif",sep=""),width = 2024, height = 2024)
#			pairs(sample_data[sample,c(1,2,3,6,4,5,7,8)],col=alpha(clustering[sample,2],0.3),pch=".",lower.panel=NULL,cex=3)
#		dev.off()
#		test=list()

#for(o in 1:8){
#	for(p in c(1:8)[-o]){
#		for(q in c(1:8)[-c(o,p)]){
#			if(!is.numeric(test[[paste(as.character(sort(c(o,p,q))),collapse = '')]])){
#				tiff(file=paste(colnames(sample_data)[o],"_",colnames(sample_data)[p],"_",colnames(sample_data)[q],"_",i,"_",j,"_ctrl.tif",sep=""),width = 2024, height = 2024)
#						scatterplot3d(sample_data[,o],sample_data[,p],sample_data[,q],mar=c(10,10,10,10),color=alpha(clustering[,2],0.3),box=FALSE,pch=".",xlab=colnames(sample_data)[o],ylab=colnames(sample_data)[p],zlab=colnames(sample_data)[q],cex.symbols=3, cex.axis=3,cex.lab=2.5)
#				dev.off()
#				test[[as.character(paste(as.character(sort(c(o,p,q))),collapse = ''))]]=test[[as.character(paste(as.character(sort(c(o,p,q))),collapse = ''))]]+1
#			}
#		}		
#	}
#}
		################################################################################################################################################################
		load(paste('~/Desktop/Projects/ips_kinome_screen/10',i,'0',j,' sample_data_ctrltrain.RData',sep=""))
		load(paste('~/Desktop/Projects/ips_kinome_screen/10',i,'0',j,' clustering_ctrltrain.RData',sep=""))
		colnames(sample_data)=c("Nuclear Area", "Nuclear Extent", "max. Oct4 Intensity", "max. Hoechst Intensity", "ø Hoechst Intensity", "ø Oct4 Intensity", "Distance from first neighbor", "Number of Neighbors")
		sample=sample(1:(length(clustering[,2])),size=50000)
		tiff(file=paste("paired_plot_",i,"_",j,".tif",sep=""),width = 2024, height = 2024)
			pairs(sample_data[sample,c(1,2,3,6,4,5,7,8)],col=alpha(clustering[sample,2],0.3),pch=".",lower.panel=NULL,cex=3)
		dev.off()
		test=list()
		for(o in 1:8){
			for(p in c(1:8)[-o]){
				for(q in c(1:8)[-c(o,p)]){
					if(!is.numeric(test[[paste(as.character(sort(c(o,p,q))),collapse = '')]])){
						tiff(file=paste(colnames(sample_data)[o],"_",colnames(sample_data)[p],"_",colnames(sample_data)[q],"_",i,"_",j,".tif",sep=""),width = 2024, height = 2024)
							scatterplot3d(sample_data[,o],sample_data[,p],sample_data[,q],mar=c(10,10,10,10),color=alpha(clustering[,2],0.3),box=FALSE,pch=".",xlab=colnames(sample_data)[o],ylab=colnames(sample_data)[p],zlab=colnames(sample_data)[q],cex.symbols=3, cex.axis=3,cex.lab=2.5)
						dev.off()
						test[[as.character(paste(as.character(sort(c(o,p,q))),collapse = ''))]]=test[[as.character(paste(as.character(sort(c(o,p,q))),collapse = ''))]]+1
					}
				}		
			}
		}
	}
}	
	
		

