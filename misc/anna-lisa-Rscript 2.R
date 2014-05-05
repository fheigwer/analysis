anna.lisa.daten=read.table("~/Desktop/anna-lisa-daten.txt",sep="\t",dec=".",header=TRUE)
for(i in 1:dim(anna.lisa.daten)[1])(
  rownames(anna.lisa.daten)[i]=paste(as.character(anna.lisa.daten[i,1]),anna.lisa.daten[i,2],sep="_")
)
anna.lisa.daten=as.matrix(anna.lisa.daten[,7:(dim(anna.lisa.daten)[2])])
anna.lisa.daten=anna.lisa.daten[!is.na(rowMeans(anna.lisa.daten)),]
for(i in 1:ncol(anna.lisa.daten)){
  anna.lisa.daten[,i]=(anna.lisa.daten[,i]-min(anna.lisa.daten[,i]))/(max(anna.lisa.daten[,i])-min(anna.lisa.daten[,i]))
}
clust=kmeans(anna.lisa.daten,10,nstart=10)
 tiff("pairs.tif",height=2024,width=2024)
 pairs(anna.lisa.daten,col=clust$cluster)
 dev.off()
