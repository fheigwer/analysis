par(mfrow=c(2,7))
for(i in 1:14){
	plot(raw_clust$centers[,i],type="b",col=c(1:6),ylab=colnames(raw_clust$centers)[i],xlab="Population")
}