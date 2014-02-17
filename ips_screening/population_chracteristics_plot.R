par(mfrow=c(2,4))
for(i in 1:8){
	plot(raw_clust$centers[,i],type="b",col=c(1:6),ylab=colnames(raw_clust$centers)[i],xlab="Population")
}