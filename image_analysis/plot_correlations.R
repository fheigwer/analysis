library(gplots)
time1=read.table(file="1_20131223_HumanDmel_test_X039_BD19.tab",header=TRUE,sep="\t",dec=".",row.names = 1)#[-c(257),]
time2=read.table(file="2_20131223_HumanDmel_test_X039_BD19.tab",header=TRUE,sep="\t",dec=".",row.names = 1)
time3=read.table(file="3_20131223_HumanDmel_test_X039_BD19.tab",header=TRUE,sep="\t",dec=".",row.names = 1)
time4=read.table(file="4_20131223_HumanDmel_test_X039_BD19.tab",header=TRUE,sep="\t",dec=".",row.names = 1)

pdf(file="1_vs_2_X039_BD19.pdf")
par(mfrow=c(4,5),oma=c(0,0,6,0))
outliers=list("schnatz"=1)
for(i in 1:20){
	smoothScatter(time1[,i],time2[,i],main=colnames(time1)[i],xlab="run 1",ylab="run 2",sub=paste("spear.corr.=",round(cor(time1[,i],time2[,i],method="spearman"),digits=2),sep=""))
	for(j in rownames(time1)[which(lm(time1[,i]~time2[,i])$residuals<quantile(lm(time1[,i]~time2[,i])$residuals,probs=c(0.05,0.95))[1])]){
		if(is.null(outliers[[j]]) ){
			outliers[[j]]=0
		}else{
			outliers[[j]]=outliers[[j]]+1
		}
	}
	for(j in rownames(time1)[which(lm(time1[,i]~time2[,i])$residuals>quantile(lm(time1[,i]~time2[,i])$residuals,probs=c(0.05,0.95))[2])]){
		if(is.null(outliers[[j]]) ){
			outliers[[j]]=0
		}else{
			outliers[[j]]=outliers[[j]]+1
		}
	}
}
test=""
bla=0
for(i in names(outliers)){
	if(outliers[[i]]>3){
		if(bla<12){
			test=paste(test,i,sep="_")
			bla=bla+1
		}else{
			test=paste(test,i,sep="\n")
			bla=0
		}
		
	}
}
title(test, outer=TRUE)
dev.off()



