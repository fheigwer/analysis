
name_table=read.table(file="/Users/b110-mm06/Desktop/Projects/ips_kinome_screen/data/Christian_IPSC/10303/wellnames.txt")
pos_wells=which(name_table[,3]=="pos")
neg_wells=which(name_table[,3]=="neg")
load("/Users/b110-mm06/Desktop/Projects/ips_kinome_screen/10403 clustering.RData")

result=c(0,0,0,0,0,0)
for(i in as.numeric(levels(factor(clustering[,1])))){
	temp=hist(clustering[which(clustering[,1]==i),2],breaks=c(0,as.numeric(levels(factor(clustering[,2])))),plot=FALSE)$count
	temp=temp #/sum(temp)*100
	result=rbind(result,temp)
}
result_sample=result[-1,]

load("/Users/b110-mm06/Desktop/Projects/ips_kinome_screen/10303 clustering_ctrl.RData")

result=c(0,0,0,0,0,0)
for(i in as.numeric(levels(factor(clustering[,1])))){
	temp=hist(clustering[which(clustering[,1]==i),2],breaks=c(0,as.numeric(levels(factor(clustering[,2])))),plot=FALSE)$count
	temp=t(temp)/sum(temp)*100
	if(i %in% pos_wells){
		rownames(temp)="pos"
	}else{		
		rownames(temp)="neg"
	}
	result=rbind(result,temp)
	
}
result_ctrl=result[-1,]
#rownames(result_ctrl[,])=rep("pos_ctrl",dim(result_ctrl[which(name_table[,3]=="pos"),])[1])
#rownames(result_ctrl[which(name_table[,3]=="neg"),])=rep("neg_ctrl",dim(result_ctrl[which(name_table[,3]=="pos"),])[1])
dev.new()
heatmap.2(rbind(result_ctrl,result_sample),scale="column")


load("/Users/b110-mm06/Desktop/Projects/ips_kinome_screen/10303 clustering.RData")

result=c(0,0,0,0,0,0)
for(i in as.numeric(levels(factor(clustering[,1])))){
	temp=hist(clustering[which(clustering[,1]==i),2],breaks=c(0,as.numeric(levels(factor(clustering[,2])))),plot=FALSE)$count
	temp=temp #/sum(temp)*100
	result=rbind(result,temp)
}
result_sample=result[-1,]

load("/Users/b110-mm06/Desktop/Projects/ips_kinome_screen/10303 clustering_ctrl.RData")

result=c(0,0,0,0,0,0)
for(i in as.numeric(levels(factor(clustering[,1])))){
	temp=hist(clustering[which(clustering[,1]==i),2],breaks=c(0,as.numeric(levels(factor(clustering[,2])))),plot=FALSE)$count
	temp=temp #/sum(temp)*100
	result=rbind(result,temp)
}
result_ctrl=result[-1,]
rownames(result_ctrl)=rep("ctrl",dim(result_ctrl)[1])
dev.new()
heatmap.2(rbind(result_ctrl,result_sample),scale="column")
