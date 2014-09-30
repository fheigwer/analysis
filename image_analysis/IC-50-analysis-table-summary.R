a=read.table("~/Desktop/Projects/time_resolved_sgi/PD-0325901-EC50-Image-based/IC-50-results/plate_1/all_1.tab",row.names=c(1),header=TRUE)
b=c()
for(i in seq(1,nrow(a),by=4)){
  b=rbind(b,colMeans(a[i:(i+3),]))
}

x=40
conc=c()
for(i in 2:24){
  x=x/2
  conc=c(conc,x)
}
conc=c(conc,0)
se <- function(x) sqrt(var(x)/length(x))

Resultmatrix_1=c()
for(feature in 1:ncol(b)){
  feature_vector=c()
  for(i in 1:24){
    j=i
    temp=c()
    while(j<=nrow(b)){
      temp=c(temp,t(b)[feature,j])
      j=j+24
    }
    feature_vector=c(feature_vector,mean(temp))
  }
  Resultmatrix_1=rbind(Resultmatrix_1,feature_vector)
}
row.names(Resultmatrix_1)=colnames(a)
row.names(sd_matrix)=colnames(a)
Resultmatrix_1=cbind(Resultmatrix_1,Resultmatrix_1[,1])
Resultmatrix_1=Resultmatrix_1[,-1]

feature_vector_1=list()
for(i in 1:24){
    j=i
    temp=c()
    while(j<=nrow(b)){
      temp=c(temp,t(b)[1,j])
      j=j+24
    }
    feature_vector_1[[i]]=temp
}

a=read.table("~/Desktop/Projects/time_resolved_sgi/PD-0325901-EC50-Image-based/IC-50-results/plate_2/all_2.tab",row.names=c(1),header=TRUE)
b=c()
for(i in seq(1,nrow(a),by=4)){
  b=rbind(b,colMeans(a[i:(i+3),]))
}

Resultmatrix_2=c()
for(feature in 1:ncol(b)){
  feature_vector=c()
  for(i in 1:24){
    j=i
    temp=c()
    while(j<=nrow(b)){
      temp=c(temp,t(b)[feature,j])
      j=j+24
    }
    feature_vector=c(feature_vector,mean(temp))
  }
  Resultmatrix_2=rbind(Resultmatrix_2,feature_vector)
}
row.names(Resultmatrix_2)=colnames(a)
row.names(sd_matrix)=colnames(a)
Resultmatrix_2=cbind(Resultmatrix_2,Resultmatrix_2[,1])
Resultmatrix_2=Resultmatrix_2[,-1]

for(i in 1:24){
  j=i
  temp=c()
  while(j<=nrow(b)){
    temp=c(temp,t(b)[1,j])
    j=j+24
  }
  feature_vector_1[[i]]=c(feature_vector_1[[i]],temp)
}


for(i in 1:nrow(Resultmatrix_1)){
  mi=min(c(Resultmatrix_1[i,],Resultmatrix_2[i,]))
  ma=max(c(Resultmatrix_1[i,],Resultmatrix_2[i,]))
  Resultmatrix_1[i,]=(Resultmatrix_1[i,]-mi)/(ma-mi) 
  Resultmatrix_2[i,]=(Resultmatrix_2[i,]-mi)/(ma-mi) 
}
  
  correlating_features=c()
  cor_matrix=cor(t(Resultmatrix_1[,c(1:19,21,22)]),t(Resultmatrix_2[,c(1:19,21,22)]),method="spearman")
  for(i in 1:nrow(Resultmatrix_1)){
    if(abs(cor_matrix[i,i])>0.8){
      correlating_features=c(correlating_features,row.names(cor_matrix)[i])
    }
  }
  
  par(mfrow=c(1,1))
  for(i in correlating_features){   
      plot(log(conc[c(1:19,21,22)]),Resultmatrix_1[i,c(1:19,21,22)],type="b",ylab=i,xlab="log[PD-0325901 nM]",ylim=c(0,1),xlim=c(3,-12))
       lines(log(conc[c(1:19,21,22)]),Resultmatrix_2[i,c(1:19,21,22)],type="b",col="red")
  }


x=40
feature_vector_2=c()
for(i in 2:24){
  x=x/2
  feature_vector_2=rbind(feature_vector_2,cbind(feature_vector_1[[i]],rep(x,times=length(feature_vector_1[[i]]))))
}
feature_vector_2=rbind(feature_vector_2,cbind(feature_vector_1[[1]],rep(0,times=length(feature_vector_1[[1]]))))

x=40
feature_vector_3=c()
for(i in 2:23){
  x=x/2
  feature_vector_3=rbind(feature_vector_3,cbind(feature_vector_1[[i]],rep(x,times=length(feature_vector_1[[i]]))))
}

par(mfrow=c(1,4))
test=as.data.frame(cbind(conc,Resultmatrix_1["cells",]))
colnames(test)=c("conc","cells")
mdl=drm(cells ~ conc, data = test, fct = LL.4(names=c("Slope","Lower Limit","Upper Limit", "ED50")))
ED(mdl,50,interval="delta")*1000
plot(mdl,sub="first replicate mean -> model")
#par(mfrow=c(1,1))
test1=as.data.frame(cbind(conc,Resultmatrix_2["cells",]))#as.data.frame(cbind(rev(conc[c(1:19,21,22)]),c(Resultmatrix_2["cells",c(1)],rev(Resultmatrix_2["cells",c(2:20,22,23)]))))
colnames(test1)=c("conc","cells")
mdl=drm(cells ~ conc, data = test1, fct = LL.4(names=c("Slope","Lower Limit","Upper Limit", "ED50")))
ED(mdl,50,interval="delta")*1000
plot(mdl,sub="second replicate mean -> model")
test1=as.data.frame(feature_vector_2)
colnames(test1)=c("cells","conc")
row.names(test1)=c()
mdl=drm(cells ~ conc, data = test1, fct = LL.4(names=c("Slope","Lower Limit","Upper Limit", "ED50")))
ED(mdl,50,interval="delta")*1000
plot(mdl,sub="all wells independend -> model")
test1=as.data.frame(feature_vector_3)
colnames(test1)=c("cells","conc")
row.names(test1)=c()
mdl=drm(cells ~ conc, data = test1, fct = LL.4(names=c("Slope","Lower Limit","Upper Limit", "ED50")))
ED(mdl,50,interval="delta")*1000
plot(mdl,sub="2-23 wells independend -> model")






