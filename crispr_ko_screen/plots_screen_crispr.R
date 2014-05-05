dat=read.table(file="Desktop/Projects/crispr_ko_screen/results/processed_data/CRISPRLib3.0_after sequencing/after_vs_after.tab",sep="\t",header=T)

hist(dat[,3],breaks=1000,ylim=c(0,300),xlim=c(0,5000),col=rgb(1,1,0,0.2))
hist(dat[,6],breaks=1000,ylim=c(0,300),xlim=c(0,5000),col=rgb(0,1,1,0.2),add=T)
random=dat[which(dat[,2]=="random"),]
summary(random)
hist(random[,3],breaks=1000,ylim=c(0,50),xlim=c(0,100),col=rgb(1,1,0,0.2))#
hist(random[,6],breaks=1000,ylim=c(0,50),xlim=c(0,100),col=rgb(0,1,1,0.2),add=T)#,ylim=c(0,50)

gene_dat=c();
gene_names=c();
for(i in levels(dat[,2])){
  gene_dat=rbind(gene_dat,colMeans(dat[which(dat[,2]==i),3:8]))
  gene_names=c(gene_names,i)
}

hist(gene_dat[,3],breaks=1000,ylim=c(0,300),xlim=c(0,500),col=rgb(1,1,0,0.2))
hist(gene_dat[,6],breaks=1000,ylim=c(0,300),xlim=c(0,500),col=rgb(0,1,1,0.2),add=T)

row.names(gene_dat)=gene_names

foldchange=(gene_dat[,4]+gene_dat[,5])/(gene_dat[,1]+gene_dat[,2])
x=seq(1:5000)

plot(log((gene_dat[,1]+gene_dat[,2])),log((gene_dat[,4]+gene_dat[,5])),xlim=c(1,10),ylim=c(1,10),xlab="before",ylab="after")
lines(log(x),log(x))
lines(log(x),log(x/2),col="green")
lines(log(x),log(x*2),col="green")
lines(log(x),log(x/4),col="red")
lines(log(x),log(x*4),col="red")
identify(log((gene_dat[,1]+gene_dat[,2])),log((gene_dat[,4]+gene_dat[,5])),n=10,labels=gene_names)


radius=(sqrt(foldchange)/2/pi)
grid=c()

for(i in seq(1,40,by=1)){
  for(j in seq(1,20,by=1)){
    grid=rbind(grid,cbind(i,j))
  }
}
count=0
x_val=max(grid[,1])+5
leg=c()
text_val=c()
textgrid=c()
for(j in seq(min(foldchange),max(foldchange),length.out=20)){
  grid=rbind(grid,cbind(x_val,count))
  leg=c(leg,sqrt(j)/2/pi)
  text_val=c(text_val,round(j,digits=1))
  textgrid=rbind(textgrid,cbind(x_val+5,count))
  count=count+1
}


#svg(filename="library_coverage.svg")
#tiff(filename="library_coverage_sorteed.tif",width=1012,height=1012)
symbols(grid[,1],grid[,2],circles=c(radius[1:800],leg),inches=0.06,fg="white",bg="red",xlab="",ylab="",cex=3,axes=F)
text(textgrid[,1],textgrid[,2],c(as.character(text_val)),cex=1.2)

