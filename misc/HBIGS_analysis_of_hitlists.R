wnt3a=read.table(file="C:\\Users\\Flo\\Desktop\\Results\\Wnt3a\\in\\topTable.txt",header=T)
wnt3=read.table(file="C:\\Users\\Flo\\Desktop\\Results\\Wnt3\\in\\topTable.txt",header=T)
control=read.table(file="C:\\Users\\Flo\\Desktop\\Results\\Control\\Dual_channel\\in\\topTable.txt",header=T)


names=wnt3a[,22]

cols=rep("black",time=length(names))

cols[which(control[,5]=="pos")]="green"

cols[which(control[,5]=="neg")]="red"

wnt3a=wnt3a[,c(4,7,8,9,10,19,20)]
wnt3a=wnt3a[which(!is.na(rowMeans(wnt3a))),]
wnt3=wnt3[,c(4,7,8,9,10,19,20)]
wnt3=wnt3[which(!is.na(rowMeans(wnt3))),]
control=control[1:1056,c(4,7,8,9,10,19,20)]
control=control[which(!is.na(rowMeans(control))),]

x=seq(0,1000,by=0.01)
par(mfrow=c(1,2))
plot(wnt3[,"normalized_r1_ch1"],wnt3[,"normalized_r2_ch1"],col=cols,main="normal scale",xlab="replicate 1",ylab="replicate 2")
lines(x,x,col="red")

plot(log(wnt3[,"normalized_r1_ch1"]),log(wnt3[,"normalized_r2_ch1"]),col=cols,main="LOG scale",xlab="replicate 1",ylab="replicate 2")
lines(log(x),log(x),col="red")


summary(lm(log(wnt3[,"normalized_r1_ch1"])~log(wnt3[,"normalized_r2_ch1"])))

cor(wnt3[,"normalized_r1_ch1"],wnt3[,"normalized_r2_ch1"],method="pearson")

avg_wnt3=(wnt3[,"normalized_r1_ch1"]+wnt3[,"normalized_r2_ch1"])/2
avg_wnt3a=(wnt3a[,"normalized_r1_ch1"]+wnt3a[,"normalized_r2_ch1"])/2
avg_control=(control[,"normalized_r1_ch1"]+control[,"normalized_r2_ch1"])/2

pairs(cbind(avg_wnt3,avg_wnt3a,avg_control),col=cols)
heatmap(cbind(avg_wnt3,avg_wnt3a,avg_control))
wnt3_rat=avg_wnt3/avg_control
wnt3a_rat=avg_wnt3a/avg_control

control_rat=avg_control/avg_control

plot(wnt3_rat[order(wnt3_rat,decreasing =TRUE)],wnt3a_rat[order(wnt3_rat,decreasing =TRUE)],xlab="Wnt3 / Control",ylab="Wnt3a / Control")

wnt3_hits=names[order(wnt3_rat,decreasing =TRUE)[1:100]]
wnt3_hits=wnt3_hits[which(!is.na(wnt3_hits))]

wnt3a_hits=names[order(wnt3a_rat,decreasing =TRUE)[1:100]]
wnt3a_hits=wnt3a_hits[which(!is.na(wnt3a_hits))]

control_hits=names[order(control_rat,decreasing =TRUE)[1:100]]
control_hits=control_hits[which(!is.na(wnt3a_hits))]

wnt3_rat_z=(wnt3_rat-mean(wnt3_rat))/sd(wnt3_rat)
wnt3a_rat_z=(wnt3a_rat-mean(wnt3a_rat))/sd(wnt3a_rat)

length(which(wnt3_rat_z>=1))
length(which(wnt3a_rat_z>=1))

length(intersect(which(wnt3_rat_z>=1),which(wnt3a_rat_z>=1)))

draw.pairwise.venn(	length(which(wnt3_rat_z>=1)), 
				length(which(wnt3a_rat_z>=1)), 
				length(intersect(which(wnt3_rat_z>=1),which(wnt3a_rat_z>=1))),
				c("Wnt3", "Wnt3a"),fill = c("blue", "red"),
lty = "blank",
cex = 2,
cat.cex = 2,
cat.pos = c(285, 105),
cat.dist = 0.09,
cat.just = list(c(-1, -1), c(1, 1)),
ext.pos = 30,
ext.dist = -0.05,
ext.length = 0.85,
ext.line.lwd = 2,
ext.line.lty = "dashed")




