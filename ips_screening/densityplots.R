par(mfrow=c(2,3))


x=density(scaled_controls[which(rownames(scaled_controls)=="pos"),1],kernel="gaussian",bw="SJ")$x
y=density(scaled_controls[which(rownames(scaled_controls)=="pos"),1],kernel="gaussian",bw="SJ")$y
y=(y-min(y))/(max(y)-min(y))
max_pos=x[which(y>=summary(y)[4])[1]]
sd_pos=sd(y)
x_neg=density(scaled_controls[which(rownames(scaled_controls)=="neg"),1],kernel="gaussian",bw="SJ")$x
y_neg=density(scaled_controls[which(rownames(scaled_controls)=="neg"),1],kernel="gaussian",bw="SJ")$y
y_neg=(y_neg-min(y_neg))/(max(y_neg)-min(y_neg))
max_neg=x_neg[which(y_neg>=summary(y_neg)[4])[1]]
sd_neg=sd(y_neg)

plot(x,y,col="green",xlim=c(-3,3),type="l",lty=1,main="Z-Score of population 1",sub=paste("F=",round(abs(max_pos-max_neg)/(abs(3*sd_pos)+abs(abs(3*sd_neg))),digits=3),sep=""))
lines(x_neg,y_neg,col="red",lty=1)



x=density(scaled_controls[which(rownames(scaled_controls)=="pos"),2],kernel="gaussian",bw="SJ")$x
y=density(scaled_controls[which(rownames(scaled_controls)=="pos"),2],kernel="gaussian",bw="SJ")$y
y=(y-min(y))/(max(y)-min(y))
max_pos=x[which(y>=summary(y)[4])[1]]
sd_pos=sd(y)
x_neg=density(scaled_controls[which(rownames(scaled_controls)=="neg"),2],kernel="gaussian",bw="SJ")$x
y_neg=density(scaled_controls[which(rownames(scaled_controls)=="neg"),2],kernel="gaussian",bw="SJ")$y
y_neg=(y_neg-min(y_neg))/(max(y_neg)-min(y_neg))
max_neg=x_neg[which(y_neg>=summary(y_neg)[4])[1]]
sd_neg=sd(y_neg)

plot(x,y,col="green",xlim=c(-3,3),type="l",lty=1,main="Z-Score of population 2",sub=paste("F=",round(abs(max_pos-max_neg)/(abs(3*sd_pos)+abs(abs(3*sd_neg))),digits=3),sep=""))
lines(x_neg,y_neg,col="red",lty=1)



x=density(scaled_controls[which(rownames(scaled_controls)=="pos"),3],kernel="gaussian",bw="SJ")$x
y=density(scaled_controls[which(rownames(scaled_controls)=="pos"),3],kernel="gaussian",bw="SJ")$y
y=(y-min(y))/(max(y)-min(y))
max_pos=x[which(y>=summary(y)[4])[1]]
sd_pos=sd(y)
x_neg=density(scaled_controls[which(rownames(scaled_controls)=="neg"),3],kernel="gaussian",bw="SJ")$x
y_neg=density(scaled_controls[which(rownames(scaled_controls)=="neg"),3],kernel="gaussian",bw="SJ")$y
y_neg=(y_neg-min(y_neg))/(max(y_neg)-min(y_neg))
max_neg=x_neg[which(y_neg>=summary(y_neg)[4])[1]]
sd_neg=sd(y_neg)

plot(x,y,col="green",xlim=c(-3,3),type="l",lty=1,main="Z-Score of population 3",sub=paste("F=",round(abs(max_pos-max_neg)/(abs(3*sd_pos)+abs(abs(3*sd_neg))),digits=3),sep=""))
lines(x_neg,y_neg,col="red",lty=1)



x=density(scaled_controls[which(rownames(scaled_controls)=="pos"),4],kernel="gaussian",bw="SJ")$x
y=density(scaled_controls[which(rownames(scaled_controls)=="pos"),4],kernel="gaussian",bw="SJ")$y
y=(y-min(y))/(max(y)-min(y))
max_pos=x[which(y>=summary(y)[4])[1]]
sd_pos=sd(y)
x_neg=density(scaled_controls[which(rownames(scaled_controls)=="neg"),4],kernel="gaussian",bw="SJ")$x
y_neg=density(scaled_controls[which(rownames(scaled_controls)=="neg"),4],kernel="gaussian",bw="SJ")$y
y_neg=(y_neg-min(y_neg))/(max(y_neg)-min(y_neg))
max_neg=x_neg[which(y_neg>=summary(y_neg)[4])[1]]
sd_neg=sd(y_neg)

plot(x,y,col="green",xlim=c(-3,3),type="l",lty=1,main="Z-Score of population 4",sub=paste("F=",round(abs(max_pos-max_neg)/(abs(3*sd_pos)+abs(abs(3*sd_neg))),digits=3),sep=""))
lines(x_neg,y_neg,col="red",lty=1)



x=density(scaled_controls[which(rownames(scaled_controls)=="pos"),5],kernel="gaussian",bw="SJ")$x
y=density(scaled_controls[which(rownames(scaled_controls)=="pos"),5],kernel="gaussian",bw="SJ")$y
y=(y-min(y))/(max(y)-min(y))
max_pos=x[which(y>=summary(y)[4])[1]]
sd_pos=sd(y)
x_neg=density(scaled_controls[which(rownames(scaled_controls)=="neg"),5],kernel="gaussian",bw="SJ")$x
y_neg=density(scaled_controls[which(rownames(scaled_controls)=="neg"),5],kernel="gaussian",bw="SJ")$y
y_neg=(y_neg-min(y_neg))/(max(y_neg)-min(y_neg))
max_neg=x_neg[which(y_neg>=summary(y_neg)[4])[1]]
sd_neg=sd(y_neg)

plot(x,y,col="green",xlim=c(-3,3),type="l",lty=1,main="Z-Score of population 5",sub=paste("F=",round(abs(max_pos-max_neg)/(abs(3*sd_pos)+abs(abs(3*sd_neg))),digits=3),sep=""))
lines(x_neg,y_neg,col="red",lty=1)

x=density(scaled_controls[which(rownames(scaled_controls)=="pos"),6],kernel="gaussian",bw="SJ")$x
y=density(scaled_controls[which(rownames(scaled_controls)=="pos"),6],kernel="gaussian",bw="SJ")$y
y=(y-min(y))/(max(y)-min(y))
max_pos=x[which(y>=summary(y)[4])[1]]
sd_pos=sd(y)
x_neg=density(scaled_controls[which(rownames(scaled_controls)=="neg"),6],kernel="gaussian",bw="SJ")$x
y_neg=density(scaled_controls[which(rownames(scaled_controls)=="neg"),6],kernel="gaussian",bw="SJ")$y
y_neg=(y_neg-min(y_neg))/(max(y_neg)-min(y_neg))
max_neg=x_neg[which(y_neg>=summary(y_neg)[4])[1]]
sd_neg=sd(y_neg)

plot(x,y,col="green",xlim=c(-3,3),type="l",lty=1,main="Z-Score of population 6",sub=paste("F=",round(abs(max_pos-max_neg)/(abs(3*sd_pos)+abs(abs(3*sd_neg))),digits=3),sep=""))
lines(x_neg,y_neg,col="red",lty=1)
