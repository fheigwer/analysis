
orig_data=Dmel_HP_20x_4tiles_ERC_HD3A_1047_S11_2014.09.12.12.35.28_results
  
  #Dmel_HP_20x_4tiles_ERC_HD3A_2005_S11_2014.09.12.16.55.17_results#
#Dmel_HP_20x_4tiles_ERC_HD3A_1045_S11_2014.09.12.10.51.01_results

rnames=aggregate(orig_data[,1],by = list(orig_data[,1]),mean)[,1]


mean_tab=aggregate(orig_data[,3:ncol(orig_data)],by = list(orig_data[,1]),mean)
row.names(mean_tab)=rnames


mean_tab=mean_tab[,-1]
mean_tab=mean_tab[,-grep("*cx$",colnames(mean_tab))]
mean_tab=mean_tab[,-grep("*cy$",colnames(mean_tab))]
mean_tab=mean_tab[,-grep("*s2$",colnames(mean_tab))]
                  

mean_tab_ctrl=apply(mean_tab, 2, function(x) x=(x/mean(x[c("C2","D2","K2","L2")])))

### very riskyy attempt the plates have different siRNAs spotted so an in plate Z-score wdoes not make sense to be compared between plates
mean_tab_z=apply(mean_tab_ctrl, 2, function(x) x=(x-mean(x))/(sd(x))) ###




library(scales)
pdf("Dmel_HP_20x_4tiles_ERC_HD3A_1047_S11_2014.09.12.12.35.28_ctrls.pdf")
par(mar=c(10,2,10,2))
x<-barplot(mean_tab_z["A1",],col=alpha("red",0.5),las=2,ylim=c(-3,3),main="dMEK")
barplot(mean_tab_z["B1",],add=T,col=alpha("green",0.5),yaxt="n",xaxt="n")
lines(y=c(2,2),x=c(min(x),max(x)),lty=1,lwd=0.5)
lines(y=c(1,1),x=c(min(x),max(x)),lty=2,lwd=0.5)
lines(y=c(-1,-1),x=c(min(x),max(x)),lty=2,lwd=0.5)
lines(y=c(-2,-2),x=c(min(x),max(x)),lty=1,lwd=0.5)

barplot(mean_tab_z["A2",],col=alpha("red",0.5),las=2,ylim=c(-3,3),main="kay")
barplot(mean_tab_z["B2",],add=T,col=alpha("green",0.5),yaxt="n",xaxt="n")
lines(y=c(2,2),x=c(min(x),max(x)),lty=1,lwd=0.5)
lines(y=c(1,1),x=c(min(x),max(x)),lty=2,lwd=0.5)
lines(y=c(-1,-1),x=c(min(x),max(x)),lty=2,lwd=0.5)
lines(y=c(-2,-2),x=c(min(x),max(x)),lty=1,lwd=0.5)
barplot(mean_tab_z["C2",],col=alpha("red",0.5),las=2,ylim=c(-3,3),main="Rluc")
barplot(mean_tab_z["D2",],add=T,col=alpha("green",0.5),yaxt="n",xaxt="n")
lines(y=c(2,2),x=c(min(x),max(x)),lty=1,lwd=0.5)
lines(y=c(1,1),x=c(min(x),max(x)),lty=2,lwd=0.5)
lines(y=c(-1,-1),x=c(min(x),max(x)),lty=2,lwd=0.5)
lines(y=c(-2,-2),x=c(min(x),max(x)),lty=1,lwd=0.5)
barplot(mean_tab_z["E2",],col=alpha("red",0.5),las=2,ylim=c(-3,3),main="Drk")
barplot(mean_tab_z["F2",],add=T,col=alpha("green",0.5),yaxt="n",xaxt="n")
lines(y=c(2,2),x=c(min(x),max(x)),lty=1,lwd=0.5)
lines(y=c(1,1),x=c(min(x),max(x)),lty=2,lwd=0.5)
lines(y=c(-1,-1),x=c(min(x),max(x)),lty=2,lwd=0.5)
lines(y=c(-2,-2),x=c(min(x),max(x)),lty=1,lwd=0.5)
barplot(mean_tab_z["G2",],col=alpha("red",0.5),las=2,ylim=c(-3,3),main="Pten")
barplot(mean_tab_z["H2",],add=T,col=alpha("green",0.5),yaxt="n",xaxt="n")
lines(y=c(2,2),x=c(min(x),max(x)),lty=1,lwd=0.5)
lines(y=c(1,1),x=c(min(x),max(x)),lty=2,lwd=0.5)
lines(y=c(-1,-1),x=c(min(x),max(x)),lty=2,lwd=0.5)
lines(y=c(-2,-2),x=c(min(x),max(x)),lty=1,lwd=0.5)
barplot(mean_tab_z["I2",],col=alpha("red",0.5),las=2,ylim=c(-3,3),main="FAK")
barplot(mean_tab_z["J2",],add=T,col=alpha("green",0.5),yaxt="n",xaxt="n")
lines(y=c(2,2),x=c(min(x),max(x)),lty=1,lwd=0.5)
lines(y=c(1,1),x=c(min(x),max(x)),lty=2,lwd=0.5)
lines(y=c(-1,-1),x=c(min(x),max(x)),lty=2,lwd=0.5)
lines(y=c(-2,-2),x=c(min(x),max(x)),lty=1,lwd=0.5)
barplot(mean_tab_z["K2",],col=alpha("red",0.5),las=2,ylim=c(-3,3),main="GFP")
barplot(mean_tab_z["L2",],add=T,col=alpha("green",0.5),yaxt="n",xaxt="n")
lines(y=c(2,2),x=c(min(x),max(x)),lty=1,lwd=0.5)
lines(y=c(1,1),x=c(min(x),max(x)),lty=2,lwd=0.5)
lines(y=c(-1,-1),x=c(min(x),max(x)),lty=2,lwd=0.5)
lines(y=c(-2,-2),x=c(min(x),max(x)),lty=1,lwd=0.5)
barplot(mean_tab_z["M2",],col=alpha("red",0.5),las=2,ylim=c(-3,3),main="msk")
barplot(mean_tab_z["N2",],add=T,col=alpha("green",0.5),yaxt="n",xaxt="n")
lines(y=c(2,2),x=c(min(x),max(x)),lty=1,lwd=0.5)
lines(y=c(1,1),x=c(min(x),max(x)),lty=2,lwd=0.5)
lines(y=c(-1,-1),x=c(min(x),max(x)),lty=2,lwd=0.5)
lines(y=c(-2,-2),x=c(min(x),max(x)),lty=1,lwd=0.5)
barplot(mean_tab_z["O2",],col=alpha("red",0.5),las=2,ylim=c(-3,3),main="ras")
barplot(mean_tab_z["P2",],add=T,col=alpha("green",0.5),yaxt="n",xaxt="n")
lines(y=c(2,2),x=c(min(x),max(x)),lty=1,lwd=0.5)
lines(y=c(1,1),x=c(min(x),max(x)),lty=2,lwd=0.5)
lines(y=c(-1,-1),x=c(min(x),max(x)),lty=2,lwd=0.5)
lines(y=c(-2,-2),x=c(min(x),max(x)),lty=1,lwd=0.5)
dev.off()