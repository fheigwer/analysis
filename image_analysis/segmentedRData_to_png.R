library(EBImage)
args=commandArgs(trailingOnly = TRUE)
options(warn=-1)
image=args[1]
identifier=image#gsub("(..).+$","\\1",image)
load(file=paste(identifier,"_segmented",".RData",sep=""))
writeImage(res, file=paste(identifier,"_segmented",".png",sep=""), type="png")