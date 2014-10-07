make_scatter_plot_ctrls<-function(size=20000,x=7,y=1,z=6,alpha=0.4){
  data("Oct4_positive_classfication",package="SingleCellIPSc.data") #k-means classes of all positive control objects
  data("Oct4_positive_feature_matrix",package="SingleCellIPSc.data") #feature vectors of all positive control objects
  data("Oct4_negative_classfication",package="SingleCellIPSc.data") #k-means classes of all negative control objects
  data("Oct4_negative_feature_matrix",package="SingleCellIPSc.data") #feature vectors of all negative control objects
  
  #draw a random subsample of the whole population of cells
  if(size>length(posclust[,2])){size<-length(posclust[,2])}
  sample=sample(1:(length(posclust[,2])),size=size)
  par(mfrow=c(1,2))
  #draw the data
  scatterplot3d(
    posmatrix[sample,x],
    posmatrix[sample,y],
    posmatrix[sample,z],
    col.lab="black",
    col.axis="black",
    col.grid="white",
    box=FALSE,
    pch=".",
    #pch=19,
    xlab=colnames(negmatrix)[7],
    ylab=colnames(negmatrix)[1],
    zlab=colnames(negmatrix)[6],
    xlim=c(-0.1,0.5),
    ylim=c(0,0.8),
    zlim=c(0,0.6),
    #color=alpha("black",0.4),
    color=alpha(posclust[sample,2],alpha)
  )
  if(size>length(negclust[,2])){size<-length(negclust[,2])}
  sample=sample(1:(length(negclust[,2])),size=size)
  scatterplot3d(
    negmatrix[sample,x],
    negmatrix[sample,y],
    negmatrix[sample,z],
    col.lab="black",
    col.axis="black",
    col.grid="white",
    box=FALSE,
    pch=".",
    #pch=19,
    xlab=colnames(negmatrix)[7],
    ylab=colnames(negmatrix)[1],
    zlab=colnames(negmatrix)[6],
    xlim=c(-0.1,0.5),
    ylim=c(0,0.8),
    zlim=c(0,0.6),
    #color=alpha("black",0.4),
    color=alpha(negclust[sample,2],alpha)
  )
}