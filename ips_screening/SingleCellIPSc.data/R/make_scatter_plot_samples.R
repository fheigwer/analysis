make_scatter_plot_samples<-function(sample_name="TAF1",size=2000,x=7,y=1,z=6,alpha=0.4){
  data("all_samples_classifications",package="SingleCellIPSc.data") #k-means classes of all control objects
  data("all_samples_feature_data",package="SingleCellIPSc.data") #feature matrix of all control objects
  data("all_samples_annotation",package="SingleCellIPSc.data") #feature matrix of all control objects
  
  #define the column, so the axis names
  colnames(sample_data_complete)=c("nuclear area", 
                                   "nuclear extent", 
                                   "maximum Oct4 intensity", 
                                   "maximum Hoechst intensity", 
                                   "mean Hoechst intensity", 
                                   "mean Oct4 intensity", 
                                   "distance from first neighbor", 
                                   "number of neighbors")
  #find all wells which contain the wanted annotation  
  pos=names_complete[,1][grep(sample_name,names_complete[,3])]
  #find and collect the objects feature vecotrs and classification that correnpond to the given well
  posmatrix=c()
  posclust=c()  
  for(i in pos){
    poswell=c(which(clustering_complete[,1]==i))  
    posmatrix=rbind(posmatrix,sample_data_complete[which(clustering_complete[,1]==i),])
    posclust=rbind(posclust,clustering_complete[which(clustering_complete[,1]==i),])
  }
  #draw a random subsample of the whole population of cells
  if(size>length(posclust[,2])){size<-length(posclust[,2])}
  sample=sample(1:(length(posclust[,2])),size=size)
  
  #draw the data
  scatterplot3d(
    posmatrix[sample,x],
    posmatrix[sample,y],
    posmatrix[sample,z],
    mar=c(10,10,10,10),
    col.lab="black",
    col.axis="black",
    col.grid="white",
    box=FALSE,
    # pch=".",
    pch=19,
    xlab=colnames(sample_data_complete)[7],
    ylab=colnames(sample_data_complete)[1],
    zlab=colnames(sample_data_complete)[6],
    xlim=c(-0.1,0.5),
    ylim=c(0,0.8),
    zlim=c(0,0.6),
    #color=alpha("black",0.4),
    color=alpha(posclust[sample,2],alpha)
  )
}