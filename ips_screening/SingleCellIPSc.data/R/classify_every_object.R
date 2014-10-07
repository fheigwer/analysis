classify_every_object<-function(plate="10301"){
  #load the different data from the package
  data("kmeans_result",package="SingleCellIPSc.data") #k-means classes of all control objects
  data("raw_controls",package="SingleCellIPSc.data") #feature matrix of all control objects
  data("raw_data_sample_10301",package="SingleCellIPSc.data") #feature matrix of all sample objects
  data("origin10301",package="SingleCellIPSc.data") #vector containing all origins of each object
  #perform a K-nearest-Neighbour classification on the basis of the formerly found k-means clustering
  clustering=cbind(origin,knn(raw_data, sample_data, raw_clust$cluster, k = 100, prob=F))
  save(clustering,file=paste(plate,"clustering.RData",sep="_"))
}