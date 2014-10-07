make_and_save_kmeans<-function(){
  data("raw_data_ctrls_with_random_cells",package="SingleCellIPSc.data")
  #discard unwanted columns they have been carefully assesed in prior analysis
  drop_v <- c(1,2,4,5,58,59,67,62,66,16,17,18,25,26,27,8,10)  
  keep <- c(1,4,20,19,25,26,47,48)
  raw_data<-raw_data[,c(keep)]
  #scale all features, feature wise between 0 and 1
  for(y in 1:dim(raw_data)[2]){raw_data[,y]=(raw_data[,y]-min(raw_data[,y]))/(max(raw_data[,y])-min(raw_data[,y]))}
  #perform kmeans clustering and save its result to the defaulte folder
  raw_clust <-kmeans(raw_data,6, nstart = 30)
  save(raw_data,file="raw_controls.RData")
  save(raw_clust,file="kmeans_result.RData")
}

