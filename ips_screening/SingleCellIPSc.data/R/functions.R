library(class)
library(scatterplot3d)
library(scales)

make_and_save_kmeans<-function(raw_data_ctrls_with_random_cells){
#discard unwanted columns they have been carefully assesed in prior analysis
drop_v <- c(1,2,4,5,58,59,67,62,66,16,17,18,25,26,27,8,10)  
keep <- c(1,4,20,19,25,26,47,48)
raw_data<-raw_data_ctrls_with_random_cells[,c(keep)]
#scale all features, feature wise between 0 and 1
for(y in 1:dim(raw_data)[2]){raw_data[,y]=(raw_data[,y]-min(raw_data[,y]))/(max(raw_data[,y])-min(raw_data[,y]))}
#perform kmeans clustering and save its result to the defaulte folder
raw_clust <-kmeans(raw_data,6, nstart = 30)
save(raw_data,file="raw_controls.RData")
save(raw_clust,file="kmeans_result.RData")
}


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
    poswellc(which(clustering_complete[,1]==i))  
    posmatrix=rbind(posmatrix,sample_data_complete[which(posclust[,1]==i),])
    posclust=rbind(posclust,clustering_complete[which(posclust[,1]==i),])
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

make_scatter_plot_ctrls<-function(sample_name="pos",size=20000,x=7,y=1,z=6,alpha=0.4){
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

make_distance_plot<-function(absolute=FALSE,relative=FALSE,feature=FALSE){
  posnew=c()
  negnew=c()
  data("neg_cluster_matrix10301",package="SingleCellIPSc.data")                 
  data("neg_feature_matrix10301",package="SingleCellIPSc.data")    
  data("pos_feature_matrix10301",package="SingleCellIPSc.data")    
  data("pos_cluster_matrix10301",package="SingleCellIPSc.data")
  for(i in unique(posclust[,1])){    
    if(absolute==TRUE){
      postats=c()
      for(j in 1:6){
        postats=c(postats,length(which(posclust[which(posclust[,1]==i),2]==j)))
      }
    }else if(relative==TRUE){
      postats=c()
      for(j in 1:6){
        postats=c(postats,length(which(posclust[which(posclust[,1]==i),2]==j)))
      }
      postats=postats/sum(postats)*100      
    }else if(feature==TRUE){
      postats=colMeans(posmatrix[which(posclust[,1]==i),c(3,6)])
    }
    posnew=rbind(posnew,postats)    
  }
  
  
  for(i in unique(negclust[,1])){    
    if(absolute==TRUE){
      negtats=c()
      for(j in 1:6){
        negtats=c(negtats,length(which(negclust[which(negclust[,1]==i),2]==j)))
      }
    }else if(relative==TRUE){
      negtats=c()
      for(j in 1:6){
        negtats=c(negtats,length(which(negclust[which(negclust[,1]==i),2]==j)))
      }
      negtats=negtats/sum(negtats)*100      
    }else if(feature==TRUE){
      negtats=colMeans(negmatrix[which(negclust[,1]==i),c(3,6)])
    }
    negnew=rbind(negnew,negtats)    
  }
  data("neg_cluster_matrix10302",package="SingleCellIPSc.data")  
  data("neg_feature_matrix10302",package="SingleCellIPSc.data")  
  data("pos_feature_matrix10302",package="SingleCellIPSc.data")  
  data("pos_cluster_matrix10302",package="SingleCellIPSc.data")   
  for(i in unique(posclust[,1])){    
    if(absolute==TRUE){
      postats=c()
      for(j in 1:6){
        postats=c(postats,length(which(posclust[which(posclust[,1]==i),2]==j)))
      }
    }else if(relative==TRUE){
      postats=c()
      for(j in 1:6){
        postats=c(postats,length(which(posclust[which(posclust[,1]==i),2]==j)))
      }
      postats=postats/sum(postats)*100      
    }else if(feature==TRUE){
      postats=colMeans(posmatrix[which(posclust[,1]==i),c(3,6)])
    }
    posnew=rbind(posnew,postats)    
  }
  
  
  for(i in unique(negclust[,1])){    
    if(absolute==TRUE){
      negtats=c()
      for(j in 1:6){
        negtats=c(negtats,length(which(negclust[which(negclust[,1]==i),2]==j)))
      }
    }else if(relative==TRUE){
      negtats=c()
      for(j in 1:6){
        negtats=c(negtats,length(which(negclust[which(negclust[,1]==i),2]==j)))
      }
      negtats=negtats/sum(negtats)*100      
    }else if(feature==TRUE){
      negtats=colMeans(negmatrix[which(negclust[,1]==i),c(3,6)])
    }
    negnew=rbind(negnew,negtats)    
  }
  data("neg_feature_matrix10303",package="SingleCellIPSc.data") 
  data("neg_cluster_matrix10303",package="SingleCellIPSc.data")  
  data("pos_cluster_matrix10303",package="SingleCellIPSc.data")               
  data("pos_feature_matrix10303",package="SingleCellIPSc.data")    
  
  for(i in unique(posclust[,1])){    
    if(absolute==TRUE){
      postats=c()
      for(j in 1:6){
        postats=c(postats,length(which(posclust[which(posclust[,1]==i),2]==j)))
      }
    }else if(relative==TRUE){
      postats=c()
      for(j in 1:6){
        postats=c(postats,length(which(posclust[which(posclust[,1]==i),2]==j)))
      }
      postats=postats/sum(postats)*100      
    }else if(feature==TRUE){
      postats=colMeans(posmatrix[which(posclust[,1]==i),c(3,6)])
    }
    posnew=rbind(posnew,postats)    
  }  
  for(i in unique(negclust[,1])){    
    if(absolute==TRUE){
      negtats=c()
      for(j in 1:6){
        negtats=c(negtats,length(which(negclust[which(negclust[,1]==i),2]==j)))
      }
    }else if(relative==TRUE){
      negtats=c()
      for(j in 1:6){
        negtats=c(negtats,length(which(negclust[which(negclust[,1]==i),2]==j)))
      }
      negtats=negtats/sum(negtats)*100      
    }else if(feature==TRUE){
      negtats=colMeans(negmatrix[which(negclust[,1]==i),c(3,6)])
    }
    negnew=rbind(negnew,negtats)    
  }
  negnew=colMeans(negnew)
  posnew=colMeans(posnew)
  samplenew=c()  
  data ("clustering_sample_10301",package="SingleCellIPSc.data") #k-means classes of all control objects     
  data ("wellnames10301",package="SingleCellIPSc.data") #k-means classes of all control objects
  data ("raw_data_sample_10301",package="SingleCellIPSc.data") #k-means classes of all control objects  
  bad=c(grep("pos",b[,3]),grep("neg",b[,3]),grep("empty",b[,3]))
  for(i in b[-bad,1]){    
      if(absolute==TRUE){
        sampletats=c()
        for(j in 1:6){
          sampletats=c(sampletats,length(which(clustering[which(clustering[,1]==i),2]==j)))
        }
      }else if(relative==TRUE){
        sampletats=c()
        for(j in 1:6){
          sampletats=c(sampletats,length(which(clustering[which(clustering[,1]==i),2]==j)))
        }
        sampletats=sampletats/sum(sampletats)*100      
      }else if(feature==TRUE){
        sampletats=colMeans(sample_data[which(clustering[,1]==i),c(3,6)])
      }
      samplenew=rbind(samplenew,sampletats)    
    }
  firstnames=b[-bad,3]
  
  data ("clustering_sample_10302",package="SingleCellIPSc.data") #k-means classes of all control objects  
  data ("raw_data_sample_10302",package="SingleCellIPSc.data") #k-means classes of all control objects  
  data ("wellnames10302",package="SingleCellIPSc.data") #k-means classes of all control objects  
  bad=c(grep("pos",b[,3]),grep("neg",b[,3]),grep("empty",b[,3]))
  for(i in b[-bad,1]){    
    if(absolute==TRUE){
      sampletats=c()
      for(j in 1:6){
        sampletats=c(sampletats,length(which(clustering[which(clustering[,1]==i),2]==j)))
      }
    }else if(relative==TRUE){
      sampletats=c()
      for(j in 1:6){
        sampletats=c(sampletats,length(which(clustering[which(clustering[,1]==i),2]==j)))
      }
      sampletats=sampletats/sum(sampletats)*100      
    }else if(feature==TRUE){
      sampletats=colMeans(sample_data[which(clustering[,1]==i),c(3,6)])
    }
    samplenew=rbind(samplenew,sampletats)    
  }
  firstnames=c(firstnames,b[-bad,3])
  data ("wellnames10303",package="SingleCellIPSc.data") #k-means classes of all control objects
  data ("clustering_sample_10303",package="SingleCellIPSc.data") #k-means classes of all control objects         
  data ("raw_data_sample_10303",package="SingleCellIPSc.data") #k-means classes of all control objects
  bad=c(grep("pos",b[,3]),grep("neg",b[,3]),grep("empty",b[,3]))
  for(i in b[-bad,1]){    
    if(absolute==TRUE){
      sampletats=c()
      for(j in 1:6){
        sampletats=c(sampletats,length(which(clustering[which(clustering[,1]==i),2]==j)))
      }
    }else if(relative==TRUE){
      sampletats=c()
      for(j in 1:6){
        sampletats=c(sampletats,length(which(clustering[which(clustering[,1]==i),2]==j)))
      }
      sampletats=sampletats/sum(sampletats)*100      
    }else if(feature==TRUE){
      sampletats=colMeans(sample_data[which(clustering[,1]==i),c(3,6)])
    }
    samplenew=rbind(samplenew,sampletats)    
  }
  row.names(samplenew)=c(firstnames,b[-bad,3]) 
  
  dist_to_pos=c()
  for(i in 1:dim(samplenew)[1]){
    dist_to_pos=c(dist_to_pos,(sum((posnew-samplenew[i,])^2))^(1/2))
  }
  names(dist_to_pos)=rownames(samplenew)
  dist_to_neg=c()
  for(i in 1:dim(samplenew)[1]){
    dist_to_neg=c(dist_to_neg,(sum((negnew-samplenew[i,])^2))^(1/2))
  }
  names(dist_to_neg)=rownames(samplenew)
  novel_behav=intersect(which(log(dist_to_pos)>log(median(dist_to_pos))),which(log(dist_to_neg)>log(median(dist_to_neg))))
  forbidden=intersect(which(log(dist_to_pos)<log(median(dist_to_pos))),which(log(dist_to_neg)<log(median(dist_to_neg))))
  true_pos=intersect(which(log(dist_to_pos)<log(median(dist_to_pos))),which(log(dist_to_neg)>log(median(dist_to_neg))))
  true_neg=intersect(which(log(dist_to_pos)>log(median(dist_to_pos))),which(log(dist_to_neg)<log(median(dist_to_neg))))
  
  #load("novel.RData")
  #load("forbidden.RData")
  #load("true_neg.RData")
  #load("true_pos.RData")
  
  plot(log(dist_to_pos), 
       log(dist_to_neg),
       type="n",
       pch=19, 
       xlab="ln ( euclidean distance to positive controls)", 
       ylab="ln ( euclidean distance to negative controls)")
  
  lines(c(log(median(dist_to_pos)), log(median(dist_to_pos))),c(-10,max(dist_to_neg)),col="red")
  lines(c(-10, max(dist_to_pos)),c(log(median(dist_to_neg)), log(median(dist_to_neg))),col="red")
  
  points(log(dist_to_pos)[novel_behav],log(dist_to_neg)[novel_behav],pch=19,col="#a6611a")
  points(log(dist_to_pos)[forbidden],log(dist_to_neg)[forbidden],pch=19,col="#dfc27d")
  points(log(dist_to_pos)[true_pos],log(dist_to_neg)[true_pos],pch=19,col="#80cdc1")
  points(log(dist_to_pos)[true_neg],log(dist_to_neg)[true_neg],pch=19,col="#018571")
  points(log(dist_to_pos)[grep("TAF1$",names(dist_to_pos))],log(dist_to_neg)[grep("TAF1$",names(dist_to_pos))],col="red",pch="x")
  points(log(dist_to_pos)[grep("SMG1$",names(dist_to_pos))],log(dist_to_neg)[grep("SMG1$",names(dist_to_pos))],col="red",pch="o")
  points(log(dist_to_pos)[grep("MAPK3$",names(dist_to_pos))],log(dist_to_neg)[grep("MAPK3$",names(dist_to_pos))],col="red",pch="+")
}



