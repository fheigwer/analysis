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



