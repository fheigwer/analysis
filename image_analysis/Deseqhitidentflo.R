DESeqhitidentflo=function(list.of.untreated.genes,list.of.treated.genes,tophits){
  
  number.of.treated.genessets=length(list.of.treated.genes)
  number.of.untreated.genessets=length(list.of.untreated.genes)
  number.of.genes=length(list.of.treated.genes[[1]]$gene)
  genenames=as.character(list.of.treated.genes[[1]]$gene)
    
  untreatedmatrix=c()
  condition=c()
  names=c()
  for(i in 1:number.of.untreated.genessets){
    untreatedmatrix=cbind(untreatedmatrix,as.integer((list.of.untreated.genes[[i]]$fullmatch+list.of.untreated.genes[[i]]$seedmatch)))
    colnames(untreatedmatrix)[i]=paste("untreatedgenes",i,sep="")
    condition=c(condition,"untreated") 
    names=c(names,paste("untreatedgenes",i,sep=""))  
  }
  treatedmatrix=c()
  for(i in 1:number.of.treated.genessets){
    treatedmatrix=cbind(treatedmatrix,as.integer((list.of.treated.genes[[i]]$fullmatch+list.of.treated.genes[[i]]$seedmatch)))
    colnames(treatedmatrix)[i]=paste("treatedgenes",(i),sep="")
    condition=c(condition,"treated")
    names=c(names,paste("treatedgenes",(i),sep=""))
  }
  dfnew=as.data.frame(cbind(untreatedmatrix,treatedmatrix),row.names=genenames)             
  
 # condition = mymeta$condition 
  cds = newCountDataSet( dfnew, condition ) 
  cds = estimateSizeFactors(cds)                                                                           # perform DESeq analysis
  sizeFactors(cds)
  cds = estimateDispersions(cds,method="per-condition",fitType="local")                                                                        
  res = nbinomTest(cds,"untreated","treated")
  
  
  resall=res[order(res$pval),]
  resdec=res[(res$log2FoldChange)<0,]
  resdec=resdec[order(resdec$pval),]
  resinc=res[(res$log2FoldChange)>0,]
  resinc=resinc[order(resinc$pval),]                                                                       # generate list with top 30 significant hits (increased / decreased / both)
  resallSig=resall[ 1:tophits, ] 
  resdecSig=resdec[1:tophits,]
  resincSig=resinc[1:tophits,]
  allout=list()
  allout[["all"]]=resallSig
  allout[["decreased"]]=resdecSig
  allout[["increased"]]=resincSig
  
  print(resallSig[,c(1,6,7)])
  
  return(allout)
}
