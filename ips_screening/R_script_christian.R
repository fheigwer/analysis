 a=as.data.frame(read.table("stuff_more.txt",sep="\t"))
 a=as.matrix(a)
 score=as.numeric(a[,2])
  count=as.numeric(a[,1])
 score_spec=score[count==1]
 score_spec=score_spec[!is.na(score_spec)]
 score_un_spec=score[count>4]
 score_un_spec=score_un_spec[!is.na(score_un_spec)]
 
 count_un=count[count>4]
 count_un= count_un[!is.na(count_un)]
 plot(density(score_un_spec),col="red")
 lines(density(score_spec),col="green")
 t.test(score_un_spec, score_spec)