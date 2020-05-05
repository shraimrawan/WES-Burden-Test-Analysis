#The output of the burden test from pseq will be fed into this R script 
bt= read.table("burden102219.txt", stringsAsFactors = FALSE, header = TRUE)

r=1:nrow(bt)

finalmat=matrix(0,nrow(bt),2)
d=matrix(nrow = 0, ncol=3)

for (val in r) {
  a=unlist(strsplit(bt[val,1],";"))
  b=paste(a, collapse = "")
  c=as.numeric(unlist(strsplit(b,"/|\\(|)")))
  d=matrix(c, ncol=3, byrow = TRUE)
  nres=d[,1] %*% d[,3]
  res=d[,2] %*% d[,3]
  result=c(nres, res)
  finalmat[val,]=result
}

colnames(finalmat)=c("nonresponder","responder")


write.table(finalmat, file = "C:/Users/shraimr/Documents/hdvariantcounts_102219.txt", sep = "\t", row.names=FALSE, quote = FALSE)
