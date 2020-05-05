rm (list = ls())
library(vcfR)
data=read.vcfR("out.vcf")
info=data@fix[,8]
r=1:length(info)
#finalmat=matrix(0,,5)
id=colnames(data@gt)
id=id[-1]
finalmat=matrix(0,length(info),5)
ind=matrix(0,length(info))


for (val in r ){
  a=unlist(strsplit(info[val],";"))
  rn=grep("ANN",a)
  if (length(rn)==0) {
    next
  }
  else { b=unlist(strsplit(a[rn],"\\|"))
  d=c(b[2],b[3],b[4],b[7],b[8])
  finalmat[val,]=d
  ind[val]=val
  }
}

var=as.data.frame(finalmat)
var$chrom=data@fix[,1]
var$position=data@fix[,2]
var$ref=data@fix[,4]
var$alt=data@fix[,5]


genotype=data@gt

get_genotype = function(x){
  y=unlist(strsplit(x,":"))[1]
}

genotype2 = genotype[,-1]
test2 = as.data.frame(matrix(0,nrow(genotype2),ncol(genotype2)))

for (i in 1:ncol(genotype2)){
  geno = apply(as.data.frame(genotype2[,i]),1,function(x) {unlist(strsplit(x,":"))[1]})
  test2[,i] = as.data.frame(geno)
}


write.table(var, file = "C:/Users/shraimr/Documents/snpeff.xls", sep = "\t", row.names=FALSE, quote = FALSE)
write.table(test2, file = "C:/Users/shraimr/Documents/snpeff2.xlsx", sep = "\t", row.names=FALSE, quote = FALSE)

gt=read.table("snpeff2.xlsx", header = TRUE)
colnames(test2)=id
test=as.matrix(test2)
test=gsub("0/0","0",test)
test=gsub("0/1","1",test)
test=gsub("1/1","2", test)
test=gsub("./.","NA",test)
test=gsub("\\.","NA",test)

finaldata=cbind(var, test)

write.table(finaldata, file = "C:/Users/shraimr/Documents/Remicade Project/SNPeff/snpeffdata_102219.txt", sep = "\t", row.names=FALSE, quote = FALSE)

