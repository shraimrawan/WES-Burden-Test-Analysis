#tap1=read.table("tap1.txt", header=TRUE, sep="\t")
alldata=read.table("snpeffdata_102219.txt", header=TRUE, sep="\t")
id=colnames(alldata)

groups=read.table("remicade_groups.txt", header=TRUE)
res=which(groups[,2]=="RR")
res_id=groups[res,1]
ind1=match(res_id, id)

nonres=which(groups[,2]=="RN")
nonres_id=groups[nonres,1]
ind2=match(nonres_id, id)



#serping1 gene
serp=which(alldata$V3=="SERPING1")
serpd=alldata[serp,]
serpd=unique(serpd)

vartype=unique(serpd$V1)
vartable=matrix(0,length(vartype),2)
rownames(vartable)=vartype
colnames(vartable)=c("nonresponder","responder")

for (val in 1:length(vartype)){
  v1=which(serpd[,1]==vartype[val])
  nres1=sum(as.numeric(as.matrix(serpd[v1,ind2])), na.rm = TRUE)
  res1=sum(as.numeric(as.matrix(serpd[v1,ind1])), na.rm = TRUE)
  r1=c(nres1,res1)
  vartable[val,]=r1
}

varimp=unique(serpd$V2)
varimpactable=matrix(0,length(varimp),2)
rownames(varimpactable)=varimp
colnames(varimpactable)=c("nonresponder","responder")

for (val in 1:length(varimp)){
  v1=which(serpd[,2]==varimp[val])
  nres1=sum(as.numeric(as.matrix(serpd[v1,ind2])), na.rm = TRUE)
  res1=sum(as.numeric(as.matrix(serpd[v1,ind1])), na.rm = TRUE)
  r1=c(nres1,res1)
  varimpactable[val,]=r1
}




#tap1 Gene
tap1=which(alldata$V3=="TAP1")
tap1d=alldata[tap1,]
tap1d=unique(tap1d)

vartype=unique(tap1d$V1)
vartable=matrix(0,length(vartype),2)
rownames(vartable)=vartype
colnames(vartable)=c("nonresponder","responder")

for (val in 1:length(vartype)){
  v1=which(tap1d[,1]==vartype[val])
  nres1=sum(as.matrix(tap1d[v1,ind2]), na.rm = TRUE)
  res1=sum(as.matrix(tap1d[v1,ind1]), na.rm = TRUE)
  r1=c(nres1,res1)
  vartable[val,]=r1
}

varimp=unique(tap1d$V2)
varimpactable=matrix(0,length(varimp),2)
rownames(varimpactable)=varimp
colnames(varimpactable)=c("nonresponder","responder")

for (val in 1:length(varimp)){
  v1=which(tap1d[,2]==varimp[val])
  nres1=sum(as.matrix(tap1d[v1,ind2]), na.rm = TRUE)
  res1=sum(as.matrix(tap1d[v1,ind1]), na.rm = TRUE)
  r1=c(nres1,res1)
  varimpactable[val,]=r1
}

#FCGR3A gene
fcg=which(alldata$V3=="FCGR3A")
fcgd=alldata[fcg,]
fcgd=unique(fcgd)

vartype=unique(fcgd$V1)
vartable=matrix(0,length(vartype),2)
rownames(vartable)=vartype
colnames(vartable)=c("nonresponder","responder")

for (val in 1:length(vartype)){
  v1=which(fcgd[,1]==vartype[val])
  nres1=sum(as.numeric(as.matrix(fcgd[v1,ind2])), na.rm = TRUE)
  res1=sum(as.numeric(as.matrix(fcgd[v1,ind1])), na.rm = TRUE)
  r1=c(nres1,res1)
  vartable[val,]=r1
}

varimp=unique(fcgd$V2)
varimpactable=matrix(0,length(varimp),2)
rownames(varimpactable)=varimp
colnames(varimpactable)=c("nonresponder","responder")

for (val in 1:length(varimp)){
  v1=which(fcgd[,2]==varimp[val])
  nres1=sum(as.numeric(as.matrix(fcgd[v1,ind2])), na.rm = TRUE)
  res1=sum(as.numeric(as.matrix(fcgd[v1,ind1])), na.rm = TRUE)
  r1=c(nres1,res1)
  varimpactable[val,]=r1
}

#ADAR gene
adar=which(alldata$V3=="ADAR")
adard=alldata[adar,]
adard=unique(adard)

vartype=unique(adard$V1)
vartable=matrix(0,length(vartype),2)
rownames(vartable)=vartype
colnames(vartable)=c("nonresponder","responder")

for (val in 1:length(vartype)){
  v1=which(adard[,1]==vartype[val])
  nres1=sum(as.numeric(as.matrix(adard[v1,ind2])), na.rm = TRUE)
  res1=sum(as.numeric(as.matrix(adard[v1,ind1])), na.rm = TRUE)
  r1=c(nres1,res1)
  vartable[val,]=r1
}

varimp=unique(adard$V2)
varimpactable=matrix(0,length(varimp),2)
rownames(varimpactable)=varimp
colnames(varimpactable)=c("nonresponder","responder")

for (val in 1:length(varimp)){
  v1=which(adard[,2]==varimp[val])
  nres1=sum(as.numeric(as.matrix(adard[v1,ind2])), na.rm = TRUE)
  res1=sum(as.numeric(as.matrix(adard[v1,ind1])), na.rm = TRUE)
  r1=c(nres1,res1)
  varimpactable[val,]=r1
}

#ITGAM gene
itgam=which(alldata$V3=="ITGAM")
itgamd=alldata[itgam,]
itgamd=unique(itgamd)

vartype=unique(itgamd$V1)
vartable=matrix(0,length(vartype),2)
rownames(vartable)=vartype
colnames(vartable)=c("nonresponder","responder")

for (val in 1:length(vartype)){
  v1=which(itgamd[,1]==vartype[val])
  nres1=sum(as.numeric(as.matrix(itgamd[v1,ind2])), na.rm = TRUE)
  res1=sum(as.numeric(as.matrix(itgamd[v1,ind1])), na.rm = TRUE)
  r1=c(nres1,res1)
  vartable[val,]=r1
}

varimp=unique(itgamd$V2)
varimpactable=matrix(0,length(varimp),2)
rownames(varimpactable)=varimp
colnames(varimpactable)=c("nonresponder","responder")

for (val in 1:length(varimp)){
  v1=which(itgamd[,2]==varimp[val])
  nres1=sum(as.numeric(as.matrix(itgamd[v1,ind2])), na.rm = TRUE)
  res1=sum(as.numeric(as.matrix(itgamd[v1,ind1])), na.rm = TRUE)
  r1=c(nres1,res1)
  varimpactable[val,]=r1
}

#BCL10 gene
bcl=which(alldata$V3=="BCL10")
bcld=alldata[bcl,]
bcld=unique(bcld)

vartype=unique(bcld$V1)
vartable=matrix(0,length(vartype),2)
rownames(vartable)=vartype
colnames(vartable)=c("nonresponder","responder")

for (val in 1:length(vartype)){
  v1=which(bcld[,1]==vartype[val])
  nres1=sum(as.numeric(as.matrix(bcld[v1,ind2])), na.rm = TRUE)
  res1=sum(as.numeric(as.matrix(bcld[v1,ind1])), na.rm = TRUE)
  r1=c(nres1,res1)
  vartable[val,]=r1
}

varimp=unique(bcld$V2)
varimpactable=matrix(0,length(varimp),2)
rownames(varimpactable)=varimp
colnames(varimpactable)=c("nonresponder","responder")

for (val in 1:length(varimp)){
  v1=which(bcld[,2]==varimp[val])
  nres1=sum(as.numeric(as.matrix(bcld[v1,ind2])), na.rm = TRUE)
  res1=sum(as.numeric(as.matrix(bcld[v1,ind1])), na.rm = TRUE)
  r1=c(nres1,res1)
  varimpactable[val,]=r1
}

#ATP6AP1
atp=which(alldata$V3=="ATP6AP1")
atpd=alldata[atp,]
atpd=unique(atpd)

vartype=unique(atpd$V1)
vartable=matrix(0,length(vartype),2)
rownames(vartable)=vartype
colnames(vartable)=c("nonresponder","responder")

for (val in 1:length(vartype)){
  v1=which(bcld[,1]==vartype[val])
  nres1=sum(as.numeric(as.matrix(bcld[v1,ind2])), na.rm = TRUE)
  res1=sum(as.numeric(as.matrix(bcld[v1,ind1])), na.rm = TRUE)
  r1=c(nres1,res1)
  vartable[val,]=r1
}

varimp=unique(atpd$V2)
varimpactable=matrix(0,length(varimp),2)
rownames(varimpactable)=varimp
colnames(varimpactable)=c("nonresponder","responder")

for (val in 1:length(varimp)){
  v1=which(bcld[,2]==varimp[val])
  nres1=sum(as.numeric(as.matrix(bcld[v1,ind2])), na.rm = TRUE)
  res1=sum(as.numeric(as.matrix(bcld[v1,ind1])), na.rm = TRUE)
  r1=c(nres1,res1)
  varimpactable[val,]=r1
}



