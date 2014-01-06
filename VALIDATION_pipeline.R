rm(list=ls())
#set.seed(54321)

library(predictionet)
library(e1071)

source("script.R")
saveres.res<-"results"
nrandom<-1000
## eight core genes from biocarta2007
goi.core <- c("CDK5", "HRAS", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3", "NGFR", "RAF1")
## number of genes
genen=339


################### GENERATING RESULTS #######################
### generate results with PN prior for the three data sets
run.analysis("PN","kd",kdoverfit=F,randn=nrandom)
run.analysis("PN","expo",randn=nrandom)
run.analysis("PN","jorissen",randn=nrandom)


### generate results with GM priors for the three data sets
for(i in 2:8){
	 run.analysis("GM","kd",index.mypriors=i,kdoverfit=F,randn=nrandom)
	 run.analysis("GM","expo",index.mypriors=i,randn=nrandom)
	 run.analysis("GM","jorissen",index.mypriors=i,randn=nrandom)
}


###################### POST PROCESSING #######################

if(!file.exists(saveres.res)) { system(sprintf("mkdir %s", saveres.res)) }

### store all results
myres<-extract.results("saveres_")

### generate Figure 2
result.heatmap(myres,save.plot=T)

### generate Figure 3
plot.bars.fscore.colbyrank(myres,save.plot=T)

### generate Table 3
names.priors<-c("PN",paste("GM",seq(2,8),sep=""))
data.list<- c("kd","expo","jorissen")
myranks<-matrix(NA,nrow=8,ncol=3,dimnames=list(goi.core,data.list))
for(i in 1:3){
	 tmp<-apply(myres[[1]][,3,i,],1,order,decreasing=T)
	 myranks[,i]<-names.priors[tmp[1,]]
}
save(myranks,file="results/best.prior.RData")
