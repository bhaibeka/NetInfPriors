run.analysis<-function(name.priors,name.data,maxpa=10,mdist=2,randn=1000,vec.mypriorsw=c(0,0.25,0.5,0.75,0.95,1),mythreshp=0.1,mypriors=NULL,index.mypriors=NULL,kdoverfit=T ){	
	load("new_data/pn_priors_counts.RData")
	
	if(name.priors=="GM"){
		load("new_data/GM_priors.RData")
		pn.priors.counts<-gm.priors.counts[[index.mypriors]]	
	}else{
		print("use PN priors")
	}
				
	goi<-goi[1:genen]
	
	saveres <- paste("saveres_",name.data,sep="")
	
	load("new_data/rubio2011_marray_frma_symbol.RData")
	data <- data[ , goi, drop=FALSE]
	annot <- annot[colnames(data), ,drop=FALSE]
	demo <- demo[rownames(data), ,drop=FALSE]

	if(name.data=="expo"){
		print("load expo data")
		load("new_data/expo_frma_symbol.RData") 
	}else if(name.data=="jorissen"){
		print("load jorissen data")
		load("new_data/jorissen_frma_symbol.RData") 
	}
	
	data <- data[ , goi, drop=FALSE]
	perts <- matrix(0, nrow=nrow(data), ncol=ncol(data), dimnames=dimnames(data))
	if(!kdoverfit){
		  for(i in 1:nrow(demo)) {
			   if(demo[i, "KD"] %in% colnames(perts)) { perts[i, demo[i, "KD"]] <- 1 }
		  }
	}
	
	if(!file.exists(saveres)) { system(sprintf("mkdir %s", saveres)) }
	print(paste("save folder:",saveres))

	load("new_data/marray_kd_genes_affected.RData")
	gaff.kd<-gaff.all.rnk
	
	
	print("generate mat.pv")
	mat.pv<-matrix(0,nr=length(goi.core),nc=genen,dimnames=list(paste(goi.core,".KD",sep=""),colnames(data)))
	for(i in 1:length(goi.core)){
		mat.pv[i,]<-gaff.kd[[i]][goi,"pv"]
		mat.pv[i,]<-p.adjust(mat.pv[i,],method="fdr")
	}
	colnames(data) <- goi<-colnames(perts)<-colnames(mat.pv) <- colnames(pn.priors.counts)<-rownames(pn.priors.counts) <- gsub(pattern="[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]", replacement=".", x=toupper(goi))
	
	for(npw in 1:length(vec.mypriorsw)){
		mypriorsw<-vec.mypriorsw[npw]
		print(paste("priors weight",mypriorsw))
		
		if(kdoverfit){
			   if(name.priors=="GM"){
				    myfile <- sprintf("%s/net0_GM%i_priorsw%i.RData", saveres,index.mypriors, round(mypriorsw*100))
			   }else{
				    myfile <- sprintf("%s/net0_PN_priorsw%i.RData", saveres, round(mypriorsw*100))			
			   }
			   print("infer network")
			   if(!file.exists(myfile)){			
				    net0 <- netinf2(data=data, perturbations=perts, priors=pn.priors.counts, priors.count=TRUE, priors.weight=mypriorsw, method="regrnet", maxparents=maxpa)
				    save(list="net0", file=myfile)	
				    print(paste("results saved in",myfile))
			   }else{
				    load(myfile)
			   }
		
		}else{
			for(kk in 1:length(goi.core)){
				if(name.priors=="GM"){
				    myfile_i <- sprintf("%s/net0_GM%i_priorsw%i_kd%s.RData", saveres,index.mypriors, round(mypriorsw*100), names(gaff.kd)[kk])
				}else{
				    myfile_i <- sprintf("%s/net0_PN_priorsw%i_kd%s.RData", saveres, round(mypriorsw*100), names(gaff.kd)[kk])
				}
				    
				message(sprintf("Knock down of %s", names(gaff.kd)[kk]))
						      
				## select the kds
				kdix <- perts[ ,names(gaff.kd)[kk]] == 1
					     
				if(!file.exists(myfile_i)){
					print("run computation")
					if(kdoverfit){
					     kdiix <- rep(FALSE, nrow(perts)) 
					}else{
					     kdiix <- kdix 
					}
					net0 <- netinf2(data=data[!kdiix, , drop=FALSE], perturbations=perts[!kdiix, , drop=FALSE], priors=pn.priors.counts, priors.count=TRUE, priors.weight=mypriorsw, method="regrnet", maxparents=maxpa)
					save(list="net0", file=myfile_i)
				}
			}
		}
		
		print("start validation")		
		print("generate random networks")
		if(kdoverfit){
			topo.random<-shuffle.nodes(net0$topology,randn)
		}
		for(kk in 1:length(goi.core)){
			   print(paste("validation KD",goi.core[kk]))
			   if(!kdoverfit){
				    if(name.priors=="GM"){
					     myfile_i<-sprintf("%s/net0_GM%i_priorsw%i_kd%s.RData", saveres,index.mypriors, round(mypriorsw*100), names(gaff.kd)[kk])
				    }else{
					     myfile_i <- sprintf("%s/net0_PN_priorsw%i_kd%s.RData", saveres, round(mypriorsw*100), names(gaff.kd)[kk])
				    }
				    load(myfile_i)
				    print(myfile_i)
    				    topo.random<-shuffle.nodes(net0$topology,randn)

			   }
			   if(name.priors=="GM"){
				    myfile_i2<- sprintf("%s/fscore_randomshuffle_GM%i_priorsw%i_kd%s.RData", saveres,index.mypriors, round(mypriorsw*100), names(gaff.kd)[kk])
			   }else{
				    myfile_i2 <- sprintf("%s/fscore_randomshuffle_PN_priorsw%i_kd%s.RData", saveres, round(mypriorsw*100), names(gaff.kd)[kk])
			   }

			   if(!file.exists(myfile_i2)){
						
				    res.random.KD.i.mdist<-get.all.genes.na.random(topo.random,mat.pv,maxdistance=mdist,pval.cut=mythreshp,mykd=kk)
				    res.random.fscore<-2*(res.random.KD.i.mdist$mat.tp)/(2*(res.random.KD.i.mdist$mat.tp)+(res.random.KD.i.mdist$mat.fp)+(res.random.KD.i.mdist$mat.fn))
			
				    tmp.score<-get.number.na(net0,mat.pv,mythreshp,maxdistance=mdist,mykd=kk)   
				    res.fscore<-2*tmp.score$cnt.tp/(2*tmp.score$cnt.tp+tmp.score$cnt.fp+tmp.score$cnt.fn	)
				    print(res.fscore)
				    save(res.fscore,res.random.fscore,file=myfile_i2)
			   }
		  }
	}
}

shuffle.nodes<-function(mynet,nrandn){	
### function generating nrand random networks by randomly shuffling the network's node names
	net.random<-array(NA,c(dim(mynet)[[1]],dim(mynet)[[2]],nrandn),dimnames=list(dimnames(mynet)[[1]],dimnames(mynet)[[2]],seq(1,nrandn)))
	x<-rownames(mynet)
	for(i in 1:nrandn){
		tmp<-mynet	
		mysort<-x[order(runif(length(x)))]
		dimnames(tmp)<-list(mysort,mysort)
		net.random[,,i]<-tmp[x,x]
	}
	return(net.random)	  
}

get.all.genes.na.random<-function(res,mat.pval,maxdistance=NULL,pval.cut,mykd){
### wrapper for get.number.na when an array of topologies is provided in res
	mat.tp<-mat.fp<-mat.fn<-array(0,c(dim(res)[3]),dimnames=list(seq(1,dim(res)[3])))
	
	for(j in 1:dim(res)[3]){
		print(paste("random topology number",j))
			tmp<-get.number.na(list("topology"=res[,,j]),mat.pval,pval.cut,maxdistance=maxdistance,mykd=mykd)
			mat.tp[j]<-tmp$cnt.tp
			mat.fp[j]<-tmp$cnt.fp
			mat.fn[j]<-tmp$cnt.fn
	}
	
	return(list("mat.tp"=mat.tp,"mat.fp"=mat.fp,"mat.fn"=mat.fn))
	
}
get.number.na<-function(res,mat.pval,thres,maxdistance=NULL,mykd){
### function evaluating the childhood of mykd up to a distance maxdistance
### classifies genes as true positives (tp), false positives (fp) and false negatives (fn)
	names.core<-c("CDK5","HRAS","MAP2K1","MAP2K2","MAPK1","MAPK3","NGFR","RAF1")
	
	hops<-get.hops.foradj(res$topology)
	diag(hops)<-NA
	
	names.aff<-names(which(mat.pval[mykd ,]<thres))
	cnt.aff<-length(names.aff)
	
	if(is.null(maxdistance)){
		cnt.fp<-sum(!is.na(hops[names.core[mykd],names(which(mat.pval[mykd,]>thres))]))
		cnt.tp<-sum(!is.na(hops[names.core[mykd] ,names.aff]))
		cnt.fn<-length(names.aff)-sum(!is.na(hops[names.core[mykd] ,names.aff]))
	}else{			
		ind<-which(!is.na(hops[names.core[mykd],names(which(mat.pval[mykd,]>=thres))]))
		cnt.fp<-length(which(hops[names.core[mykd],names(ind)]<=maxdistance))
		ind<-which(!is.na(hops[names.core[mykd],names((which(mat.pval[mykd,]<thres)))]))			
		cnt.tp<-length(which(hops[names.core[mykd],names(ind)]<=maxdistance))
		cnt.fn<-cnt.aff-cnt.tp
	}
			
	return(list("cnt.tp"=cnt.tp,"cnt.fn"=cnt.fn,"cnt.fp"=cnt.fp)) 
}
get.hops.foradj<-function(am){

	am[am == 0] <- NA
	hops <- t(allShortestPaths(t(am))[[1]])
	dimnames(hops) <- dimnames(am)
	return(hops)
}

extract.results<-function(base.dir,vec.mypriorsw=c(0,0.25,0.5,0.75,0.95,1), vec.data=c("kd","expo","jorissen")){
	 goi.core <- c("CDK5", "HRAS", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3", "NGFR", "RAF1")
	 
	 myres<-myres.pval<-array(NA,dim=c(8,6,3,8),dimnames=list(goi.core,vec.mypriorsw,vec.data,c("PN",paste("GM",seq(2,8),sep=""))))
	 
	 for(i in 1:8){
		  for(j in 1:6){
			   for(k in 1:3){
				    saveres<-paste(base.dir,vec.data[k],sep="")

				    for(l in 1:8){
					     if(l==1){
						   index.mypriors<-"PN"
					     }else{
						   index.mypriors<-paste("GM",l,sep="")
					     }
					     myfile<-paste(saveres,"/fscore_randomshuffle_",index.mypriors,"_priorsw",round(vec.mypriorsw[j]*100),"_kd",goi.core[i],".RData",sep="")
					     if(file.exists(myfile)){
						      load(myfile)
						      myres[i,j,k,l]<-res.fscore						      
						      myres.pval[i,j,k,l]<-sum(res.fscore<=res.random.fscore)/1000 
					     }
				    }
			}
		}
	}
	
		
	 return(list(myres,myres.pval))

}

### generate result plots and tables

result.heatmap<-function(myres,index.priorw=6,index.data=1,save.plot=F){
	 
	myres.fscore<-t(myres[[1]][,index.priorw,index.data,])
	myres.pval<-t(myres[[2]][,index.priorw,index.data,])
	
	tmp<-(myres.pval<0.05)+(myres.pval<0.1)
	
	tmp.heat<-heatmap(tmp)
	myorder1<-rownames(myres.fscore)[rev(tmp.heat$rowInd)]
	myorder2<-colnames(myres.fscore)[rev(tmp.heat$colInd)]

	myres.fscore<-myres.fscore[myorder1,myorder2]
	myres.pval<-myres.pval[myorder1,myorder2]
	mycols <-colorRampPalette(c("blue3","cyan","aquamarine","yellow","orange","red"))(101)
	
	if(save.plot==TRUE){
		pdf(paste("results/heatmap_shufflerandom.pdf",sep=""),width=4.5,height=5)
	}
	
	
	plot(0:8, 0:8, type = "n",axes=F,xlab="",ylab="")
	
	axis(1,labels=myorder2,at=seq(0.5,8,1),line=F,tick=F,las=2)
	axis(2,labels=rev(myorder1),at=seq(0.5,8,1),line=F,tick=F,las=1)
		
	myres.fscore=myres.fscore/max(myres.fscore)

	mysignif <-(myres.pval<0.05)+(myres.pval<0.1)
	myangle<-c(-45,0,45)
	mydensity<-c(35,35,80)

	for(j in 1:8){
		for (i in 1:8){
			polygon(c(j-1,j-1,j,j),c(8-i,8-i+1,8-i+1,8-i),density=mydensity[mysignif[i,j]+1],col=mycols[round(myres.fscore[i,j]*100)+1], angle=myangle[mysignif[i,j]+1],border=NA)
		}	
	}
	for(j in 0:8){
		abline(h=j,col="lightgrey",lwd=2)
		abline(v=j,col="lightgrey",lwd=2)
	}

	if(save.plot==TRUE){		
		dev.off()
	}
 }
 
 plot.bars.fscore.colbyrank<-function(myres,data.index=1,prior.index=3,save.plot=FALSE){
	res.fscore<-myres[[1]]

	if(save.plot==TRUE){
		pdf(paste("results/barplot_fscore_bysourcecol.pdf",sep=""),width=12,height=4)
	}
	goi.core <- c("CDK5", "HRAS", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3", "NGFR", "RAF1")
	mypriorsources<-dimnames(res.fscore)[[4]][-c(2,10)]

	par(mfrow=c(1,8))
	par(mar=c(3,4,2,0.1))
	
	myymax<-0.5

	mycol<-rainbow(9, s=0.5, v=0.9)[-1]

	for(i in 1:8){
		tmp<-sort(res.fscore[i,prior.index,data.index,],decreasing=T,index.return=T)
		barplot(tmp$x,main=paste(dimnames(res.fscore)[[1]][i],".KD",sep=""),col=mycol[tmp$ix],las=2,ylim=c(0,myymax),ylab=("fscore"))
	}
	legend("topright",paste(mypriorsources),fill=(mycol),cex=0.95,bty="n")
	if(save.plot==TRUE){	
		dev.off()
	}
}

 




