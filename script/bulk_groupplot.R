#data file outputpath

library(ggplot2)
library(ggpubr)

args<-commandArgs(T)

pathdata<-args[1]
data<-read.table(pathdata)
outputpath<-args[2]

readsnum=c()
for(file in data$V1){
subnum=as.numeric(system(paste0("wc -l ",file,"| awk '{print$1}'"),intern = TRUE))
readsnum=c(readsnum,subnum)
}
idlen=sort(readsnum)[which(sort(readsnum)>1000)[1]]

pro<-data.frame(Proportion=NA,group=NA,type=NA)
pro<-pro[-1,]
div<-data.frame(Diversity=NA,group=NA,type=NA)
div<-div[-1,]
cuttable<-pro
cutted<-pro
names(cutted)[1]<-"Num"
uncutted<-cutted
len<-data.frame(del_length=NA,group=NA,type=NA)
len<-len[-1,]

del_num_cal<-function(array){
	indels<-unique(strsplit(array,"_|&")[[1]])
	if(identical(indels,"NONE")){
	 return(0)
	}else{
	del=unlist(lapply(indels[grep("D",indels)],function(x){as.numeric(strsplit(x,"D")[[1]][1])}))
	ins=unlist(lapply(indels[grep("I",indels)],function(x){as.numeric(strsplit(x,"I")[[1]][1])}))
	return((sum(del)-sum(ins)))
	}
}

for(i in seq(data$V1)){
	reads<-read.table(data$V1[i])$V1
	dataid<-data$V3[i]
	type=data$V2[i]
	fullarray=data$V4[i]
	fullarray_num=length(strsplit(data$V4[i],"_")[[1]])
	edit_num<-length(reads[reads!=fullarray])
	pro<-rbind(pro,data.frame(Proportion=edit_num/length(reads),
								group=dataid,
								type=type))
	cuttable<-rbind(cuttable,data.frame(Proportion=length(grep("NONE",reads))/length(reads),
										group=dataid,
										type=type))
	cutted<-rbind(cutted,data.frame(Num=unlist(lapply(reads,
												      function(x,y){
														y-length(which(unlist(strsplit(x,"_"))=="NONE"))
													  },
													  y=fullarray_num)),
									group=dataid,
									type=type))	
    uncutted<-rbind(uncutted,data.frame(Num=unlist(lapply(reads,
												      function(x){
														length(which(unlist(strsplit(x,"_"))=="NONE"))
													  })),
									group=dataid,
									type=type))
	len<-rbind(len,data.frame(del_length=unlist(lapply(reads,del_num_cal)),
							group=dataid,
							type=type))
	if(length(reads)<idlen){
		uselen=length(reads)
	}else{
		uselen=idlen
	}
	reads<-reads[sample(seq(reads),uselen)]
	div<-rbind(div,data.frame(Diversity=length(unique(as.character(reads))),
								group=dataid,
								type=type))
}
div$Diversity_scale<-(div$Diversity-min(div$Diversity))/(max(div$Diversity)-min(div$Diversity))

pro$group<-factor(pro$group,levels=unique(data$V3))
div$group<-factor(div$group,levels=unique(data$V3))
cuttable$group<-factor(cuttable$group,levels=unique(data$V3))
cutted$group<-factor(cutted$group,levels=unique(data$V3))
uncutted$group<-factor(uncutted$group,levels = unique(data$V3))
len$group<-factor(len$group,levels = unique(data$V3))

pro$type<-factor(pro$type,levels=unique(data$V2))
div$type<-factor(div$type,levels=unique(data$V2))
cuttable$type<-factor(cuttable$type,levels=unique(data$V2))
cutted$type<-factor(cutted$type,levels=unique(data$V2))
uncutted$type<-factor(uncutted$type,levels = unique(data$V2))
len$type<-factor(len$type,levels = unique(data$V2))

summ<-list(div,pro,cuttable,cutted,uncutted,len)
saveRDS(summ,paste0(outputpath,"/summary.rds"))


#tbc
