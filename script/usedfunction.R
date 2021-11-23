library(ggplot2)
library(Biostrings)
library(stringdist)
library(rlist)
library(dplyr) 
library(circlize)
library(RColorBrewer)
library(parallel)
library(pheatmap)
library(reshape2)
library(ggtree)
library(data.tree)

ReadFasta = function(filename){
  sv = read.table(filename)
  scarfull = DNAString(sv[2,1])
  return(c("scarfull" = scarfull))
}

ReadCutsite = function(segref,reftype=NULL){
    colnames(segref) = c("indx","start","end")   
    scar = NULL
    type = NULL
    if(is.null(reftype)){
      for (i in 2:nrow(segref)) {
        scar = c(scar,segref[i,]$start:segref[i,]$end)
        type = c(type,rep(segref[i,]$indx,(segref[i,]$end-segref[i,]$start)+1))
      }
      scarseg = data.frame("scar" = scar,"type" = as.character(type))   
    }else{
      endsite<-NA
      for (i in 2:nrow(segref)) {
        endsite<-c(endsite,segref[["end"]][i]+segref[["start"]][1])
      }
      endsite[nrow(segref)]<-segref[["end"]][1]
      scarseg = data.frame("scar" = c(1:endsite[nrow(segref)]),"type" = NA)
      #endsite<-is.na(endsite)
      for (i in rev(endsite)) {
        if(!is.na(i)){
          scarseg$type[1:i]<-which(endsite==i)-1
         }else{
          break
         }
      }      
    }
     return(scarseg)
}

align_to_range = function(p,s,cut){
  pn <- strsplit(p,"")[[1]]
  sn <- strsplit(s,"")[[1]]
  lenp <- length(pn)
  index <- 1
  i <- 0
  del_flag <- F
  in_flag <- F
  del_start<-c()
  del_end<-c()
  ins_start<-c()
  ins_width<-c()
  ins_editseq<-c()
  while(i < lenp){
    i <- i + 1
    if(sn[[i]] == '-'){
      if(!del_flag){
        del_flag <- T
        del_start<-c(del_start,index)
        width <- 1
      }
      else{
        width <- width + 1
      }
    }
    else{
      if(del_flag){
        del_flag <- F
        del_end<-c(del_end,index)
        #print(paste("del stop width", width))
      }
    }
    if(pn[[i]] == '-'){
      if(!in_flag){
        in_flag <- T
        ins_start<-c(ins_start,index)
        width <- 1
        
        editseq<-sn[[i]]
      }
      else{
        width <- width + 1
        
        editseq<-paste0(editseq,sn[[i]])
      }
    }
    else{
      if(in_flag){
        in_flag <- F
        ins_width<-c(ins_width,width)
        ins_editseq<-c(ins_editseq,editseq)
        #print(paste("in stop width", width))
      }
    }
    if(pn[[i]] != '-')
      index <- index + 1
  }
  if(del_flag){
    del_flag <- F
    del_end<-c(del_end,index)
  }
  if(in_flag){
    in_flag <- F
    ins_width<-c(ins_width,width)
    ins_editseq<-c(ins_editseq,editseq)
    #print(paste("in stop width", width))
  }
  
  ins_start = ins_start-cut
  ins_end = ins_start + ins_width
  del_start = del_start-cut
  del_end = del_end-cut
  
  ins<-IRanges(start = ins_start,end = ins_end)
  mcols(ins)$seq <- ins_editseq
  del<-IRanges(start = del_start,end = del_end)
  
  
  return(list("del" = del,"ins" = ins))
  
}
FindScar = function(data,scarfull,scar,align_score=NULL,type=NULL,indel.coverage=NULL,cln){
  mat  =  nucleotideSubstitutionMatrix(match = 1, mismatch = -3)
  testFunction  =  function (data_in) {
    return(tryCatch(data_in, error=function(e) "unknown"))
  }
  if(is.null(align_score)){
    align_score = length(scarfull)- 2*(scar["end"][1,] - scar["start"][1,])-6
  }else{
    align_score = align_score
  }
  if(is.null(type)){
    type="None"
  }else{
    type=type
  }
  #type="none"
  find_barcode<-function(data){
    #tycpe="none"
    s3<-DNAString(as.character(data))
    alig<-pairwiseAlignment(scarfull,s3,substitutionMatrix = mat,type="global-local",gapOpening = 6, gapExtension= 1)
    if(is.null(indel.coverage)){
      scarshort = subseq(as.character(scarfull[[1]]),scar["start"][1,],scar["end"][1,])
    }else{
      scarshort = as.character(scarfull[[1]])
    }
    if(score(alig)<=align_score){
      r_read<-"unknown"
      r_scar<-"unknown"
      del=NA
      ins=NA
    }else{
      scar_pos = matchPattern(scarshort,as.character(pattern(alig)),with.indels = T,max.mismatch = 10)
      r_read = as.character(subject(alig))
      if(length(scar_pos)!= 0){
        r_scar = testFunction(subseq(as.character(subject(alig)),start=start(scar_pos),end=end(scar_pos)))
      }else{
        r_scar = "unknown"
      }
      stopifnot(is(alig, "PairwiseAlignments"), length(data) == 1L)
      
      p <- as.character(alignedPattern(alig)[[1L]])
      s <- as.character(alignedSubject(alig)[[1L]])
      if(is.null(indel.coverage)){
        delins = align_to_range(p,s,scar["start"][1,])
      }else{
        delins = align_to_range(p,s,1)
      }
      del = delins$del
      ins = delins$ins
      if(TRUE %in% (del@start<0) | TRUE %in% (ins@start<0)){
        r_scar = "unknown"
      }
    }
    fin_dat<-data.frame(new.reads=r_read, scar.BC=r_scar,type=type)
    return(list("del" = del,"ins" = ins,fin_dat))
  }
  
  cl = makeCluster(cln)
  clusterEvalQ(cl,library(Biostrings))
  environment(align_score) <- .GlobalEnv
  environment(scarfull) <- .GlobalEnv
  environment(scar) <- .GlobalEnv
  environment(data) <- .GlobalEnv
  environment(mat) <- .GlobalEnv
  environment(type) <- .GlobalEnv
  environment(find_barcode) <- .GlobalEnv
  environment(testFunction) <- .GlobalEnv
  environment(indel.coverage) <- .GlobalEnv
  clusterExport(cl,c('align_score','scarfull','scar','data','mat','indel.coverage','find_barcode','testFunction',"align_to_range","type"),envir = environment())
  scar_BC = parLapply(cl,data$Read.Seq,find_barcode)
  stopCluster(cl)
  
  #output
  data_2<-do.call("rbind",sapply(scar_BC,function(x){return(x[3])}))
  data<-cbind(data,data_2)
  scar_BC<-scar_BC[data$scar.BC!="unknown"]
  data<-data[data$scar.BC!="unknown",]
  saveRDS(list(scar_BC,data),"reads_metadata.rds")
  write.table(data,"all_UMI_reads_scar_full.txt",quote=F,sep="\t",row.names=F)
  data$scar_f<-gsub("[-]", "",as.character(data$scar.BC))
  data_v1<-data
  data_v1<-data_v1[data_v1$scar.BC!="unknown",]
  write.table(data_v1,"UMI_reads_scar_full.txt",quote=F,sep="\t",row.names=F)
  saveRDS(scar_BC,"indel.rds")
  return(list("INDEL" = scar_BC,"Scar" = data_v1))
}


change_form_stat<-function(indel){
  indel<-indel[c(1,2)]
  indel<-unlist(indel)
  if(length(indel)==0){
    return("unknown")
  }else{
    ins<-data.frame(indel[2])
    #ins$seq=as.character(indel[[2]]@elementMetadata$seq)
    site_ins<-apply(ins,1,function(x){c(x[1]:x[2])})
    if(dim(ins)[1]==1){
      site_ins<-list(as.numeric(site_ins))
    }
    cutsite_ins<-lapply(site_ins,function(x){unique(scarref$type[scarref$scar %in% x])})
    tag_ins<-apply(ins,1,function(x){paste0(x[3],"I+",x[1],x[4])})
    tag_ins<-lapply(tag_ins,function(x){rep(x,length(cutsite_ins[[which(tag_ins==x)]]))})
    tag_ins<-unlist(tag_ins)
    cutsite_ins<-unlist(cutsite_ins)
    del<-data.frame(indel[1])
    site_del<-apply(del,1,function(x){c(x[1]:x[2])})
    if(dim(del)[1]==1){
      site_del<-list(as.numeric( site_del))
    }
    cutsite_del<-lapply(site_del,function(x){unique(scarref$type[scarref$scar %in% x])})
    tag_del<-apply(del,1,function(x){paste0((x[3]-1),"D+",x[1])})
    tag_del<-lapply(tag_del,function(x){rep(x,length(cutsite_del[[which(tag_del==x)]]))})
    tag_del<-unlist(tag_del)
    cutsite_del<-unlist(cutsite_del)
    tag<-c(tag_del,tag_ins)
    cutsite<-c(cutsite_del,cutsite_ins)
    tag_all<-rep("NONE",length(unique(scarref$type)))
    if(length(tag)==0){
      return(paste(tag_all,collapse = "_"))
    }else{
      for(x in c(1:length(tag))){
        if(tag_all[as.numeric(cutsite[x])]=="NONE"){
          tag_all[as.numeric(cutsite[x])]<-tag[x]
        }else{
          tag_all[as.numeric(cutsite[x])]<-paste(tag_all[as.numeric(cutsite[x])],tag[x],sep="&")
        }
      }
    }
    return(paste(tag_all,collapse = "_"))
  }
}
INDELChangeForm = function(scarinfo,scarref,cln){
  cl<-makeCluster(cln)
  environment(change_form_stat) <- .GlobalEnv
  environment(scarinfo) <- .GlobalEnv
  environment(scarref) <- .GlobalEnv
  clusterExport(cl,c('scarinfo','change_form_stat',"scarref"), envir = environment())
  scar_form_p<-parLapply(cl,scarinfo$INDEL,change_form_stat)
  stopCluster(cl)
  scar_form<-unlist(scar_form_p)
  scar_form<-gsub(" ","",scar_form)
  data<-scarinfo$Scar 
  data$scar_form<-scar_form
  write.csv(data,"indel_pattern.csv",quote=F,row.names = F)
  return(data)  
}


INDELIdents = function(scarinfo,scarref,scarfull,scar,method.use=NULL,indel.coverage=NULL,cln){
  data<-scarinfo$Scar
  Cell.BC<-data.frame(table(data$Cell.BC))
  Cell.BC<-Cell.BC[Cell.BC$Freq>1,]
  data_1<-data[data$Cell.BC %in% Cell.BC$Var1,]
  data_1<-data_1[,c("Cell.BC","UMI","scar_f","scar_form")]
  mat  =  nucleotideSubstitutionMatrix(match = 1, mismatch = -3)
  # x = as.character(Cell.BC$Var1)[2]
  # dat =data_1
  max_reads_stat = function(x,dat,method.use=NULL){
    temreads = dat[dat$Cell.BC==x,]
    read_data = data.frame(table(as.character(temreads$scar_f)))
    read_data = read_data[order(-read_data$Freq),]
    scar_data = data.frame(table(as.character(temreads$scar_form)))
    scar_data = scar_data[order(-scar_data$Freq),]
    if(scar_data$Freq[1]==1){
      #can't define
      pattern="unkown"
      reads_num=0
      reads_pro=0
      umi_num=0
      umi_pro=0
      #UMI=0
      del=NA
      ins=NA
    }else{
      if(method.use=="consensus"){
        scarstrdist = stringdistmatrix(as.character(read_data$Var1),as.character(read_data$Var1))
        scarindex = which(apply(scarstrdist,1,function(x){sum(read_data$Freq[which(x<9)])>(sum(read_data$Freq)/3)}))
        if(length(scarindex) > 0){
          fin_read = consensusString(temreads$scar_f[temreads$scar_f %in% as.character(read_data$Var1)[scarindex]])
          reads_pro = round(length(which(temreads$scar_f %in%as.character(read_data$Var1)[scarindex]))/sum(read_data$Freq),4)
          reads_num = length(which(temreads$scar_f %in%as.character(read_data$Var1)[scarindex]))
          umi_num = length(unique(temreads$UMI[which(temreads$scar_f %in%as.character(read_data$Var1)[scarindex])]))
          umi_pro = round(umi_num/length(unique(temreads$UMI)),4)
          fin_read_cons = gsub("\\?","",fin_read)
          s1 = DNAString(fin_read_cons)
          aligc = pairwiseAlignment(scarfull,s1,substitutionMatrix = mat,type="global-local",gapOpening = 6, gapExtension = 1)
          stopifnot(is(aligc, "PairwiseAlignments"), length(x) == 1L)
          p <- as.character(alignedPattern(aligc)[[1L]])
          s.cons <- as.character(alignedSubject(aligc)[[1L]])
          if(is.null(indel.coverage)){
            indel = align_to_range(p,s.cons,scar["start"][1,])
          }else{
            indel = align_to_range(p,s.cons,1)
          }
          pattern =  change_form_stat(indel)
        }else{
          pattern = "unknown"
          reads_num=0
          reads_pro=0
          umi_num=0
          umi_pro=0
          #UMI=0
          #del=NA
          #ins=NA
          indel=NULL
        }
        #print("consensus")
      }else if(method.use=="reads.num"){
        pattern = as.character(scar_data$Var1[1])
        reads_num = scar_data$Freq[1]
        reads_pro = round(reads_num/sum(scar_data$Freq),4)
        umi_num = length(unique(temreads$UMI[which(temreads$scar_form %in% pattern)]))
        umi_pro = round(umi_num/length(unique(temreads$UMI)),4)
        #print("reads.num")
      }else{
        #method.use=="umi.num"
        read_data_umi = unique(temreads)
        read_data_umi = data.frame(table(read_data_umi$scar_form))
        read_data_umi = read_data_umi[order(-read_data_umi$Freq),]
        pattern = as.character(read_data_umi$Var1[1])
        umi_num = read_data_umi$Freq[1]
        umi_pro = round(umi_num/length(unique(temreads)$UMI),4)
        reads_num = scar_data[which(scar_data$Var1 == pattern),"Freq"]
        reads_pro = round(reads_num/sum(scar_data$Freq),4)
        #print("umi.num")
      }
    }
    fin_line = data.frame("Cell.BC" = x,
                          "pattern" = pattern,
                          "reads_num" = reads_num,
                          "reads_pro" = reads_pro,
                          "umi_num" = umi_num,
                          "umi_pro" = umi_pro,
                          stringsAsFactors = F)
    if(method.use=="consensus"){
      return(list(indel,fin_line))
    }else{
      return(fin_line)
    }
  }
  cl<-makeCluster(cln)
  clusterEvalQ(cl,library(Biostrings))
  clusterEvalQ(cl,library(stringdist))
  environment(data_1) <- .GlobalEnv
  environment(max_reads_stat) <- .GlobalEnv
  environment(Cell.BC) <- .GlobalEnv
  environment(scarref) <- .GlobalEnv
  environment(align_to_range) <- .GlobalEnv
  environment(change_form_stat) <- .GlobalEnv
  environment(scarfull) <- .GlobalEnv
  environment(mat) <- .GlobalEnv
  environment(scar) <- .GlobalEnv
  environment(method.use) <- .GlobalEnv
  environment(indel.coverage) <- .GlobalEnv
  clusterExport(cl,c('data_1','indel.coverage','method.use','max_reads_stat','Cell.BC',"align_to_range","scarref","change_form_stat","scarfull","mat","scar"), envir = environment())
  data_con<-parLapply(cl,as.character(Cell.BC$Var1),function(x)tryCatch(max_reads_stat(x,dat=data_1,method.use=method.use),error=function(e) NULL))
  stopCluster(cl)
  data_con<-data_con[!sapply(data_con,is.null)]
  
  if(method.use=="consensus"){
    INDEL_ranges_man <-list()
    for(i in seq(length(data_con))){
      INDEL_ranges_man <- c(INDEL_ranges_man,list(data_con[[i]][[1]]))
    }
    data_con_sub<-data.frame(Cell.BC=NA,pattern=NA,reads_num=NA, reads_pro=NA, umi_num=NA, umi_pro=NA)
    data_con_sub<-data_con_sub[-1,]
    for(i in seq(length(data_con))){
      data_con_sub<-rbind(data_con_sub,data_con[[i]][[2]])
    }
    data_con<-data_con_sub
    INDEL_ranges_man<-INDEL_ranges_man[data_con$pattern!="unkown"]
    data_con<-data_con[data_con$pattern!="unkown",]
    write.csv(data_con,"final_scarform.csv",quote=F,row.names = F)
    saveRDS(INDEL_ranges_man,"indel_ident.rds")
  }else{
    data_con<-do.call("rbind",data_con)
    data_con<-data_con[data_con$pattern!="unkown",]
    write.csv(data_con,"final_scarform.csv",quote=F,row.names = F)

    INDEL_ranges <- scarinfo$INDEL
    #INDEL_ranges <-INDEL_ranges[data$Cell.BC %in% Cell.BC$Var1]
    INDEL_ranges_man<-list()
    for(scarform in data_con$pattern){
      index=which(data$scar_form==scarform)[1]
      INDEL_ranges_man<-c(INDEL_ranges_man,list(INDEL_ranges[[index]][c(1,2)]))
    }
    saveRDS(INDEL_ranges_man,"indel_ident.rds")
  }
  return(list(indel=INDEL_ranges_man,info=data_con))
}

                                                                        
TagDataProcess = function(data,Cells=NULL,prefix=NULL){  
  TagStat = function(x) {
    x = as.character(x)
    umi = x[5]
    CB = x[1]
    x = x[2]
    x = unlist(strsplit(x,"_|&"))
    x = x[!x%in%c("NONE")]
    x = unique(x)
    if(length(x) == 0){
      return(NA)
    }else{
      return(data.frame(Cell.BC = CB,Reads_num = umi,Tag = x))
    }
  } 
 
  #tag = NULL
  #common.CB = NULL
  #if((dim(data)[1]>1)){
   # for (i in 1:(dim(data)[1]-1)) {
    #  common.CB = c(common.CB,intersect(data$Cell.BC[i], data$Cell.BC[i+1]))
    #}
  #}else{
    #common.CB = data$Cell.BC[1]
  #}
  
  #for (i in 1:dim(data)[1]) {
    if(!is.null(Cells)){
      data=data[data$Cell.BC %in% Cells$Cell.BC,]
    }
    #data[[i]]=data[[i]][data[[i]]$Cell.BC %in% common.CB,]
    tagi = apply(data,1,TagStat)
    tagi = do.call("rbind",tagi)
    tagi = na.omit(tagi)
    tabi = data.frame(table(tagi$Tag)/length(as.character(unique(tagi$Cell.BC))))
    #black list filter
    if(!is.null(prefix)){
      tagi$Tag = paste(prefix[i], tagi$Tag, sep = "")
    }else{
      tagi$Tag=tagi$Tag
    }
    #tag = rbind(tag,tagi)
  #}
  return(tagi)
}
                                                                        

BuildTagTree = function(tag,Cells=NULL){
  #tag stat
  Tag = data.frame(table(tag$Tag))
  tag$Cell.num = Tag$Freq[match(tag$Tag,Tag$Var1)]
  tag_tab = acast(tag,Cell.BC~Tag)
  tag_tab[!is.na(tag_tab)] = 1
  tag_tab[is.na(tag_tab)] = 0
  #tag integrate 
  cell_tab = data.frame(table(tag$Cell.BC))
  cell_tab = cell_tab[order(-cell_tab$Freq),]
  tags_all = lapply(as.character(cell_tab$Var1),function(x){sort(as.character(tag$Tag[tag$Cell.BC == x]))})
  tags_paste = sapply(tags_all,function(x){paste(x,collapse = "_")})
  tags_tab = data.frame(table(tags_paste))
  tags_tab$num = unlist(lapply(as.character(tags_tab$tags_paste),function(x){length(strsplit(x,split = "_")[[1]])}))
  tags_tab = tags_tab[order(-tags_tab$num,-tags_tab$Freq),]
  tags_uni = lapply(as.character(tags_tab$tags_paste),function(x){strsplit(x,split = "_")[[1]]})
  Tag_1 = Tag[Tag$Var1 %in% unlist(tags_uni[tags_tab$num == 1]),]
  Tag_1 = Tag_1[order(-Tag_1$Freq),]
  tags_uni[tags_tab$num == 1] = as.list(as.character(Tag_1$Var1))
  tags_tab[tags_tab$num == 1,] = tags_tab[tags_tab$num == 1,][match(as.character(Tag_1$Var1),
                                                                    as.character(tags_tab$tags_paste[tags_tab$num==1])),]
  
  #node build
  cluster_stat = function(i){
    x = tags_uni[[i]]
    n = tags_tab$num[i]
    tags_belone = NA
    for(y_ind in which(tags_tab$num<n)){
      y=tags_uni[[y_ind]]
      if(length(intersect(x,y))>0 & length(setdiff(y,x))==0){
        tags_belone = y_ind
        break
      }else{
        next
      }
    }
    return(tags_belone)
  }
  
  belons = sapply(which(tags_tab$num > 1),cluster_stat)
  belons[tags_tab$num == 1]= NA
  
  Tags = as.list(which(tags_tab$num == 1))
  names(Tags) = which(tags_tab$num == 1)
  nodes = as.list(c(1:length(tags_tab$num)))
  
  for(i in c(1:length(tags_tab$num))){
    n = belons[[i]]
    while(!is.na(n)){
      nodes[[i]] = c(n,nodes[[i]])
      n = belons[[n]]
    }
  }
  nodes_len = sapply(nodes,length)
  
  for(i in c(1:length(tags_tab$num))){
    nodes[[i]] = as.character(tags_tab$tags_paste)[nodes[[i]]]
    if(length(nodes[[i]])<max(nodes_len)){
      nodes[[i]][(length(nodes[[i]])+1):max(nodes_len)] = NA
    }else{
      next
    }
  }
  nodes=data.frame(do.call("rbind",nodes))
  names(nodes) = paste("N",as.character(1:ncol(nodes)),sep = "")
  nodes$pathString = do.call(paste, c("N0",nodes, sep="/"))
  nodes$pathString = gsub("/NA","",nodes$pathString)
  
  #save tree figure and rds
  population = as.Node(nodes)
  saveNetwork(diagonalNetwork(ToListExplicit(population, unname = TRUE),
                              margin = 10,fontSize = 8,width=15*dim(nodes)[2] ,height = 30*dim(nodes)[1]),
              file = "tree.html")
  
  saveRDS(population,"tree.rds")
  
  #save celltype tab
  cell_tab$tags = tags_paste
  if(!is.na(Cells)){
    cell_tab$celltype = Cells$Cell.type[match(as.character(cell_tab$Var1),Cells$Cell.BC)]
  }
  write.csv(cell_tab,"cell_tab.csv",row.names = F,quote = F)
  return(list(population,cell_tab))
  
}
                                                                        


