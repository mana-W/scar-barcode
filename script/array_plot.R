
library(IRanges)

x = "24D+3_3I+70AAA_12D+100"

plotarray <- function(x,endsite){
  x <- unique(unlist(strsplit(x,"_|&")))
  x <- x[x!="NONE"]
  xd <- x[grep("D",x)]
  xi <- x[grep("I",x)]
  
  drange <- function(i){
    is <- as.numeric(sapply(i,function(y){strsplit(y,"D[+]")[[1]][2]}))
    iw <- as.numeric(sapply(i,function(y){strsplit(y,"D[+]")[[1]][1]}))
    return(IRanges(start = is,width = iw,type = "Del"))
  }
  xd <- drange(xd)
  irange <- function(i){
    is <- as.numeric(gsub("A|T|C|G", "", unlist(sapply(i,function(y){strsplit(y,"I[+]")[[1]][2]}))))
    iw <- as.numeric(sapply(i,function(y){strsplit(y,"I[+]")[[1]][1]}))
    return(IRanges(start = is,width = iw,type = "Ins"))
  }
  xi <- irange(xi)
  x <- c(xd,xi)
  all <- IRanges(start = 0,end = endsite,type = "All")
  back <- setdiff(all,x)
  mcols(back)$type <- "All"
  x <- c(x, back)
  coldata <- data.frame(col = c("red","blue","gray77"),type=c("Del","Ins","All"))
  cols <- coldata$col[match(x@elementMetadata$type,coldata$type)]
  plot(0:endsite, 0:endsite, xlim = c(0,endsite), ylim = c (0,3), type = "n")
  rect(start(x), 1, end(x), 2,col = cols,border = NA)
}

plotarray(x,180)
