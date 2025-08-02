#getVisiumHDspotCol.R
#A function to extract average R, G, B color intensities of each bin from the immunofluorescence image paired with the Visium HD data.

library(Seurat)
library(ggplot2)
library(jpeg)
library(BiocParallel)
getVisiumHDspotCol_value <- function(X,vName){
  return(as.list(sapply(vName,function(x)return(get(x)(X)))))
}
getVisiumHDspotCol <- function(D,strImg,r=NULL,core=4){
  oriIMG <- jpeg::readJPEG(strImg)
  ## process
  H <- nrow(oriIMG)
  W <- ncol(oriIMG)
  RGB <- c("R","G","B")
  vName <- c("mean","median","sd")
  cInfo <- lapply(D@images,function(img){
    if(is.null(r)) r <- round(img@boundaries$centroids@radius)
    res <- bplapply(1:nrow(img@boundaries$centroids@coords),function(i){
      x <- round(img@boundaries$centroids@coords[i,"x"])
      y <- round(img@boundaries$centroids@coords[i,"y"])
      #if(x<1 || x>nrow(oriIMG) || y<1 || y>ncol(oriIMG)) return()
      selX <- max(1,x-r):min(H,x+r)
      selY <- max(1,y-r):min(W,y+r)
      if(min(diff(selX))<0) selX <- NULL
      if(min(diff(selY))<0) selY <- NULL
      if(length(dim(oriIMG))==2){
        allCNL <- getVisiumHDspotCol_value(oriIMG[selX,selY],vName)
      }else{
        allCNL <- list()
        for(k in 1:length(RGB)){
          allCNL <- c(allCNL,setNames(getVisiumHDspotCol_value(oriIMG[selX,selY,k],vName),
                                      paste(RGB[k],vName,sep="_")))
        }
      }
      return(data.frame(allCNL))
    },BPPARAM = MulticoreParam(workers=core))#min(length(Dlist),parallelly::availableCores()-2)
    res <- do.call(rbind,res)
    rownames(res) <- img@boundaries$centroids@cells
    return(res)
  })
  names(cInfo) <- NULL
  return(do.call(rbind,cInfo))
}
