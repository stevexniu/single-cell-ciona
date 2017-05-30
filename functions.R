pcScree <- function (object, pc.genes = NULL, num.pcs = 15, cumulative = FALSE)
{ 
  data.use = object@scale.data
  if (is.null(pc.genes)){
    pc.genes = object@var.genes
  } else{
    pc.genes = unique(pc.genes[pc.genes %in% rownames(data.use)])
  }
  pc.genes.var = apply(data.use[pc.genes, ], 1, var)
  pc.data = data.use[pc.genes[pc.genes.var > 0], ]
  pca.obj = prcomp(pc.data)
  pve = 100*pca.obj$sdev[1:num.pcs]^2/sum(pca.obj$sdev^2)
  if (cumulative){
    plot(cumsum(pve), type="o", ylab="Cumulative PVE", xlab="Principal Component", col="red")
  } else {
    plot(pve,type = "s", ylab = "PVE",
         xlab="Principal Component", col="blue")
  }
}

RegressOut<-function(object,latent.vars,genes.regress=NULL,do.scale=TRUE,do.spline=FALSE,
                     spline.span=0.75,do.plot=FALSE,nCol=2,num.genes=10,...) 
{
  library(pbapply)
  set.ifnull=function(x,y) {
    if(is.null(x)) x=y
    return(x)}
  genes.regress=set.ifnull(genes.regress,rownames(object@scale.data))
  genes.regress=ainb(genes.regress,rownames(object@scale.data))
  latent.data=fetch.data(object,latent.vars)
  exp.data=t(object@data[genes.regress,])
  regression.mat=cbind(latent.data,exp.data)
  new.data=t(pbsapply(genes.regress, function(x) {
    regression.mat.2=latent.data
    regression.mat.2[,"GENE"] = regression.mat[,x];
    fmla=as.formula(paste("GENE ", " ~ ", paste(latent.vars,collapse="+"),sep=""));
    return(lm(fmla,data = regression.mat.2)$residuals)
  }))
  
  if(do.spline==TRUE){
    new.data=t(pbsapply(genes.regress, function(x) {
      regression.mat.2=latent.data
      regression.mat.2[,"GENE"] = regression.mat[,x];
      fmla=as.formula(paste("GENE ", " ~ ", paste(latent.vars,collapse="+"),sep=""));
      return(loess(fmla,span = spline.span, data = regression.mat.2)$residuals)
    }))
  }
  
  if (do.plot==TRUE){
    PCs=unique(as.numeric(unlist(strsplit(gsub("[^0-9]", "", unlist(latent.vars)), ""))))
    genes.plot=pcTopGenes(object,PCs,num.genes)
    nCol=nCol
    par(mfrow=c(nCol,length(genes.plot)/nCol))
    for(gene in genes.plot){
      fm=as.formula(paste("`",gene,"`","~", paste(latent.vars,collapse="+"),sep = ""))
      plot(regression.mat[,latent.vars[1]],regression.mat[,gene],xlab = latent.vars[1],ylab=as.character(gene))
      if(do.spline==TRUE){
        loess(fm,span = spline.span, data = regression.mat)
      }else{abline(lm(fm,data = regression.mat),col="red",lwd=2)}
    }
    par(mfrow=c(1,1))
  }
  
  if (do.scale==TRUE) {
    new.data=t(scale(t(new.data)))
  }
  new.data[is.na(new.data)]=0
  object@scale.data=new.data
  return(object)
  
}

feature.plot.pseudo <- function (object, features.plot, dim.1 = 1, dim.2 = 2, cells.use = NULL, 
                                 pt.size = 1, cols.use = heat.colors(10), pch.use = 16, reduction.use = "tsne", 
                                 use.imputed = FALSE, nCol = NULL,name.x=NULL,name.y=NULL,yaxt=NULL,plot.legend=FALSE,...) 
{
  set.ifnull=function(x,y) {
    if(is.null(x)) x=y
    return(x)
  }
  
  translate.dim.code <- function (reduction.use) 
  {
    return.code = "PC"
    if (reduction.use == "ica") 
      return.code = "IC"
    if (reduction.use == "tsne") 
      return.code = "tSNE_"
    if (reduction.use == "mds") 
      return.code = "MDS"
    return(return.code)
  }
  
  cells.use = set.ifnull(cells.use, colnames(object@data))
  dim.code = "PC"
  if (is.null(nCol)) {
    nCol = 2
    if (length(features.plot) > 6) 
      nCol = 3
    if (length(features.plot) > 9) 
      nCol = 4
  }
  num.row = floor(length(features.plot)/nCol - 1e-05) + 1
  par(mfrow = c(num.row, nCol))
  dim.code = translate.dim.code(reduction.use)
  dim.codes = paste(dim.code, c(dim.1, dim.2), sep = "")
  data.plot = fetch.data(object, dim.codes)
  x1 = paste(dim.code, dim.1, sep = "")
  x2 = paste(dim.code, dim.2, sep = "")
  data.plot$x = data.plot[, x1]
  data.plot$y = data.plot[, x2]
  data.plot$pt.size = pt.size
  data.use = data.frame(t(fetch.data(object, features.plot, 
                                     cells.use = cells.use, use.imputed = use.imputed)))
  for (i in features.plot) {
    data.gene = na.omit(data.frame(data.use[i, ]))
    data.cut = as.numeric(as.factor(cut(as.numeric(data.gene), 
                                        breaks = length(cols.use))))
    data.col = rev(cols.use)[data.cut]
    name.x=set.ifnull(name.x,"tSNE_1")
    name.y=set.ifnull(name.y,"tSNE_2")
    plot(data.plot$x, data.plot$y, col = data.col, cex = pt.size, 
         pch = pch.use, main = i, xlab = name.x, ylab = name.y,yaxt=yaxt,...)
    if(plot.legend) legend("topright",legend = c("High",rep("",8),"Low"),border = F,fill = heat.colors(10),y.intersp = 0.2,bty="n",cex = 0.8)
  }
  rp()
}

jackStrawPlot.new<-function (object, PCs = 1:5, nCol = 3, score.thresh = 1e-05, 
                             plot.x.lim = 0.1, plot.y.lim = 0.3) 
{
  pAll = object@jackStraw.empP
  pAll <- pAll[, PCs, drop = F]
  pAll$Contig <- rownames(pAll)
  pAll.l <- melt(pAll, id.vars = "Contig")
  colnames(pAll.l) <- c("Contig", "PC", "Value")
  qq.df <- NULL
  score.df <- NULL
  for (i in PCs) {
    q <- qqplot(pAll[, i], runif(1000), plot.it = FALSE)
    pc.score = prop.test(c(length(which(pAll[, i] <= score.thresh)), 
                           floor(nrow(pAll) * score.thresh)), c(nrow(pAll), 
                                                                nrow(pAll)))$p.val
    if (length(which(pAll[, i] <= score.thresh)) == 0) 
      pc.score = 1
    if (is.null(score.df)) 
      score.df <- data.frame(PC = paste("PC", i, sep = ""), 
                             Score = pc.score)
    else score.df <- rbind(score.df, data.frame(PC = paste("PC", 
                                                           i, sep = ""), Score = pc.score))
    if (is.null(qq.df)) 
      qq.df <- data.frame(x = q$x, y = q$y, PC = paste("PC", 
                                                       i, sep = ""))
    else qq.df <- rbind(qq.df, data.frame(x = q$x, y = q$y, 
                                          PC = paste("PC", i, sep = "")))
  }
  pAll.l$PC.Score <- paste(score.df$PC, sprintf("%1.3g", score.df$Score))
  gp <- ggplot(pAll.l, aes(sample = Value)) + stat_qq(distribution = qunif) + 
    facet_wrap("PC.Score", ncol = nCol) + labs(x = "Theoretical [runif(1000)]", 
                                               y = "Empirical") + xlim(0, plot.y.lim) + ylim(0, plot.x.lim) + 
    coord_flip() + geom_abline(intercept = 0, slope = 1, 
                               linetype = "dashed", na.rm = T) + theme_bw()
  return(gp)
}

tsne.pseudo <- function (object, do.label = FALSE, do.col.ret=FALSE,label.pt.size = 1, label.cex.text = 1,
                         label.cols.use = NULL, name.x=NULL,name.y=NULL,font.use=NULL, ...) 
{
  set.ifnull <- function(x,y) {
    if(is.null(x)) x=y
    return(x)
  }
  
  .local <- function (object, do.label = FALSE, do.col.ret=FALSE,label.pt.size = 1, 
                      label.cex.text = 1, label.cols.use = NULL, cells.use = NULL,
                      ...) 
  {
    cells.use = set.ifnull(cells.use, object@cell.names)
    cells.use = ainb(cells.use, object@cell.names)
    name.x = set.ifnull(name.x,"tSNE_1")
    name.y = set.ifnull(name.y,"tSNE_2")
    if (do.label == TRUE) {
      label.cols.use = set.ifnull(label.cols.use, rainbow(length(levels(object@ident[cells.use]))))
      set.seed(1)
      label.cols.use = sample(label.cols.use)
      if (do.col.ret == TRUE){
        print(label.cols.use) 
      }
      if (length(unique(object@tsne.rot[,2]))==1){
        plot(object@tsne.rot[cells.use, 1], object@tsne.rot[cells.use,2],col = label.cols.use[as.integer(object@ident[cells.use])], 
             pch = 16, xlab = name.x, ylab = name.y, cex = label.pt.size,yaxt='n',cex.lab=1.2,font.lab=2,...)
        k.centers = t(sapply(levels(object@ident), function(x) apply(object@tsne.rot[which.cells(object, 
                                                                                                 x), ], 2, median)))
        points(k.centers[, 1], k.centers[, 2],cex = 1.3, 
               col = "white", pch = 16)
        text(k.centers[, 1], k.centers[, 2],levels(object@ident), 
             cex = label.cex.text,font=font.use)
      } else {
        plot(object@tsne.rot[cells.use, 1], object@tsne.rot[cells.use, 
                                                            2], col = label.cols.use[as.integer(object@ident[cells.use])], 
             pch = 16, xlab = name.x, ylab = name.y, cex = label.pt.size,cex.lab=1.2, font.lab=2,...)
        k.centers = t(sapply(levels(object@ident), function(x) apply(object@tsne.rot[which.cells(object, 
                                                                                                 x), ], 2, median)))
        points(k.centers[, 1], k.centers[, 2], cex = 1.3, 
               col = "white", pch = 16)
        text(k.centers[, 1], k.centers[, 2], levels(object@ident), 
             cex = label.cex.text,font=font.use)
      }
    }
    else {
      if (length(unique(object@tsne.rot[,2]))==1){
        plot(object@tsne.rot[cells.use, 1], object@tsne.rot[cells.use,2],col = label.cols.use[as.integer(object@ident[cells.use])], 
             pch = 16, xlab = name.x, ylab = name.y, cex = label.pt.size,yaxt='n',cex.lab=1.2,font.lab=2,...)
      } else {
        plot(object@tsne.rot[cells.use, 1], object@tsne.rot[cells.use, 
                                                            2], col = label.cols.use[as.integer(object@ident[cells.use])], 
             pch = 16, xlab = name.x, ylab = name.y, cex = label.pt.size,cex.lab=1.2, font.lab=2,...)
      }
    }
  }
  .local(object, do.label, do.col.ret, label.pt.size, label.cex.text, label.cols.use, 
         ...) 
}

genes.plot.pseudo <- function (object,genes.use, cell.ids = NULL, pseudo="pseudo.time",pch.use = 16, cex.use = 1.5, do.induct=TRUE,
                               do.scale=FALSE,use.imputed = FALSE, method = "logit", spline.span = 0.75,name.x=NULL,do.line=FALSE,col.use=NULL, bin = 25, nState = 2,
                               do.smooth=FALSE,knn=10,inset=c(0,0),ylim=NULL,conf.int=FALSE,do.label=F,do.ret.mat=F,lwd=1,do.ret.pt=F,use.loess=F,...) 
{
  set.ifnull=function(x,y) {
    if(is.null(x)) x=y
    return(x)
  }
  
  nn.count.add <- function(mat, k, dmat) {
    dmat[dmat == 0] <- Inf
    knn <- t(apply(dmat, 1, function(x) rank(x) <= k))
    mat <- mat + mat %*% t(knn)
    return(mat)
  }
  
  cell.ids = set.ifnull(cell.ids, object@cell.names)
  ident.use = as.factor(object@ident[cell.ids])
  xlab=paste(name.x,"Pseudotime")
  ret.mat=c()
  turn.mat=c()
  turn.pt=list(on=c(),off=c())
  if(do.induct & method == "hmm") {
    library(RHmm)
    object.pseudo = pseudoBin(object, bin.size = bin)}
    
  if (do.smooth){
    pt.ncc=fetch.data(object,pseudo)
    dmat = as.matrix(dist(pt.ncc))
    count.mat = object@scale.data
    object@scale.data = nn.count.add(count.mat, knn, dmat)
    data.use = data.frame(t(fetch.data(object, c(pseudo, genes.use), 
                                       cells.use = cell.ids, use.imputed = use.imputed,use.scaled = TRUE)))
  }
  
  if(do.scale){
    data.use = data.frame(t(fetch.data(object, c(pseudo, genes.use), 
                                       cells.use = cell.ids, use.imputed = use.imputed,use.scaled = TRUE)))
  }else{
    data.use = data.frame(t(fetch.data(object, c(pseudo, genes.use), 
                                       cells.use = cell.ids, use.imputed = use.imputed,use.scaled = FALSE)))
  }
  
  pseudo.order=order(data.use[pseudo,cell.ids])
  pseudo.data=unlist(data.use[pseudo,cell.ids][pseudo.order])
  object@scale.data=object@scale.data[,pseudo.order]
  ylim=set.ifnull(ylim,range(data.use[genes.use,cell.ids]))
  
  plot(0,0,xlim = range(pseudo.data),ylim = ylim,type = "l",xlab=xlab,ylab="log(FPKM)",
       cex = cex.use, pch = pch.use, font.lab=2,cex.lab=1.2,bty="l",...)
  col.use=set.ifnull(col.use,rainbow(length(genes.use)))
  if(do.line){
    for (i in 1:length(genes.use)){
      g2 = as.numeric(data.use[genes.use[i], cell.ids][pseudo.order])
      points(pseudo.data, g2,col=col.use[i],pch=20)
      lines(pseudo.data, g2,col=col.use[i],lwd=lwd)
    }
  }else{
    for (i in 1:length(genes.use)){
      g2 = as.numeric(data.use[genes.use[i], cell.ids][pseudo.order])
      loess.fit = loess(g2 ~ pseudo.data, span = spline.span)
      lines(pseudo.data, loess.fit$fitted,type="l",col=col.use[i],lwd=lwd)
      ret.mat=rbind(ret.mat,loess.fit$fitted)
    
    if(conf.int){
      prid = predict(loess.fit,se=T)
      lines(pseudo.data[order(pseudo.data)], (prid$fit - qt(0.975,prid$df) * prid$se)[order(pseudo.data)],col=col.use[i],lty=2)
      lines(pseudo.data[order(pseudo.data)], (prid$fit + qt(0.975,prid$df) * prid$se)[order(pseudo.data)],col=col.use[i],lty=2)
    }
      
    if(do.induct & method == "logit"){
      if(do.smooth){
        gene.norm=apply(matrix(object@scale.data[genes.use[i],cell.ids]),2, FUN = function(X) (X - min(X))/diff(range(X)))
      }
      
      if(use.loess){
        logit.use <- loess.fit$fitted
        gene.norm=(logit.use - min(logit.use))/diff(range(logit.use))}
      
      if(!use.loess & !do.smooth){
        gene.norm=apply(object@data[genes.use[i],cell.ids],1, FUN = function(X) (X - min(X))/diff(range(X)))
      }
      
      gene.norm[which(gene.norm>0.5)]=1
      gene.norm[which(gene.norm<0.5)]=0
      model <- glm(gene.norm[pseudo.order]~pseudo.data,family=binomial(link='logit'))
      model.fit=abs(fitted(model)-0.5)
      turn=min(model.fit)
      turn.point=pseudo.data[which(model.fit==turn)]
      turn.mat=c(turn.mat,turn.point)
      abline(v=turn.point,lwd=2,col=col.use[i],lty="longdash")
    }
      if(method == "hmm"){
        if(use.loess){hmm.use <- loess.fit$fitted} else {hmm.use <- unlist(object@data[genes.use[i],])}
        if(max(hmm.use) <= 2){
          turn=1
          turn.point=pseudo.data[turn]
          names(turn.point)=genes.use[i]
          abline(v=turn.point,lwd=2,col=col.use[i],lty="longdash")
          turn.pt$off=c(turn.pt$off,turn.point)
        } else {
          hmm.fit=HMMFit(hmm.use,nStates = nState)
          hmm.vit=viterbi(hmm.fit,hmm.use)
          if(length(unique(hmm.vit$states)) == 1){
            turn=1
            turn.point=pseudo.data[turn]
            names(turn.point)=genes.use[i]
            abline(v=turn.point,lwd=2,col=col.use[i],lty="longdash")
            if(loess.fit$fitted[turn] >=1){
              turn.pt$on=c(turn.pt$on,turn.point)
            } else{
              turn.pt$off=c(turn.pt$off,turn.point)
            }
          } else {
            counter = 1
            turn = NULL
            for (j in 1:(length(hmm.use)-1)){
              if(hmm.vit$states[j] - hmm.vit$states[j+1] != 0){
                turn.point=pseudo.data[j]
                names(turn.point)=genes.use[i]
                abline(v=turn.point,lwd=2,col=col.use[i],lty="longdash")
                if(hmm.use[j]-hmm.use[j+1]>0){
                  turn.pt$off = c(turn.pt$off, turn.point)
                } else {turn.pt$on=c(turn.pt$on, turn.point)}
              }
            }
          }
        }
      }
  }
  }
  if(do.label){
    legend("topright",legend = genes.use, 
           lty=1, col=col.use, cex=.75,xpd=T,inset = inset)
  }
  if(do.ret.mat){ 
    rownames(ret.mat)=genes.use
    colnames(ret.mat)=cell.ids[pseudo.order]
    return(ret.mat)}
  if(do.ret.pt){
    if(method=="hmm") {return(turn.pt)} else {
      names(turn.mat)=genes.use
      return(turn.mat)
    }
  }
}

genePlot.pseudo <- function (object,gene,pseudo="pseudo.time", cell.ids = NULL, col.use = NULL, 
                             pch.use = 16, cex.use = 1.5, use.imputed = FALSE, do.ident = FALSE, 
                             do.spline = FALSE,do.line=FALSE,spline.span = 0.75,name.x=NULL, do.logit=FALSE,
                             do.smooth=FALSE,knn=10,use.scale=FALSE,do.return=FALSE,conf.int=FALSE,...) 
{
  set.ifnull=function(x,y) {
    if(is.null(x)) x=y
    return(x)
  }
  nn.count.add <- function(mat, k, dmat) {
    dmat[dmat == 0] <- Inf
    knn <- t(apply(dmat, 1, function(x) rank(x) <= k))
    mat <- mat + mat %*% t(knn)
    return(mat)
  }
  
  if (do.smooth){
    pt.ncc=fetch.data(object,pseudo)
    dmat = as.matrix(dist(pt.ncc))
    count.mat = object@scale.data
    object@scale.data = nn.count.add(count.mat, knn, dmat)
    data.use = data.frame(t(fetch.data(object, c(pseudo, gene), 
                                       cells.use = cell.ids, use.imputed = use.imputed,use.scaled = TRUE)))
  }else{
    data.use = data.frame(t(fetch.data(object, c(pseudo, gene), 
                                       cells.use = cell.ids, use.imputed = use.imputed,use.scaled = use.scale)))
  }
  
  cell.ids = set.ifnull(cell.ids, object@cell.names)
  name.x = set.ifnull(name.x, pseudo)
  pseudo.order=order(data.use[pseudo,cell.ids])
  g1 = as.numeric(data.use[pseudo, cell.ids])
  g2 = as.numeric(data.use[gene, cell.ids])
  ident.use = as.factor(object@ident[cell.ids])
  pseudo.data=as.matrix(fetch.data(object,pseudo))
  
  if(do.logit){
    if(do.smooth){
      gene.norm=apply(matrix(object@scale.data[gene,colnames(object@scale.data)]),2, FUN = function(X) (X - min(X))/diff(range(X)))
    }
    
    if(do.spline){
      loess.fit = loess(g2 ~ g1, span = spline.span)
      gene.norm = (loess.fit$fitted - min(loess.fit$fitted))/diff(range(loess.fit$fitted))
    }
    
    if(!do.smooth & !do.spline){
      gene.norm=apply(object@data[gene,colnames(object@data)],1, FUN = function(X) (X - min(X))/diff(range(X)))
    }
    
    gene.norm[which(gene.norm>0.5)]=1
    gene.norm[which(gene.norm<0.5)]=0
    model <- glm(gene.norm[pseudo.order]~g1,family=binomial(link='logit'))
    model.fit=abs(fitted(model)-0.5)
    turn=min(model.fit)
    turn.point=pseudo.data[which(model.fit==turn),1]
  }
  
  if (length(col.use) > 1) {
    col.use = col.use[as.numeric(ident.use)]
  }
  else {
    col.use = set.ifnull(col.use, as.numeric(ident.use))
  }
  
  g1 = as.numeric(data.use[pseudo, cell.ids])
  g2 = unlist(object@data[gene, cell.ids])
  gene.cor = round(cor(g1, g2), 2)
  plot(g1, g2, xlab = name.x, ylab = gene, col = col.use, cex = cex.use, 
       main = "", pch = pch.use, font.lab=2,cex.lab=1.2, ...)
  
  if (do.logit){
    abline(v=turn.point,lwd=4,col="purple",lty="longdash")
  }
  
  if (do.spline) {
    loess.fit = loess(g2 ~ g1, span = spline.span)
    lines(g1[order(pseudo.data)], loess.fit$fitted[order(pseudo.data)],col = "black",lwd=4)
    if(conf.int){
      prid = predict(loess.fit,se=T)
      lines(g1[order(g1)], (prid$fit - qt(0.975,prid$df) * prid$se)[order(g1)],col = "black",lwd=1.2,lty=2)
      lines(g1[order(g1)], (prid$fit + qt(0.975,prid$df) * prid$se)[order(g1)],col = "black",lwd=1.2,lty=2)
    }
  }
  
  if(do.line){
    lines(g1[order(pseudo.data)], g2[order(pseudo.data)],col = "purple",lwd=3)
  }
  
  if (do.ident) {
    return(identify(g1, g2, labels = cell.ids))
  }
  
  if (do.return){
    if (do.spline){
      return(data.frame(g1,loess.fit$fitted))}else{
        return(data.frame(g1,g2))
      }
  }
}

vlnPlot.FPKM <- function (object, features.plot, name=NULL, nCol = NULL, ylab.max = 12, 
                          do.ret = FALSE, do.sort = FALSE, size.x.use = 16, size.y.use = 20, 
                          size.title.use = 20, use.imputed = FALSE, adjust.use = 1, do.scale=FALSE,
                          size.use = 1, cols.use = NULL, group.by = NULL,name.y=NULL,ratio.plot=1,name.x=NULL,x.axis=TRUE) 
{ 
  set.ifnull=function(x,y) {
    if(is.null(x)) x=y
    return(x)
  }
  
  plt.Vln.FPKM<-function (gene, data, cell.ident, ylab.max = 12, do.ret = FALSE, 
                          do.sort = FALSE, size.x.use = 16, size.y.use = 16, size.title.use = 20, 
                          adjust.use = 1, size.use = 1, cols.use = NULL,x.axis=TRUE) 
  {
    data$gene = as.character(rownames(data))
    data.use = data.frame(data[gene, ])
    if (length(gene) == 1) {
      data.melt = data.frame(rep(gene, length(cell.ident)))
      colnames(data.melt)[1] = "gene"
      data.melt$value = as.numeric(data[1, 1:length(cell.ident)])
      data.melt$id = names(data)[1:length(cell.ident)]
    }
    if (length(gene) > 1) 
      data.melt = melt(data.use, id = "gene")
    data.melt$ident = cell.ident
    noise <- rnorm(length(data.melt$value))/1e+05
    data.melt$value = as.numeric(as.character(data.melt$value)) + 
      noise
    if (do.sort) {
      data.melt$ident = factor(data.melt$ident, levels = names(rev(sort(tapply(data.melt$value, 
                                                                               data.melt$ident, mean)))))
    }
    name.y=set.ifnull(name.y,"Expression level (log FPKM)\n")
    p = ggplot(data.melt, aes(factor(ident), value))
    p2 = p + geom_violin(scale = "width", adjust = adjust.use, 
                         trim = TRUE, aes(fill = factor(ident))) + ylab(name.y)
    if (!is.null(cols.use)) {
      p2 = p2 + scale_fill_manual(values = cols.use)
    }
    name=set.ifnull(name,gene)
    name.x=set.ifnull(name.x,"\nCell Type")
    p3 = p2 + theme(legend.position = "top") + guides(fill = guide_legend(title = NULL))  
      #geom_jitter(height = 0, size = size.use) + xlab(name.x)
    if (!x.axis) {
      p4 = p3 + theme(axis.title.x = element_text(face = "bold", 
                                                  colour = "#990000", size = size.x.use), axis.text.x = element_text(angle = 90, 
                                                                                                                     vjust = 0.5, size = 18)) + theme_bw() + theme(
                                                                                                                       plot.background = element_blank(),
                                                                                                                       panel.grid.major = element_blank(),
                                                                                                                       panel.grid.minor = element_blank()) + scale_x_continuous(labels = "")
    } else {
      p4 = p3 + theme(axis.title.x = element_text(face = "bold", 
                                                  colour = "#990000", size = size.x.use), axis.text.x = element_text(angle = 90, 
                                                                                                                     vjust = 0.5, size = 18)) + theme_bw() + theme(
                                                                                                                       plot.background = element_blank(),
                                                                                                                       panel.grid.major = element_blank(),
                                                                                                                       panel.grid.minor = element_blank())}
    gene=name
    
    p5 = (p4 + theme(axis.title.y = element_text(face = "bold", 
                                                 colour = "#990000", size = size.y.use,vjust=1.5), axis.text.y = element_text(angle = 90, 
                                                                                                                              vjust = 0.5, size = 18)) + ggtitle(gene) + theme(plot.title = element_text(size = size.title.use, 
                                                                                                                                                                                                         face = "bold")))
    p5 = p5 + theme(axis.text=element_text(size=18,vjust = 2),axis.title.x=element_text(size=18,face = "bold"))+theme(legend.position = "none") +
      coord_fixed(ratio=ratio.plot)
    if (do.ret == TRUE) {
      return(p5)
    }
    else {
      print(p5)
    }
  }
  
  pltVln<-function (gene, data, cell.ident, ylab.max = 12, do.ret = FALSE, 
                    do.sort = FALSE, size.x.use = 16, size.y.use = 16, size.title.use = 20, 
                    adjust.use = 1, size.use = 1, cols.use = NULL) 
  {
    data$gene = as.character(rownames(data))
    data.use = data.frame(data[gene, ])
    if (length(gene) == 1) {
      data.melt = data.frame(rep(gene, length(cell.ident)))
      colnames(data.melt)[1] = "gene"
      data.melt$value = as.numeric(data[1, 1:length(cell.ident)])
      data.melt$id = names(data)[1:length(cell.ident)]
    }
    if (length(gene) > 1) 
      data.melt = melt(data.use, id = "gene")
    data.melt$ident = cell.ident
    noise <- rnorm(length(data.melt$value))/1e+05
    data.melt$value = as.numeric(as.character(data.melt$value)) + 
      noise
    if (do.sort) {
      data.melt$ident = factor(data.melt$ident, levels = names(rev(sort(tapply(data.melt$value, 
                                                                               data.melt$ident, mean)))))
    }
    p = ggplot(data.melt, aes(factor(ident), value))
    p2 = p + geom_violin(scale = "width", adjust = adjust.use, 
                         trim = F, show.legend = F, aes(fill = factor(ident))) + geom_boxplot(width=0.15)+ ylab("Expression level (logFPKM)")
    if (!is.null(cols.use)) {
      p2 = p2 + scale_fill_manual(values = cols.use)
    }
    p3 = p2 + ylim(2000,7000)
      #geom_jitter(height = 0, size = size.use) + xlab(NULL)
    p4 = p3 + theme(axis.title.x = element_text(face = "bold", 
                                                colour = "#990000", size = size.x.use), axis.text.x = element_text(angle = 90, 
                                                                                                                   vjust = 0.5, size = 5)) + theme_bw() + theme(
                                                                                                                     plot.background = element_blank(),
                                                                                                                     panel.grid.major = element_blank(),
                                                                                                                     panel.grid.minor = element_blank())
    p5 = (p4 + theme(axis.title.y = element_text(face = "bold", 
                                                 colour = "#990000", size = size.y.use), axis.text.y = element_text(angle = 90, 
                                                                                                                    vjust = 0.5, size = 12)) + ggtitle(gene) + theme(plot.title = element_text(size = size.title.use, 
                                                                                                                                                                                               face = "bold")))
    if (do.ret == TRUE) {
      return(p5)
    }
    else {
      print(p5)
    }
  }
  
  if (is.null(nCol)) {
    nCol = 2
    if (length(features.plot) > 6) 
      nCol = 3
    if (length(features.plot) > 9) 
      nCol = 4
  }
  if(do.scale){ 
    data.use = data.frame(t(fetch.data(object, features.plot, 
                                                   use.imputed = use.imputed,use.scaled = T)))
  }else{
  data.use = data.frame(t(fetch.data(object, features.plot, 
                                     use.imputed = use.imputed)))
  }
  ident.use = object@ident
  if (!is.null(group.by)) 
    ident.use = as.factor(fetch.data(object, group.by)[, 
                                                       1])
  pList = lapply(features.plot, function(x) pltVln(x, data.use[x, 
                                                                     ], ident.use, ylab.max, TRUE, do.sort, size.x.use, size.y.use, 
                                                         size.title.use, adjust.use, size.use, cols.use))
  if (do.ret) {
    return(pList)
  }
  else {
    multiplotList(pList, cols = nCol)
    rp()
  }
}

boxPlot.FPKM <- function (object, features.plot, name=NULL, nCol = NULL, ylab.max = 12, 
                          do.ret = FALSE, do.sort = FALSE, size.x.use = 20, size.y.use = 20, 
                          size.title.use = 20, use.imputed = FALSE, adjust.use = 1, 
                          size.use = 1, cols.use = NULL, group.by = NULL,name.y=NULL,ratio.plot=1,name.x=NULL,do.jitter=FALSE) 
{
  set.ifnull=function(x,y) {
    if(is.null(x)) x=y
    return(x)
  }
  
  plt.Box.FPKM<-function (gene, data, cell.ident, ylab.max = 12, do.ret = FALSE, 
                          do.sort = FALSE, size.x.use = 16, size.y.use = 16, size.title.use = 20, 
                          adjust.use = 1, size.use = 1, cols.use = NULL) 
  {
    data$gene = as.character(rownames(data))
    data.use = data.frame(data[gene, ])
    if (length(gene) == 1) {
      data.melt = data.frame(rep(gene, length(cell.ident)))
      colnames(data.melt)[1] = "gene"
      data.melt$value = as.numeric(data[1, 1:length(cell.ident)])
      data.melt$id = names(data)[1:length(cell.ident)]
    }
    if (length(gene) > 1) 
      data.melt = melt(data.use, id = "gene")
    data.melt$ident = cell.ident
    noise <- rnorm(length(data.melt$value))/1e+05
    data.melt$value = as.numeric(as.character(data.melt$value)) + 
      noise
    if (do.sort) {
      data.melt$ident = factor(data.melt$ident, levels = names(rev(sort(tapply(data.melt$value, 
                                                                               data.melt$ident, mean)))))
    }
    name.y=set.ifnull(name.y,"Expression level (log FPKM)\n")
    p = ggplot(data.melt, aes(factor(ident), value))
    p2 = p + geom_boxplot(aes(fill = factor(ident))) + ylab(name.y) + xlab(name.x)
    if (!is.null(cols.use)) {
      p2 = p2 + scale_fill_manual(values = cols.use)
    }
    name=set.ifnull(name,gene)
    name.x=set.ifnull(name.x,"\nCell Type")
    if(do.jitter){
      p3 = p2 + theme(legend.position = "top") + guides(fill = guide_legend(title = NULL)) + xlab(name.x) +  geom_jitter(height = 0, size = size.use)
    } else {
      p3 = p2 + theme(legend.position = "top") + guides(fill = guide_legend(title = NULL)) + xlab(name.x)
    }
    p4 = p3 + theme(axis.title.x = element_text(face = "bold", 
                                                colour = "#990000", size = size.x.use), axis.text.x = element_text(angle = 90, 
                                                                                                                   vjust = 0.5, size = 18)) + theme_bw() + theme(
                                                                                                                     plot.background = element_blank(),
                                                                                                                     panel.grid.major = element_blank(),
                                                                                                                     panel.grid.minor = element_blank())
    if(name==""){gene=NULL}
    
    p5 = (p4 + theme(axis.title.y = element_text(face = "bold", 
                                                 colour = "#990000", size = size.y.use,vjust=1.5), axis.text.y = element_text(angle = 90, 
                                                                                                                              vjust = 0.5, size = 18)) + ggtitle(gene) + theme(plot.title = element_text(size = size.title.use, 
                                                                                                                                                                                                         face = "bold")))
    p5 = p5 + theme(axis.text=element_text(size=18,vjust = 2),axis.title.x=element_text(size=18,face = "bold"))+theme(legend.position = "none") +
      coord_fixed(ratio=ratio.plot)
    if (do.ret == TRUE) {
      return(p5)
    }
    else {
      print(p5)
    }
  }
  
  
  if (is.null(nCol)) {
    nCol = 2
    if (length(features.plot) > 6) 
      nCol = 3
    if (length(features.plot) > 9) 
      nCol = 4
  }
  data.use = data.frame(t(fetch.data(object, features.plot, 
                                     use.imputed = use.imputed)))
  ident.use = object@ident
  if (!is.null(group.by)) 
    ident.use = as.factor(fetch.data(object, group.by)[, 
                                                       1])
  pList = lapply(features.plot, function(x) plt.Box.FPKM(x, data.use[x, 
                                                                     ], ident.use, ylab.max, TRUE, do.sort, size.x.use, size.y.use, 
                                                         size.title.use, adjust.use, size.use, cols.use))
  if (do.ret) {
    return(pList)
  }
  else {
    multiplotList(pList, cols = nCol)
    rp()
  }
}

pseudo.gene.cluster <- function (object, cell.ids=NULL, pseudo="pseudo.time",genes.use,do.smooth=F,do.spline=F, 
                                 spline.span = 0.75, knn=10,...) 
{  
  set.ifnull=function(x,y) {
    if(is.null(x)) x=y
    return(x)
  }
  
  nn.count.add <- function(mat, k, dmat) {
    dmat[dmat == 0] <- Inf
    knn <- t(apply(dmat, 1, function(x) rank(x) <= k))
    mat <- mat + mat %*% t(knn)
    return(mat)
  }
  
  if (do.smooth){
    pt.ncc=fetch.data(object,pseudo)
    dmat = as.matrix(dist(pt.ncc))
    count.mat = object@scale.data
    object@scale.data = nn.count.add(count.mat, knn, dmat)
  }
  
  pseudo.data=as.matrix(fetch.data(object,pseudo))
  pseudo.order=order(pseudo.data)
  pseudo.data=pseudo.data[pseudo.order]
  cell.ids = set.ifnull(cell.ids, object@cell.names)
  cluster=list(up=NULL,down=NULL)
  turn.mat=c()
  ret.mat=c()
  
  for (gene in genes.use){
    if(do.spline){
      loess.fit = loess(unlist(object@data[gene, pseudo.order]) ~ pseudo.data, span = spline.span)
      data.use=loess.fit$fitted
    } else {
      data.use=unlist(object@data[gene,])
    }
    gene.norm=(data.use - min(data.use))/diff(range(data.use))
    gene.norm[which(gene.norm>0.5)]=1
    gene.norm[which(gene.norm<0.5)]=0
    model <- glm(gene.norm[pseudo.order]~pseudo.data,family=binomial(link='logit'))
    model.fit=fitted(model)
    turn=which.min(abs(model.fit-0.5))
    turn.point=pseudo.data[turn]
    attr(gene,"names")=turn.point

    if(turn==1 | turn==length(gene.norm)){
      if(turn==1){
        if(model.fit[turn+1]>model.fit[turn]){cluster[[1]]=append(cluster[[1]],gene)}
        if(model.fit[turn+1]<model.fit[turn]){cluster[[2]]=append(cluster[[2]],gene)}
      }
      if(turn==length(gene.norm)){
        if(model.fit[turn]>model.fit[turn-1]){cluster[[1]]=append(cluster[[1]],gene)}
        if(model.fit[turn]<model.fit[turn-1]){cluster[[2]]=append(cluster[[2]],gene)}
      }
    } else {
      if (model.fit[turn]>model.fit[turn-1] || model.fit[turn+1]>model.fit[turn]){
        cluster[[1]]=append(cluster[[1]],gene)
      }else{
        cluster[[2]]=append(cluster[[2]],gene)
      }
    }
  }
  cluster1=list(up=NULL,down=NULL)
  cluster1$up=cluster[[1]][order(as.double(names(cluster[[1]])),decreasing = T)]
  cluster1$down=cluster[[2]][order(as.double(names(cluster[[2]])),decreasing = F)]
  return(cluster1)
}

scale.pseudo<-function(object,pseudo="pseudo.time",cell.ids=NULL,do.PC=FALSE,do.scale=FALSE)
{ 
  set.ifnull=function(x,y) {
    if(is.null(x)) x=y
    return(x)
  }
  
  if(do.scale){
    cell.ids=set.ifnull(cell.ids,colnames(object@data))
    pseudo.data=matrix(t(fetch.data(object,pseudo,cells.use = cell.ids)))
    pseudo.scale=apply(pseudo.data,2, FUN = function(X) (X - min(X))/diff(range(X)))
    object@data.info[,"pseudo.time"]=pseudo.scale 
  }
  
  if(do.PC){
    if(do.scale){object@pca.rot[,1]=pseudo.scale}
    else{object@pca.rot[,1]=object@data.info[,"pseudo.time"]}
  }
  return(object)
}
 
