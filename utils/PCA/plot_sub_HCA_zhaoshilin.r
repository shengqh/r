###以下为常用命令
#heatmap_3(data, Colv=NA,Rowv=NA,scale="none",legend=2)
#heatmap(data,scale="none",col=col,Colv=NA,distfun = function(x) as.dist(1 - cor(t(x))), hclustfun = function(x) hclust(x, method="av"))
#
#
#
##对称颜色可以用以下函数生成
#col=colorRampPalette(c('green','black','red'))(1024)
#col = heat.colors(1024)
#等等
#
#
#对数值不到(-1,1)的数据，想以(-1,1)为基础产生颜色，则可以用
#col=colorRampPalette(c('green','black','red'))(1024)
#col=col_part(c(-1,1),data_part,col)
#否则，以自身为基础产生颜色，用
#col=col_process(data)
#
#######################




library(Heatplus)

sub_HCA_result<-function(data_for_HCA,result_HCA,cut.off.cor,col_total,las=1,legend.plot=T,sort_col=TRUE,file="sub_HCA_result") {
	#根据cor分开聚类结果
	dend_cut<- cut(result_HCA$Rowv, h=1-cut.off.cor)
	length(dend_cut$lower)->number
	#输出到文件
	pdf(paste(file,".pdf",sep=""),width=18, height=6)
	par(mfrow=c(1,3),mar=c(6, 3.1, 3.1, 4.1),cex=1.1); 
	FN<-file(paste(file,".txt",sep=""),"w")
	for(i in 1:length(dend_cut$lower)){
		cat("Cluster_",i,"\n", file=FN,sep="")
		cat(labels(dend_cut$lower[[i]]), file=FN,sep="\n")
		result<-data_for_HCA[order.dendrogram(dend_cut$lower[[i]]),,drop=F]
		if (sort_col) {result<-result[,order.dendrogram(result_HCA$Colv),drop=F]}
		plot.submatrix(result,paste("Cluster: ",i,"\nTotal: ",nrow(result),sep=""), legend.plot=legend.plot,las=las)
		plot(dend_cut$lower[[i]],hori=T, edge.root=T)
		col=col_part(data_for_HCA,result,col_total)
		image(t(result),col=col,xaxt="n",yaxt="n")
		axis(1, at=seq(0,1,1/(ncol(result)-1)),labels=colnames(result), las=las)
	}
	close(FN)
	dev.off()
}



#调用的画图函数
#self-made function of plot submatrix
plot.submatrix<-function(data,name,legend.plot=T,las=1) {	
	plot(c(1,(length(data[1,]))),range(data),type="n",xaxt="n",xlab="Fractions index",ylab="")# set up the plot
	axis(1, at=1:ncol(data),labels=colnames(data), las=las) 
	colors <- rainbow(length(data[,1])) # add lines
	colors <-c("black",rainbow(length(data[,1])-1)) 	
	for (i in 1:length(data[,1])) {#colors<-c("red","black")
		lines(c(1:(length(data[1,]))),data[i,],type="l", lwd=2,col=colors[i])
	}	
	title(name)	# add a title and legend
	if (legend.plot) {legend(1,max(data), c(row.names(data)),cex=0.6, lty="solid",lwd=1,col=colors)}
}


#不对称颜色的处理
col_process<-function(data,middle_value=0,col_small='green',col_middle='black',col_large='red',col_number=1024) {
	col<-colorRampPalette(c(col_small,col_middle,col_large))(col_number)
	if (abs(max(data,na.rm=T))>=abs(min(data,na.rm=T))) {
		cut.off<-round(quantile(1:col_number,probs=1-(abs(max(data,na.rm=T))+abs(min(data,na.rm=T)))/(2*abs(max(data,na.rm=T)))))
		col<-col[cut.off:length(col)]
	} else {
		cut.off<-round(quantile(1:col_number,probs=(abs(max(data,na.rm=T))+abs(min(data,na.rm=T)))/(2*abs(min(data,na.rm=T)))))
		col<-col[1:cut.off]
	}
	return(col)
}

#在总体中取出部分颜色
col_part<-function(data_all,data_part,col) {
	min_all<-min(data_all,na.rm=T)
	max_all<-max(data_all,na.rm=T)
	min_part<-min(data_part,na.rm=T)
	max_part<-max(data_part,na.rm=T)
	cut_off_low<-round(quantile(1:length(col),(min_part-min_all)/(max_all-min_all)))
	cut_off_high<-round(quantile(1:length(col),(max_part-min_all)/(max_all-min_all)))
	col=col[cut_off_low:cut_off_high]
	return(col)
}




###heatmap_3
heatmap_3<-function (x, Rowv = NULL, Colv = NULL, distfun = function(x) as.dist(1 - cor(t(x))), hclustfun = hclust, 
    method.hca="average", add.expr, scale = c("row", "column", "none"), margin=c(5,0,0,5), cexRow = 0.2 + 
        1/log10(nr), cexCol = 0.2 + 1/log10(nc), na.rm = TRUE, 
    do.dendro = c(TRUE, TRUE), legend = 0, legfrac = 8, col = heat.colors(1024), 
    trim, ...) 
{
    scale <- if (missing(scale)) 
        "none" else match.arg(scale)
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1) 
        stop("`x' must have at least 2 rows and 2 columns")
    r.cex <- cexRow
    c.cex <- cexCol
    if (missing(Rowv)) 
        Rowv <- rowMeans(x, na.rm = na.rm)
    if (missing(Colv)) 
        Colv <- colMeans(x, na.rm = na.rm)
    if (!identical(Rowv, NA)) {
        if (!inherits(Rowv, "dendrogram")) {
            hcr <- hclustfun(distfun(x),method=method.hca)
            ddr <- as.dendrogram(hcr)
            ddr <- reorder(ddr, Rowv)
        }
        else ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
    } else {
        rowInd = 1:nr
        do.dendro[1] = FALSE
    }
    if (!identical(Colv, NA)) {
        if (!inherits(Colv, "dendrogram")) {
            hcc <- hclustfun(distfun(t(x)),method=method.hca)
            ddc <- as.dendrogram(hcc)
            ddc <- reorder(ddc, Colv)
        }
        else ddc <- Colv
        colInd <- order.dendrogram(ddc)
    } else {
        colInd = 1:nc
        do.dendro[2] = FALSE
    }
    x <- x[rowInd, colInd]
    if (scale == "row") {
        x <- sweep(x, 1, rowMeans(x, na.rm = na.rm))
        sd <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sd, "/")
    }
    else if (scale == "column") {
        x <- sweep(x, 2, colMeans(x, na.rm = na.rm))
        sd <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sd, "/")
    }
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    if (!missing(trim)) {
        trim = min(trim[1], 1 - trim[1])
        lo = quantile(x, trim, na.rm = na.rm)
        hi = quantile(x, 1 - trim, na.rm = na.rm)
        x[x < lo] = lo
        x[x > hi] = hi
    }
    do.xaxis = !is.null(colnames(x))
    do.yaxis = !is.null(rownames(x))
    margin[1] = if (do.xaxis) 
        margin[1]
        else 2
    margin[2] = if (do.dendro[1]) 
        margin[2]
    else 2
    margin[3] = if (do.dendro[2]) 
        margin[3]
    else 2
    margin[4] = if (do.yaxis) 
        margin[4]
    else 2
    if (do.dendro[1] & do.dendro[2]) {
        ll = matrix(c(0, 3, 2, 1), 2, 2, byrow = TRUE)
        ll.width = c(1, 4)
        ll.height = c(1, 4)
    }
    else if (do.dendro[1]) {
        ll = matrix(c(2, 1), 1, 2, byrow = TRUE)
        ll.width = c(1, 4)
        ll.height = 4
    }
    else if (do.dendro[2]) {
        ll = matrix(c(2, 1), 2, 1, byrow = FALSE)
        ll.width = 4
        ll.height = c(1, 4)
    }
    else {
        ll = matrix(1, 1, 1)
        ll.width = 1
        ll.height = 1
    }
    if (legend %in% 1:4) {
        plotnum = max(ll) + 1
        nc = ncol(ll)
        nr = nrow(ll)
        if (legend == 1) {
            ll = rbind(ll, if (nc == 1) 
                plotnum
            else c(0, plotnum))
            ll.height = c(ll.height, sum(ll.height)/(legfrac - 
                1))
            leg.hor = TRUE
        }
        else if (legend == 2) {
            ll = cbind(if (nr == 1) 
                plotnum
            else c(0, plotnum), ll)
            ll.width = c(sum(ll.width)/(legfrac - 1), ll.width)
            leg.hor = FALSE
        }
        else if (legend == 3) {
            ll = rbind(if (nc == 1) 
                plotnum
            else c(0, plotnum), ll)
            ll.height = c(sum(ll.height)/(legfrac - 1), ll.height)
            leg.hor = TRUE
        }
        else if (legend == 4) {
            ll = cbind(ll, if (nr == 1) 
                plotnum
            else c(0, plotnum))
            ll.width = c(ll.width, sum(ll.width)/(legfrac - 1))
            leg.hor = FALSE
        }
    }
    layout(ll, width = ll.width, height = ll.height, respect = TRUE)
    par(mar = margin)
    image(1:ncol(x), 1:nrow(x), t(x), axes = FALSE, xlim = c(0.5, 
        ncol(x) + 0.5), ylim = c(0.5, nrow(x) + 0.5), xlab = "", 
        ylab = "", col = col, ...)
    if (do.xaxis) {
        axis(1, 1:ncol(x), las = 2, line = -0.5, tick = 0, labels = colnames(x), 
            cex.axis = c.cex)
    }
    if (do.yaxis) {
        axis(4, 1:nrow(x), las = 2, line = -0.5, tick = 0, labels = rownames(x), 
            cex.axis = r.cex)
    }
    if (!missing(add.expr)) 
        eval(substitute(add.expr))
    if (do.dendro[1]) {
        mm = margin
        mm[2] = 3
        mm[4] = 0
        par(mar = mm)
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    if (do.dendro[2]) {
        mm = margin
        mm[1] = 0
        mm[3] = 3
        par(mar = mm)
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    if (legend %in% 1:4) {
        dummy.x <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), 
            length = length(col))
        dummy.z <- matrix(dummy.x, ncol = 1)
        if (leg.hor) {
            par(mar = c(2, margin[2], 2, margin[4]))
            image(x = dummy.x, y = 1, z = dummy.z, yaxt = "n", 
                col = col)
        }
        else {
            par(mar = c(margin[1], 2, margin[3], 2))
            image(x = 1, y = dummy.x, z = t(dummy.z), xaxt = "n", 
                col = col)
        }
    }
    invisible(list(rowInd = rowInd, colInd = colInd, Rowv = if (do.dendro[1]) ddr,
     Colv = if (do.dendro[2]) ddc))
}