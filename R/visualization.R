

# plot UMAP for the raw data
plot.umap = function(x, labels, main, colors, pos = "topright", pad=0.1, cex=1, pch=19, add=FALSE, legend.suffix="", cex.main=1.5, cex.legend=1.5) {
  layout = x
  if (is(x, "umap")) {
    layout = x$layout
  }
  
  xylim = range(layout)
  xylim = xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
  if (!add) {
    par(mar=c(0.2,0.7,1.2,0.7), ps=10)
    plot(xylim, xylim, type="n", axes=F, frame=F, xlab = "UMAP1")
    rect(xylim[1], xylim[1], xylim[2], xylim[2], border="#aaaaaa", lwd=0.25)
  }
  points(layout[,1], layout[,2], col=colors[as.integer(labels)],
         cex=cex, pch=pch)
  mtext(side=3, main, cex=cex.main)
  
  
  labels.u = unique(labels)
  legend.pos = pos
  legend.text = as.character(labels.u)
  if (add) {
    legend.pos = "bottomleft"
    legend.text = paste(as.character(labels.u), legend.suffix)
  }
  
  legend(legend.pos, legend=legend.text, inset=0.03,
         col=colors[as.integer(labels.u)],
         bty="n", pch=pch, cex=cex.legend)
}

# functions for plotting the heatmaps

clu_sep <- function(CZ_best){
  v <- c()
  for (i in unique(CZ_best)){ 
    v = c(v, which(CZ_best == i))
  }
  return(v)
}


line_sep <- function(CZ_best){
  v <- c()
  clusters <- sort(unique(CZ_best))  
  v[1] <- sum(CZ_best == clusters[1])
  
  for (i in 2:length(clusters)){
    v[i] = v[i-1] + sum(CZ_best == clusters[i])
  }
  
  return (v)
}

clu_num <- function(CZ_best){
  v <- c()
  clusters <- sort(unique(CZ_best)) 
  for (i in clusters){
    v = c(v, sum(CZ_best == i)) 
  }
  return(v)
}


# strenghten the signals
strong_signal <- function(raw_data, CX_best, rowsep, n_pick) {
  n_cluster <- length(unique(CX_best))
  b_vec <- c(0, cumsum(rowsep))
  
  for (i in 1:n_cluster) {
    j <- (b_vec[i] + 1):b_vec[i + 1]
    data <- raw_data[j, , drop = FALSE]
    

    actual_n_pick <- min(n_pick, nrow(data))
    
    if (actual_n_pick < n_pick) {
      warning("Not enough rows to sample from in cluster ", i, ". Using ", actual_n_pick, " rows instead.")
    }
    
    for (k in 1:nrow(data)) {
      select <- sample(1:nrow(data), size = actual_n_pick, replace = FALSE)
      data[k, ] <- colMeans(data[select, , drop = FALSE], na.rm = TRUE)
    }
    
    raw_data[j, ] <- data
  }
  
  return(raw_data)
}

strong_signal_iter <- function(raw_data, CX_best, rowsep, n_pick, n_iter = 1) {
  n_cluster <- length(unique(CX_best))
  b_vec <- c(0, cumsum(rowsep))
  
  for (iter in 1:n_iter) {
    for (i in 1:n_cluster) {
      j <- (b_vec[i] + 1):b_vec[i + 1]
      data <- raw_data[j, , drop = FALSE]
      
      actual_n_pick <- min(n_pick, nrow(data))
      if (actual_n_pick < n_pick) {
        warning("Not enough rows to sample from in cluster ", i, 
                ". Using ", actual_n_pick, " rows instead.")
      }
      
      for (k in 1:nrow(data)) {
        select <- sample(1:nrow(data), size = actual_n_pick, replace = FALSE)
        data[k, ] <- colMeans(data[select, , drop = FALSE], na.rm = TRUE)
      }
      
      raw_data[j, ] <- data
    }
  }
  
  return(raw_data)
}


#reorder the feature sequence based on the vertical signal strength
feature_signal_reorder<-function(colsep,X){
  n_cluster = length(colsep)
  z_vec = c(0,colsep)
  for(i in 1:n_cluster){
    # for the j-th block
    j = (z_vec[i]+1):z_vec[i+1]
    data = X[,j]
    jj = order(as.vector(colSums(data,na.rm = TRUE)))
    X[,j] = data[,jj]
  }
  return(X)
}


heatmap_fun_atac <- function(X, scaleyellowred, colsep, rowsep_x) {
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("X must be a numeric matrix.")
  }
  
  heatmap.2(X, dendrogram = "none", col = scaleyellowred, na.color = "yellow",
            labRow = "", labCol = "",  
            Rowv = FALSE, Colv = FALSE,
            margins = c(12, 10), 
            trace = "none",
            key = TRUE, key.xlab = "gene activity score", key.title = "",
            density.info = "none",
            xlab = "features", ylab = "cells",
            sepcolor = "blue",
            lwd = 3,
            colsep = colsep,
            rowsep = rowsep_x)
}




heatmap_fun_rna <- function(X, scaleyellowred, colsep, rowsep_x) {
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("X must be a numeric matrix.")
  }
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  par(cex.axis = 1.4, font.axis = 2, cex.lab = 1.5, font.lab = 2)
  
  heatmap.2(X,
            dendrogram = "none",
            col = scaleyellowred,
            na.color = "yellow",
            labRow = "", labCol = "",  
            Rowv = FALSE, Colv = FALSE,
            margins = c(12, 10),
            trace = "none",
            key = TRUE,
            key.xlab = "log(UMI + 1)",
            key.title = "",
            keysize = 1.5,
            key.par = list(font = 2),
            density.info = "none",
            xlab = "features",
            ylab = "cells",
            sepcolor = "blue",
            lwd = 3,
            colsep = colsep,
            rowsep = rowsep_x)
}


#re-order the label to ensure

reorder_label <- function(CX_best, old_label){
 
  unique_labels <- unique(old_label)
  label_map <- setNames(seq_along(unique_labels), unique_labels)
  
 
  CX0 <- sapply(CX_best, function(x) label_map[as.character(x)])
  
  return(CX0)
}



heatmap_cent <- function(raw_data_y, CZ_best){
  raw_data_y_center2 <- raw_data_y
  raw_data_y_center2 <- scale(t(raw_data_y_center2), center=T, scale=F)
  raw_data_y_center2 <- t(raw_data_y_center2)
  raw_data_y_center2 <- scale(raw_data_y_center2, center=T, scale=F)
  thres_high <- quantile(raw_data_y_center2, 0.9, na.rm=T) 
  thres_low <- quantile(raw_data_y_center2, 0.1, na.rm=T)
  raw_data_y_center2[which(raw_data_y_center2 >= thres_high)] <- thres_high
  raw_data_y_center2[which(raw_data_y_center2 <= thres_low)] <- thres_low
 
  pdf("heatmap_output.pdf", width=12, height=10)  
  heatmap.2(raw_data_y_center2, dendrogram="none", col = scaleyellowred, na.color="gray",
            labRow = "", labCol = "",
            Rowv=FALSE, Colv = FALSE,
            margins = c(5, 5), 
            trace = "none",
            key=T, key.xlab="log2(normalized read count + 1)", key.title="",
            density.info = "none",
            xlab = "features",
            ylab = "cells",
            sepcolor="blue",
            lwd=3,
            colsep=line_sep(CZ_best),
            rowsep=rowsep_y, lty=2
  )
  dev.off() 
}

