#species=c("mouse", "human", "rat")
species=c("human", "rat")
types=c("miRNA", "tRNA")
dirs=c("miRNA_overlap", "tRNA")

for(spe in species){
  for(index in c(1:length(types))){
    dir = paste0("H:/shengquanhu/projects/vangard/VANGARD00055_guoyan_mirna_v2/", spe, "/topN_bowtie1_genome_cutadapt_pm_count_", dirs[index] , "_position/result")
    file = paste0("VANGARD00055_", spe, "_", types[index], ".position")
    
    setwd(dir)
    
    positions=read.table(file, sep="\t", header=T)
    
    files=unique(positions$File)
    
    cols=c(rgb(139, 35, 35, 0,maxColorValue=255), 
           rgb(139, 35, 35, 255*0.2,maxColorValue=255), 
           rgb(139, 35, 35,, 255*0.4,maxColorValue=255), 
           rgb(139, 35, 35,255*0.6,maxColorValue=255), 
           rgb(139, 35, 35,255*0.8,maxColorValue=255), 
           rgb(139, 35, 35,255,maxColorValue=255))
    
    for(file in files){
      #  file="2570-KCV-01-19"
      fp=positions[positions$File==file,]
      features=as.character(unique(fp$Feature))
      if(length(features) > 100){
        features = features[1:100]
        fp=fp[fp$Feature %in% features,]
      }
      counts=lapply(features, function(x){
        idx=which(fp$Feature==x)
        unlist(fp$Count[idx])[1]
      })
      counts=round(unlist(counts))
      strands=lapply(features, function(x){
        idx=which(fp$Feature==x)
        unlist(fp$Strand[idx])[1]
      })
      strands=unlist(strands)
      fm=data.frame(feature=features, count=counts, strand=strands)
      tnames=apply(as.matrix(fm),1,function(x){
        paste0(x[1], "(", as.numeric(x[2]),"):",x[3])
      })
      index = nrow(fm) + 1 - as.numeric(rownames(fm))
      names(index) = features
      fpf=as.character(fp$Feature)
      
      fp$Index = index[fpf]
      fp$Color=cols[round(fp$Percentage * 5)+1]
      
      png(paste0(file, ".png"), width=8000, height=6000, res=300)
      plot(fp$Position, fp$Index, pch=19, cex=fp$Percentage, xlim=c(-30, 120), col=fp$Color, xlab="Position of tRNA", ylab="tRNA name and Count", yaxt="n", bty="n", xaxt="n", main=file)
      axis(side=1, at=c(-10:120), las=2)
      lines(x=c(1,1), y=c(0, nrow(fm)),col="blue")
      text(rep(-28,nrow(fm)), index, tnames, adj=c(0,0.5))
      dev.off()
    }
  }
}