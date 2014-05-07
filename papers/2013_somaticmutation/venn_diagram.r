library("VennDiagram")

getintersect<-function(data, names){
  result<-data
  lapply(names, function(name){
    result<<-result[result[,name] != '',]
  })
  nrow(result)
}

drawvenn<-function(data, cnames){
  draw.quad.venn(getintersect(matched, cnames[1]),
                 getintersect(matched, cnames[2]),
                 getintersect(matched, cnames[3]),
                 getintersect(matched, cnames[4]),
                 
                 getintersect(matched, cnames[c(1:2)]),
                 getintersect(matched, cnames[c(1,3)]),
                 getintersect(matched, cnames[c(1,4)]),
                 getintersect(matched, cnames[c(2,3)]),
                 getintersect(matched, cnames[c(2,4)]),
                 getintersect(matched, cnames[c(3,4)]),

                 getintersect(matched, cnames[c(1,2,3)]),
                 getintersect(matched, cnames[c(1,2,4)]),
                 getintersect(matched, cnames[c(1,3,4)]),
                 getintersect(matched, cnames[c(2,3,4)]),
                 
                 getintersect(matched, cnames[c(1,2,3,4)]),
                 
                 category=cnames
  )
}

matched<-read.table("H:/shengquanhu/projects/somaticmutation/all_muTect_varscan2_rsmc.site.tsv",sep='\t',header=T)
matched$DNA_M_V<-(matched$DNA_MUTECT != '' & matched$DNA_VARSCAN2 != '')
matched$DNA_M_V[!matched$DNA_M_V]<-''

drawvenn(matched, colnames(matched)[c(5:8)])
drawvenn(matched, colnames(matched)[c(5,9:11)])

drawvenn(matched, colnames(matched)[c(19,9:11)])

         names<-c("DNA_TCGA", "DNA_MUTECT", "DNA_VARSCAN2", "DNA_RSMC")

draw.quad.venn(tcga,
               dna_mutect,
               dna_varscan2,
               dna_rsmc,

               getintersect(matched, c("DNA_TCGA", "DNA_MUTECT")),
               getintersect(matched, c("DNA_TCGA", "DNA_VARSCAN2")),
               getintersect(matched, c("DNA_TCGA", "DNA_RSMC")),
               getintersect(matched, c("DNA_MUTECT", "DNA_VARSCAN2")),
               getintersect(matched, c("DNA_MUTECT", "DNA_RSMC")),
               getintersect(matched, c("DNA_VARSCAN2", "DNA_RSMC")),
               
               getintersect(matched, c("DNA_TCGA", "DNA_MUTECT", "DNA_VARSCAN2")),
               getintersect(matched, c("DNA_TCGA", "DNA_MUTECT", "DNA_RSMC")),
               getintersect(matched, c("DNA_TCGA", "DNA_VARSCAN2", "DNA_RSMC")),
               getintersect(matched, c("DNA_MUTECT", "DNA_VARSCAN2", "DNA_RSMC")),

               getintersect(matched, c("DNA_TCGA", "DNA_MUTECT", "DNA_VARSCAN2", "DNA_RSMC")),
               
               category=c("DNA_TCGA", "DNA_MUTECT", "DNA_VARSCAN2", "DNA_RSMC")
)
