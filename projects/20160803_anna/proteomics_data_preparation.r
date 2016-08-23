info<-read.csv("Brittain Metadata RTI.csv")
d0<-read.delim("proteomics_expression_data.txt", header=T, stringsAsFactors = F)
d0<-d0[,c(3,12:ncol(d0))]
d0$ID_<-gsub("\\s+","", d0$ID_)
d0$ID_<-gsub("C15T0061212","C15-T0061212", d0$ID_)

#d0$ID_[!(d0$ID_ %in% info$ID)]

d1<-merge(info, d0, by.x="ID", by.y="ID_", incomparables = NA)

# patient PPH189CM5264 condition  correction: 07/15/2016

d1$PAH<-ifelse(d1$ID=="PPH189CM5264", 1, d1$PAH)

# d1[1:3, 1:20]

mat0<-d1[, 11:ncol(d1)]

d2<-cbind(d1[,1:10], log2(data.matrix(mat0)))

d3<-d2

# 1. PVH vs. PAH

d4.PVH<-subset(d3, PVH==1); d4.PVH$PVH_PAH<-'PVH'
d4.PAH<-subset(d3, PAH==1); d4.PAH$PVH_PAH<-'PAH'
d4<-rbind(d4.PVH, d4.PAH)

# 2. PAH vs. NoPH

d5.NoPH<-subset(d3, NoPH==1); d5.NoPH$PAH_NoPH<-'NoPH'
d5.PAH<-subset(d3, PAH==1); d5.PAH$PAH_NoPH<-'PAH'
d5<-rbind(d5.NoPH, d5.PAH)

# 3. Combined vs. PVH 

d6.Combined<-subset(d3, Combined==1); d6.Combined$PVH_Combined<-'Combined'
d6.PVH<-subset(d3, PVH==1); d6.PVH$PVH_Combined<-'PVH'
d6<-rbind(d6.Combined, d6.PVH)

# 4. Combined vs. PAH…I know power will be very limited in the comparisons using the combined group.

d7.Combined<-subset(d3, Combined==1); d7.Combined$PAH_Combined<-'Combined'
d7.PAH<-subset(d3, PAH==1); d7.PAH$PAH_Combined<-'PAH'
d7<-rbind(d7.Combined, d7.PAH)
