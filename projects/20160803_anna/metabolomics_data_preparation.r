info<-read.csv("Brittain Metadata RTI.csv")
d0<-read.csv("Hemnes combined final.csv")
lod<-read.csv("LOD.csv")


d1<-merge(info, d0, by.x="ID", by.y="Sample.Identification", incomparables = NA)

# patient PPH189CM5264 condition  correction: 07/15/2016

d1$PAH<-ifelse(d1$ID=="PPH189CM5264", 1, d1$PAH)

# d1[1:3, 1:20]

mat0<-d1[, 16:ncol(d1)]

# by column, adjust for LOD levels by replacing with LOD/2 - use calc or SOP value??? currently using the nonzero one of the two
for (i in 1:ncol(mat0)){
  mat0[,i]<-replace(as.vector(mat0[, i]), as.vector(mat0[, i])=="<LOD", (lod$LOD.sum[i])/2)
}

uniquevalue<-sapply(mat0, function(x) length(unique(x)))
mat1<-subset(mat0, select=uniquevalue>1)

d2<-cbind(d1[,1:15], log2(data.matrix(mat1)))

#====================================================================
## univariable

library(rms)

## multiple imputation for age and BMI

a<-aregImpute(~Age + Gender + BMI + NoPH + PVH + Combined + PAH + CTD.PAH + UMC, n.impute=20, nk=0, data=d2) 

#==============================================================================

d3<-d2
#if you didn’t already, I would lump the CTD_PAH with the PAH subjects for these initial analyses.


d3$PAH_lump<-ifelse((d3$PAH==1 | d3$CTD.PAH==1), 1, 0)


# 1. PVH vs. PAH; n=58

d4.PVH<-subset(d3, PVH==1); d4.PVH$PVH_PAH_lump<-'PVH'
d4.PAH<-subset(d3, PAH_lump==1); d4.PAH$PVH_PAH_lump<-'PAH'
d4<-rbind(d4.PVH, d4.PAH)


# 2. PAH vs. NoPH

d5.NoPH<-subset(d3, NoPH==1); d5.NoPH$PAH_NoPH<-'NoPH'
d5.PAH<-subset(d3, PAH_lump==1); d5.PAH$PAH_NoPH<-'PAH'
d5<-rbind(d5.NoPH, d5.PAH)

# 3. Combined vs. PVH 

d6.Combined<-subset(d3, Combined==1); d6.Combined$PVH_Combined<-'Combined'
d6.PVH<-subset(d3, PVH==1); d6.PVH$PVH_Combined<-'PVH'
d6<-rbind(d6.Combined, d6.PVH)

# 4. Combined vs. PAH…I know power will be very limited in the comparisons using the combined group.

d7.Combined<-subset(d3, Combined==1); d7.Combined$PAH_Combined<-'Combined'
d7.PAH<-subset(d3, PAH_lump==1); d7.PAH$PAH_Combined<-'PAH'
d7<-rbind(d7.Combined, d7.PAH)

# 5.  PAH vs. PVH comparison without lumping CTD_PAH in?

d8.PAH<-subset(d3, PAH==1); d8.PAH$PAH_PVH<-'PAH'
d8.PVH<-subset(d3, PVH==1); d8.PVH$PAH_PVH<-'PVH'
d8<-rbind(d8.PAH, d8.PVH)
