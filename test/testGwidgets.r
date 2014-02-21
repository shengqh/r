library("gWidgets")
library("RnaSeqSampleSize")

win <- gwindow("RNASeq Sample Size Estimation 1.0.0",width=1024, height=768, visible=FALSE)
topgroup<-ggroup(cont=win,horizontal=F)
fontlist <- list(size=14)
redfontlist <- list(unlist(fontlist), color="red")

parameters <- gframe("Parameters", container=topgroup, horizontal=F)
font(parameters)<-fontlist
tbl <- glayout(container = parameters)

tbl[1,1] <- glabel("m:",container = tbl)
tbl[1,2] <- (mEdit<-gedit("10000", container=tbl, font.attr=fontlist))
tbl[1,3] <- glabel("total number of genes for testing",container = tbl, font.attr=fontlist)

tbl[2,1] <- glabel("m1:",container = tbl, font.attr=fontlist)
tbl[2,2] <- (m1Edit<-gedit("100", container=tbl, font.attr=fontlist))
tbl[2,3] <- glabel("expected number of prognostic genes",container = tbl, font.attr=fontlist)

tbl[3,1] <- glabel("power:",container = tbl, font.attr=fontlist)
tbl[3,2] <- (powerEdit<-gedit("0.8", container=tbl, font.attr=fontlist))
tbl[3,3] <- glabel("power to detect prognostic genes",container = tbl, font.attr=fontlist)

tbl[4,1] <- glabel("f:",container = tbl)
tbl[4,2] <- (fEdit<-gedit("0.1", container=tbl))
tbl[4,3] <- glabel("FDR level",container = tbl)

tbl[5,1] <- glabel("w:",container = tbl)
tbl[5,2] <- (wEdit<-gedit("1.0", container=tbl))
tbl[5,3] <- glabel("ratio of normalization factors between two groups",container = tbl)

tbl[6,1] <- glabel("rho:",container = tbl)
tbl[6,2] <- (rhoEdit<-gedit("2.0", container=tbl))
tbl[6,3] <- glabel("minimum fold changes for prognostic genes in control group",container = tbl)

tbl[7,1] <- glabel("lambda0:",container = tbl)
tbl[7,2] <- (lambda0Edit<-gedit("5", container=tbl))
tbl[7,3] <- glabel("average read counts for prognostic gene in control group",container = tbl)

tbl[8,1] <- glabel("phi0:",container = tbl)
tbl[8,2] <- (phi0Edit<-gedit("1", container=tbl))
tbl[8,3] <- glabel("dispersion for prognostic genes in control group",container = tbl)

lapply(1:8, function(i){
  lapply(1:3, function(j){
    font(tbl[i,j])<-fontlist
  })
})

samplesize <- gframe("Calculation Result", container=topgroup, horizontal=F)
resultEdit<-gtext("Estimated sample size=", container=samplesize, expand=TRUE, font.attr=fontlist)
enabled(resultEdit)<-FALSE

description <- gframe("Description", container=topgroup, horizontal=F, expand=TRUE)
dText<-gtext("", container=description, expand=TRUE, font.attr=fontlist)
enabled(dText)<-FALSE

doCalculation<-function(h, ...){
  svalue(resultEdit)<-"Estimated sample size = ???"
  font(resultEdit)<-fontlist
  ret<-sample_size(as.numeric(svalue(mEdit)), 
                   as.numeric(svalue(m1Edit)), 
                   as.numeric(svalue(powerEdit)),
                   as.numeric(svalue(fEdit)),
                   as.numeric(svalue(wEdit)),
                   as.numeric(svalue(rhoEdit)),
                   as.numeric(svalue(lambda0Edit)),
                   as.numeric(svalue(phi0Edit)))
  svalue(dText)<-paste0("We are planning a study to identify differential gene expression between two groups. Prior data indicate that ",
    "the minimum average read counts among the prognostic genes in the control group is ", svalue(lambda0Edit), 
    ", the maximum dispersion is ", svalue(phi0Edit),
    ", and the ratio of the geometric mean of normalization factors is ", svalue(wEdit),
    ". Suppose that the total number genes for testing is ", svalue(mEdit),
    " and the top ", svalue(m1Edit), 
    " genes are prognostic. If the desired minimum fold changes is ", svalue(rhoEdit),
    ", we will need to study ", ret,
    " subjects to be able to reject the null hypothesis that the fold changes is 1 with probability (power) ", svalue(powerEdit),
    " (", as.numeric(svalue(m1Edit)) * as.numeric(svalue(powerEdit)), "/", svalue(m1Edit),
    ") using exact test. The FDR associated with this test of this null hypothesis is ", svalue(fEdit) ,".")
  svalue(resultEdit)<-paste0("Estimated sample size = ", ret)
  font(resultEdit)<-redfontlist
}

button.group<-ggroup(container=topgroup)
addSpring(button.group)
calcbutton<-gbutton("Calculate", container=button.group, handler=doCalculation)
copybutton<-gbutton("Copy Description", container=button.group, handler = function(h, ...) writeClipboard(svalue(dText)) )
closebutton<-gbutton("Close", container=button.group, handler = function(h,...) dispose(win))

font(calcbutton)<-fontlist
font(copybutton)<-fontlist
font(closebutton)<-fontlist

visible(win)<-TRUE
