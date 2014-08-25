library("DupChecker")

geoDownload("GSE1478", targetDir = "h:/DupChecker/", overwrite=TRUE)
file.rename("h:/DupChecker/GSE1478", "h:/DupChecker/GSE1478-dup")
geoDownload("GSE1478", targetDir = "h:/DupChecker/", overwrite=TRUE)

datafile<-buildFileTable("h:/DupChecker/")
datafile

result<-validateFile(datafile)
result
