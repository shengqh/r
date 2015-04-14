setwd("H:/shengquanhu/projects/rcpa/TurboRaw2Mgf/ions_0.01_final")
plot(data$Frequency, data$MedianIntensity)
dataHighFrequency<-data[data$Frequency > 0.2,]
plot(dataHighFrequency$Frequency, dataHighFrequency$MedianIntensity)
hist(dataHighFrequency$MedianIntensity, breaks=50)
