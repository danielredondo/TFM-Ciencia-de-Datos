# ----- Ruta de trabajo -----

# Windows
setwd("C:/Users/dredondo/Dropbox/Transporte_interno/Máster/Ciencia de Datos/TFM/Analisis_preliminar/descarga_20200509")
# Mac
#setwd("/Users/daniel/Dropbox/Transporte_interno/Máster/Ciencia de Datos/TFM/Analisis_preliminar/descarga_2020050")

# ----- Carga de paquetes -----

library(BiocManager) # Para instalar KnowSeq
library(KnowSeq)     # Para trabajar con genes - instalado a partir de .tar.gz de SWAD
library(tictoc)      # Para medir tiempos con tic() y toc() a lo MATLAB 
library(R.utils)     # Para gunzip (descompresión de .gz)
library(dplyr)       # Para select, filter, pipes, ...
library(beepr)       # Para avisar con beeps cuando acaba un proceso
library(caret)       # Para ml
library(e1071)       # Para svm
library(reshape)     # Para melt
library(gplots)      # Para heatmaps

# ----- Preprocesamiento para adecuar ficheros a KnowSeq -----

# Lectura de sample_sheet
samplesInfo <- read.table("gdc_sample_sheet.2020-05-09.tsv", sep = "\t", header = T)

# Se elimina ".counts.gz" de samplesInfo, porque countsToMatrix luego añade la extensión
samplesInfo[,2] <- as.character(samplesInfo[,2])
for(i in 1:nrow(samplesInfo)){
  samplesInfo[i,2] <- substr(samplesInfo[i,2], 1, nchar(samplesInfo[i, 2]) - 10)
}

# Comprobación 
head(samplesInfo)

# Definición de parámetros
Run <- samplesInfo$File.Name
Path <- samplesInfo$File.ID
Class <- samplesInfo$Sample.Type

table(Class)

# Los casos metastásicos se recodifican a "Primary Tumor", aunque no sea lo ideal...
Class <- gsub("Metastatic", "Primary Tumor", Class)
table(Class)

# Creación de dataframe
SamplesDataFrame <- data.frame(Run, Path, Class)

# Exportación a CSV de SamplesDataFrame
setwd(dir = "gdc_download_20200509_173000.089392//")
write.csv(SamplesDataFrame, file = "SamplesDataFrame.csv")
