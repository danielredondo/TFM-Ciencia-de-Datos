# ----- Semilla aleatoria para reproducibilidad -----
rm(list = ls())
set.seed(31415)

# Número de genes a seleccionar por los métodos de selección de características
numero_de_genes <- 10

# ----- Ruta de trabajo -----

# Windows
#setwd("C:/Users/dredondo/Dropbox/Transporte_interno/Máster/Ciencia de Datos/")

# Mac
setwd("/Users/daniel/Dropbox/Transporte_interno/Máster/Ciencia de Datos/")

# Carpeta de datos
setwd("TFM/Analisis_higado/data/")

# ----- Carga de paquetes -----

# Instalación de KnowSeq: (es una versión fija de GitHub)
# devtools::install_github("CasedUgr/KnowSeq", ref = "f59cb9e1cb02702697c208cf2c61c45d6e0b7a08", force = TRUE) 
# Si hay problemas del tipo "Error: (converted from warning)" : Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
library(KnowSeq)     # Para trabajar con datos de transcriptómica de GDC Portal
library(dplyr)       # Para select, filter, pipes, ...
library(tictoc)      # Para medir tiempos con tic() y toc() a lo MATLAB 
library(beepr)       # Para avisar con beeps cuando acaba un proceso
library(caret)       # Para machine learning 
library(e1071)       # Para SVM
library(gplots)      # Para heatmaps
library(reshape2)    # Para melt
library(ggalluvial)  # Para diagrama de Sankey

# -----  Sobreescribir la función dataPlot con una nueva función que pinta líneas discontinuas -----

source("../funciones_personalizadas_knowseq/dataPlot.R")

# ----- Preprocesamiento para adecuar ficheros a KnowSeq -----

# Lectura de sample_sheet
samplesInfo <- read.table("gdc_sample_sheet.2020-06-11.tsv", sep = "\t", header = T)

# Se elimina ".counts.gz" de samplesInfo, porque countsToMatrix luego añade la extensión
samplesInfo[,2] <- as.character(samplesInfo[,2])
for(i in 1:nrow(samplesInfo)){
  samplesInfo[i,2] <- substr(samplesInfo[i,2], 1, nchar(samplesInfo[i, 2]) - 10)
}

# Comprobación de importación correcta
head(samplesInfo)

# Definición de parámetros
Run <- samplesInfo$File.Name
Path <- samplesInfo$File.ID
Class <- samplesInfo$Sample.Type

table(Class)

# Los casos que no son "Primary Tumor" o "Solid Tissue Normal" se eliminan de:
#    - El fichero gdc_sample_sheet
#    - El fichero gdc_manifest
#    - La carpeta de gdc_download 
Path[which(Class %in% c("Metastatic", "Recurrent Tumor"))]

# Creación de dataframe
SamplesDataFrame <- data.frame(Run, Path, Class)

# Exportación a CSV de SamplesDataFrame
setwd(dir = "gdc_download/")
write.csv(SamplesDataFrame, file = "03_SamplesDataFrame.csv")

# ----- Unificación de ficheros en formato matriz -----

tic("countsToMatrix") # 85 segundos
countsInformation <- countsToMatrix("03_SamplesDataFrame.csv", extension = "counts")
toc()

# ----- La matriz de cuentas y las etiquetas se guardan -----

countsMatrix <- countsInformation$countsMatrix
labels <- countsInformation$labels

# ----- Se descargan los nombres de los genes ----- 

myAnnotation <- getGenesAnnotation(rownames(countsMatrix))

# ----- Cálculo de matriz de expresión de genes -----

tic("calculateGeneExpressionValues") # 82 segundos
expressionMatrix <- calculateGeneExpressionValues(countsMatrix, myAnnotation, genesNames = TRUE)
toc()

# ----- Controlando por el efecto batch -----

tic("batchEffectRemoval") # 230 segundos
svaMod <- batchEffectRemoval(expressionMatrix, as.factor(labels), method = "sva")
toc()

# ----- Total - Extracción de DEG (Expresión Diferencial de Genes) -----

tic("DEGsExtraction") # 13 segundos
DEGsInformation <- DEGsExtraction(expressionMatrix, as.factor(labels),
                                  # p-valor
                                  pvalue = 0.001,
                                  # Ajuste por batchEffect
                                  svaCorrection = TRUE, svaMod = svaMod)
toc()

# Número de genes extraídos: 2274
print(nrow(DEGsInformation$DEGsMatrix))

topTable <- DEGsInformation$Table
DEGsMatrix <- DEGsInformation$DEGsMatrix

write.csv2(labels, file = "../saved_files/labels.csv", row.names = F)
write.csv2(DEGsMatrix, file = "../saved_files/DEGsMatrix.csv")

# ----- Sólo 200 genes - Extracción de DEG (Expresión Diferencial de Genes) -----

tic("DEGsExtraction") # 13 segundos
DEGsInformation <- DEGsExtraction(expressionMatrix, as.factor(labels),
                                  # p-valor
                                  #pvalue = 0.001,
                                  number = 200,
                                  # Ajuste por batchEffect
                                  svaCorrection = TRUE, svaMod = svaMod)
toc()

# Número de genes extraídos: 200
print(nrow(DEGsInformation$DEGsMatrix))

topTable <- DEGsInformation$Table
DEGsMatrix <- DEGsInformation$DEGsMatrix

write.csv2(labels, file = "../saved_files/higado_200genes_labels.csv", row.names = F)
write.csv2(DEGsMatrix, file = "../saved_files/higado_200genes_DEGsMatrix.csv")
