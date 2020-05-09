# ----- Ruta de trabajo -----

# Windows
setwd("C:/Users/dredondo/Dropbox/Transporte_interno/Máster/Ciencia de Datos/TFM/Análisis preliminar/")
# Mac
#setwd("/Users/daniel/Dropbox/Transporte_interno/Máster/Ciencia de Datos/TFM/Análisis preliminar")

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

# 1. Descargar fichero de GDC Portal tras seguir guión 1 con cáncer de páncreas

# 2. Descomprimir fichero .tar (a mano): Extraer en carpeta nueva

# 3. Descomprimir distintos ficheros comprimidos

# ----- Descompresión masiva de ficheros descargados de GDC (sólo la primera vez) -----

setwd("descarga_20200509/gdc_download_20200509_173000.089392/")

# Descomprimir todos los ficheros
nombres <- list.files()
numero_ficheros <- length(nombres) - 1 # Para quitar el manifest.txt

for(i in 1:numero_ficheros){
  setwd(paste0(getwd(), "/", nombres[i]))
  nombre_fichero_comprimido <- list.files(pattern = "*.gz")
  gunzip(nombre_fichero_comprimido)
  # Medidor
  if(i %% 5 == 0) print(paste0(i, "/", numero_ficheros))
  setwd("../")
}

# Volver a carpeta de trabajo original
setwd("../..")