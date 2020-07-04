# ----- Ruta de trabajo -----

# Windows
#setwd("C:/Users/dredondo/Dropbox/Transporte_interno/Máster/Ciencia de Datos/TFM/Analisis_final_higado/")

# Mac
setwd("/Users/daniel/Dropbox/Transporte_interno/Máster/Ciencia de Datos/TFM/Analisis_final_higado/")

# ----- Carga de paquetes -----

library(R.utils)     # Para gunzip (descompresión de .gz)

# 1. Descargar fichero de GDC Portal tras seguir guión 1 con cáncer de hígado

# 2. Descomprimir fichero .tar (a mano): Extraer en carpeta nueva

# 3. Descomprimir distintos ficheros comprimidos

# ----- Descompresión masiva de ficheros descargados de GDC (sólo la primera vez) -----

setwd("data/gdc_download")

# Listar todas las subcarpetas
nombres <- list.files()
numero_ficheros <- length(nombres) - 1 # Para quitar el fichero manifest.txt

for(i in 1:numero_ficheros){
  setwd(paste0(getwd(), "/", nombres[i]))
  nombre_fichero_comprimido <- list.files(pattern = "*.gz")
  # Se descomprime el fichero comprimido y se elimina
  gunzip(nombre_fichero_comprimido, remove = TRUE)
  # Indicador
  if(i %% 5 == 0) print(paste0(i, "/", numero_ficheros))
  setwd("../")
}