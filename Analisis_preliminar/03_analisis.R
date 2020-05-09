# ----- Ruta de trabajo -----
set.seed(1991)

# Windows
setwd("C:/Users/dredondo/Dropbox/Transporte_interno/Máster/Ciencia de Datos/")

# Mac
#setwd("/Users/daniel/Dropbox/Transporte_interno/Máster/Ciencia de Dato/")

setwd("TFM/Análisis preliminar/descarga_20200509/gdc_download_20200509_173000.089392/")

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

# ----- Unificación de ficheros en formato matriz -----

tic("countsToMatrix")
countsInformation <- countsToMatrix("SamplesDataFrame.csv", extension = "counts")
toc()

# ----- La matriz de cuentas y las etiquetas se guardan -----

countsMatrix <- countsInformation$countsMatrix
labels <- countsInformation$labels

# ----- Se descargan los nombres de los genes ----- 

myAnnotation <- getGenesAnnotation(rownames(countsMatrix))

# ----- Cálculo de matriz de expresión de genes -----

tic("calculateGeneExpressionValues")
expressionMatrix <- calculateGeneExpressionValues(countsMatrix, myAnnotation, genesNames = TRUE)
toc()

# ----- Boxplots de la expresión para todos los genes -----

dataPlot(expressionMatrix, labels, mode = "boxplot", colours = c("blue", "red"), toPNG = TRUE)
dataPlot(expressionMatrix, labels, mode = "orderedBoxplot", colours = c("blue", "red"), toPNG = TRUE)

# ----- Controlando por el efecto batch -----

tic("batchEffectRemoval") 
svaMod <- batchEffectRemoval(expressionMatrix, as.factor(labels), method = "sva")
toc()

# ----- Extracción de DEG (Expresión Diferencial de Genes) -----

tic("DEGsExtraction") 
DEGsInformation <- DEGsExtraction(expressionMatrix, as.factor(labels),
                                  # Default parameters
                                  lfc = 1.0, pvalue = 0.1,
                                  # Ajuste por batchEffect
                                  svaCorrection = TRUE, svaMod = svaMod)
toc()

# Número de genes extraídos
print(nrow(DEGsInformation$DEGsMatrix))

topTable <- DEGsInformation$Table
DEGsMatrix <- DEGsInformation$DEGsMatrix

head(topTable) # Si la primera columna es negativa, infrarepresentado
# positiva = gen sobreexpresado



# Se traspone la matriz
DEGsMatrixML <- t(DEGsMatrix)

# ----- Representación de DEG -----

# Boxplots para todas las muestras de los primeros genes
dataPlot(DEGsMatrix[1:11, ], labels, mode = "genesBoxplot", toPNG = TRUE)

# Mapa de calor para todas las muestras de los primeros genes
dataPlot(DEGsMatrix[1:11, ], labels, mode = "heatmap", toPNG = TRUE)

# ----- Partición entrenamiento-test -----

# Partición 50% / 50% con balanceo de clase
indices <- createDataPartition(SamplesDataFrame$Class, p = .5, list = FALSE)
particion <- list(training = DEGsMatrixML[indices, ], test = DEGsMatrixML[-indices, ])

# Conjuntos
particion.entrenamiento <- particion$training
particion.test <- particion$test

# Etiquetas
labels_train <- SamplesDataFrame$Class[indices]
labels_test  <- SamplesDataFrame$Class[-indices]

# Número de casos
# Train
table(labels_train)
# Test
table(labels_test)

# Verificar balanceo de clase en entrenamiento y test
# Train
labels_train %>% table %>% prop.table %>% round(3) * 100
# Test
labels_test %>% table %>% prop.table %>% round(3) * 100

# ----- Selección de características -----

# Método mRMR (mínima redundancia, máxima relevancia)
mrmrRanking <- featureSelection(particion.entrenamiento, labels_train, colnames(particion.entrenamiento),
                                mode = "mrmr")

# Método random forest
rfRanking <- featureSelection(particion.entrenamiento, labels_train, colnames(particion.entrenamiento),
                              mode = "rf")

# Método Disease association ranking (en base a scores obtenidos en la literatura)
daRanking <- featureSelection(particion.entrenamiento, labels_train, colnames(particion.entrenamiento),
                              mode = "da", disease = "pancreatic cancer")

# -----  Sobreescribir la función dataPlot con la nueva función que pinta líneas discontinuas -----

source("../../Funciones_actualizadas_KnowSeq/dataPlot.R")

# ----- Resultados de SVM con validación cruzada para cada método de selección de características -----

numero_folds <- 2

# mRMR
results_cv_svm_mrmr <- svm_CV(particion.entrenamiento, labels_train, names(mrmrRanking),
                              numFold = numero_folds)
results_cv_svm_mrmr$bestParameters

# random forest
results_cv_svm_rf <- svm_CV(particion.entrenamiento, labels_train, rfRanking,
                            numFold = numero_folds)
results_cv_svm_rf$bestParameters

# disease association
results_cv_svm_da <- svm_CV(particion.entrenamiento, labels_train, names(daRanking),
                            numFold = numero_folds)
results_cv_svm_da$bestParameters

# ----- Resultados de SVM gráficamente -----

# Plotting the accuracy of all the folds evaluated in the CV process
#png(filename = "MRMR.png")
dataPlot(results_cv_svm_mrmr$accMatrix[, 1:11], mode = "classResults",
         main = "mRMR - Accuracy for each fold", xlab = "Genes", ylab = "Accuracy", toPNG = TRUE)
#dev.off()

# Plotting the accuracy of all the folds evaluated in the CV process
dataPlot(results_cv_svm_rf$accMatrix[, 1:9], colours = rainbow(numero_folds), mode = "classResults",
         main = "rf - Accuracy for each fold", xlab = "Genes", ylab = "Accuracy", toPNG = TRUE)

# Plotting the accuracy of all the folds evaluated in the CV process
dataPlot(results_cv_svm_da$accMatrix[, 1:9], colours = rainbow(numero_folds), mode = "classResults",
         main = "da - Accuracy for each fold", xlab = "Genes", ylab = "Accuracy", toPNG = TRUE)

# --- Mejor método en CV ---
genes_a_usar <- c(1, 2, 3, 5)

for(i in 1:length(genes_a_usar)){
  precisiones <- c(mean(results_cv_svm_mrmr$accMatrix[, genes_a_usar[i]]),
                   mean(results_cv_svm_rf$accMatrix[, genes_a_usar[i]]),
                   mean(results_cv_svm_da$accMatrix[, genes_a_usar[i]]))
  
  names(precisiones) <- c("MRMR", "RF", "DA")
  cat(paste("Precisiones para", genes_a_usar[i], "gen(es):\n"))
  print(precisiones)  
  cat("-----------------------------\n")
}

# ----- Preprocesamiento de test para aplicar SVM -----

#save.image(file = "entorno.RData")
#load("entorno.RData")

# ----- Resultados de SVM train-test -----

# mRMR
results_svm_mrmr <- svm_test(train = particion.entrenamiento, labels_train,
                             test = particion.test, labels_test, names(mrmrRanking),
                             bestParameters = results_cv_svm_mrmr$bestParameters)

# random forest
results_svm_rf <- svm_test(train = particion.entrenamiento, labels_train,
                           test = particion.test, labels_test, rfRanking,
                           bestParameters = results_cv_svm_rf$bestParameters)

# disease association
results_svm_da <- svm_test(train = particion.entrenamiento, labels_train,
                           test = particion.test, labels_test, names(daRanking),
                           bestParameters = results_cv_svm_da$bestParameters)

# ----- Matrices de confusión -----

# Matriz de confusión
results_svm_mrmr$cfMats[[1]]$table
results_svm_rf$cfMats[[1]]$table
results_svm_da$cfMats[[1]]$table

# Gráficamente
dataPlot(results_svm_mrmr$cfMats[[1]]$table, labels_test,
         mode = "confusionMatrix")
dataPlot(results_svm_rf$cfMats[[1]]$table, labels_test,
         mode = "confusionMatrix")
dataPlot(results_svm_da$cfMats[[1]]$table, labels_test,
         mode = "confusionMatrix")

# ----- Mejor método hallado en validación cruzada (rf, 5 genes) -----

dataPlot(results_svm_rf$cfMats[[5]]$table, labels_test,
         mode = "confusionMatrix")

# ----- Descarga de información sobre enfermedades relacionadas con DEGs -----

# Se recuperan las 20 enfermedades más importantes vinculadas
# Se usan los 5 genes más relevantes usando rf, el mejor método encontrado anteriormente
diseases <- DEGsToDiseases(rfRanking[1:5], getEvidences = TRUE, size = 20)
diseases

# Extracción de todas las enfermedades relacionadas con los 5 genes
enfermedades <- c()
for(i in 1:5){
  # Se extraen las enfermedades relacionadas
  enfermedades <- c(enfermedades, diseases[[i]]$summary[, 1])
  # Se eliminan duplicados
  enfermedades <- unique(enfermedades)
}
enfermedades

# Buscar cáncer
enfermedades[grep(pattern = "cancer", x = enfermedades)]
enfermedades[grep(pattern = "neopl", x = enfermedades)]
enfermedades[grep(pattern = "panc", x = enfermedades)]

# ----- Ontología de genes -----

# labelsGo <- gsub("Solid Tissue Normal",0,labels)
# labelsGo <- gsub("Primary Tumor",1,labelsGo)
# GOsMatrix <- geneOntologyEnrichment(rownames(DEGsMatrix)[1], labelsGo)
# str(GOsMatrix, max.level = 10)


# ----- Pathway visualization (no encuentra ningún pathway, la extraigo de KEGG finalmente) -----
# Annotation_counts <- getAnnotationFromEnsembl(rownames(countsMatrix))
# Annotation_DEGsMatrix <- Annotation_counts %>% filter(external_gene_name %in% rownames(DEGsMatrix))
# 
# expressionMatrix2 <- expressionMatrix
# Annotation_expressionMatrix <- Annotation_counts %>% filter(external_gene_name %in% rownames(expressionMatrix2))
# 
# DEGsPathwayVisualization(DEGsMatrix = DEGsMatrix, DEGsAnnotation = Annotation_DEGsMatrix,
#                          expressionMatrix = expressionMatrix2, expressionAnnotation = Annotation_expressionMatrix,
#                          labels = labels)
# 