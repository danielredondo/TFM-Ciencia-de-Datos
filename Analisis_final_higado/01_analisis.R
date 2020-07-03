# ----- Semilla aleatoria para reproducibilidad -----
set.seed(31415)

# ----- Ruta de trabajo -----

# Windows
#setwd("C:/Users/dredondo/Dropbox/Transporte_interno/Máster/Ciencia de Datos/")
# Mac
setwd("/Users/daniel/Dropbox/Transporte_interno/Máster/Ciencia de Datos/")

# Carpeta de datos
setwd("TFM/Analisis_final_higado/data/")

# ----- Carga de paquetes -----

library(BiocManager) # Para instalar KnowSeq
# BiocManager::install("KnowSeq")
# devtools::install_github("CasedUgr/KnowSeq")
library(KnowSeq)     # Para trabajar con datos de transcriptómica de GDC de Portal - instalado con BiocManager::install("KnowSeq")
library(dplyr)       # Para select, filter, pipes, ...
library(tictoc)      # Para medir tiempos con tic() y toc() a lo MATLAB 
library(beepr)       # Para avisar con beeps cuando acaba un proceso
library(caret)       # Para machine learning 
library(e1071)       # Para SVM
library(gplots)      # Para heatmaps
library(reshape2)    # Para melt
library(ggalluvial)  # Para diagrama de Sankey

# -----  Sobreescribir la función dataPlot con la nueva función que pinta líneas discontinuas -----

source("../Funciones_actualizadas_KnowSeq/dataPlot.R")

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
write.csv(SamplesDataFrame, file = "SamplesDataFrame.csv")

# ----- Unificación de ficheros en formato matriz -----

tic("countsToMatrix") # 85 segundos
countsInformation <- countsToMatrix("SamplesDataFrame.csv", extension = "counts")
toc()

# ----- La matriz de cuentas y las etiquetas se guardan -----

countsMatrix <- countsInformation$countsMatrix
labels <- countsInformation$labels

# ----- Se descargan los nombres de los genes ----- 

myAnnotation <- getGenesAnnotation(rownames(countsMatrix))

# ----- Cálculo de matriz de expresión de genes -----

tic("calculateGeneExpressionValues") # 110 segundos
expressionMatrix <- calculateGeneExpressionValues(countsMatrix, myAnnotation, genesNames = TRUE)
toc()

# ----- Boxplots de la expresión para todos los genes -----

png(filename = "../../01_images/00_boxplot_all.png", width = 13, height = 8, units = "in", res = 300)
dataPlot(expressionMatrix, labels, mode = "boxplot", colours = c("red", "green"))
legend(1, 29, c("Primary Tumor", "Solid Tissue Normal"),
       fill = c("red", "green"), xpd = TRUE)
title("Boxplot de matriz de expresión para todos los genes")
dev.off()

png(filename = "../../01_images/01_boxplot_all_ordered.png", width = 13, height = 8, units = "in", res = 300)
dataPlot(expressionMatrix, labels, mode = "orderedBoxplot", colours = c("red", "green"))
legend(1, 29, c("Primary Tumor", "Solid Tissue Normal"),
       fill = c("red", "green"), xpd = TRUE)
title("Boxplot de matriz de expresión para todos los genes")
dev.off()

# ----- Controlando por el efecto batch -----

tic("batchEffectRemoval") # 263 segundos
svaMod <- batchEffectRemoval(expressionMatrix, as.factor(labels), method = "sva")
toc()

# ----- Extracción de DEG (Expresión Diferencial de Genes) -----

tic("DEGsExtraction") # 20 segundos
DEGsInformation <- DEGsExtraction(expressionMatrix, as.factor(labels),
                                  # p-valor
                                  pvalue = 0.001,
                                  # Número de genes
                                  #number = 12,
                                  # Ajuste por batchEffect
                                  svaCorrection = TRUE, svaMod = svaMod)
toc()

# Número de genes extraídos: 2268
print(nrow(DEGsInformation$DEGsMatrix))

topTable <- DEGsInformation$Table
DEGsMatrix <- DEGsInformation$DEGsMatrix

head(topTable) # Si la primera columna es negativa, infraexpresado
               # positiva = gen sobreexpresado

# Se traspone la matriz
DEGsMatrixML <- t(DEGsMatrix)

# ----- Representación de DEG -----

# Boxplots para todas las muestras de los primeros genes
png(filename = "../../01_images/02_boxplot_primeros_12_genes.png", width = 13, height = 8, units = "in", res = 300)
dataPlot(DEGsMatrix[1:12, ], labels, mode = "genesBoxplot")
dev.off()

# Mapa de calor para todas las muestras de los primeros genes
png(filename = "../../01_images/03_heatmap_primeros_12_genes.png", width = 13, height = 8, units = "in", res = 300)
DEGsMatrix_heatmap <- DEGsMatrix[, c(which(SamplesDataFrame$Class == "Primary Tumor"),
                                     which(SamplesDataFrame$Class == "Solid Tissue Normal"))]
labels_heatmap <- labels[c(which(SamplesDataFrame$Class == "Primary Tumor"),
                           which(SamplesDataFrame$Class == "Solid Tissue Normal"))]
dataPlot(DEGsMatrix_heatmap[1:12, ], labels_heatmap, mode = "heatmap")
dev.off()

# ----- Partición entrenamiento-test -----

# Nota shiny: de aquí en adelante se usa lo siguiente:
# labels
# DEGsMatrix
#tic("Full proceso shiny")
# Para Shiny:
#load("../saved_files/matriz.RData")
#elementos <- ls()
#elementos <- elementos[which(elementos != "matriz")]
#rm(list = elementos)
#rm(elementos)
#source("../../Funciones_actualizadas_KnowSeq/dataPlot.R")
# Extraer labels
#labels <- matriz[1, ] %>% as.vector
# Extraer DEGsMatrix
#DEGsMatrix <- matriz[2:nrow(matriz), ]
#filas <- rownames(DEGsMatrix)
#DEGsMatrix <- apply(DEGsMatrix, 2, as.numeric)
#rownames(DEGsMatrix) <- filas
# Crear DEGsMatrixML
#DEGsMatrixML <- t(DEGsMatrix)
#rm(matriz)
#rm(filas)

# Partición 75% / 25% con balanceo de clase
porcentaje <- .75
indices <- createDataPartition(labels, p = porcentaje, list = FALSE)
particion <- list(training = DEGsMatrixML[indices, ], test = DEGsMatrixML[-indices, ])

# Conjuntos
particion.entrenamiento <- particion$training
particion.test <- particion$test

# Etiquetas
labels_train <- labels[indices]
labels_test  <- labels[-indices]

# Número de casos
# Train
table(labels_train)
entr_tum <- table(labels_train)[1]
entr_san <- table(labels_train)[2]

# Test
table(labels_test)
test_tum <- table(labels_test)[1]
test_san <- table(labels_test)[2]

# Total
table(labels)

# Verificar balanceo de clase en entrenamiento y test
# Train
labels_train %>% table %>% prop.table %>% round(3) * 100
# Test
labels_test %>% table %>% prop.table %>% round(3) * 100
# Total
labels %>% table %>% prop.table %>% round(3) * 100

# Diagrama de Sankey
datos_sankey <- data.frame(tipo = c(paste0("Tumor\n", entr_tum + test_tum, " casos"),
                                    paste0("Tumor\n", entr_tum + test_tum, " casos"),
                                    paste0("Tejido normal\n", entr_san + test_san, " casos"),
                                    paste0("Tejido normal\n", entr_san + test_san, " casos")),
                           traintest = c(paste0("Entrenamiento\n", entr_tum, " tumores\n", entr_san, " tejido normal"),
                                         paste0("Test\n", test_tum, " tumores\n", test_san, " tejido normal"),
                                         paste0("Entrenamiento\n", entr_tum, " tumores\n", entr_san, " tejido normal"),
                                         paste0("Test\n", test_tum, " tumores\n", test_san, " tejido normal")),
                           value = c(entr_tum, test_tum, entr_san, test_san))

# Pequeño reorden para que mejorar la presentación de los datos
datos_sankey$tipo <- factor(datos_sankey$tipo,
                            levels = c(paste0("Tumor\n", entr_tum + test_tum, " casos"), paste0("Tejido normal\n", entr_san + test_san, " casos")),
                            ordered = T)

ggplot(data = datos_sankey,
       aes(axis1 = tipo, axis2 = traintest, y = value)) +
  scale_x_discrete(limits = c("Tipo de muestra", "Partición\nentrenamiento-test"),
                   expand = c(.1, .05)) +
  ylab("") +
  geom_alluvium(col = "black", alpha = 1) +
  geom_alluvium(aes(fill = tipo), alpha = .6, show.legend = FALSE) +
  geom_stratum() +
  geom_text(stat = "stratum", infer.label = TRUE, cex = 3) +
  theme_minimal() +
  ggtitle("Partición en conjuntos de entrenamiento y test",
          "Reparto 75% - 25% con balanceo de clases") +
  theme(plot.title = element_text(hjust = .5),
        plot.subtitle = element_text(hjust = .5),
        axis.text = element_text(color = "black", margin = margin(t = -30), size = 12),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()) 

ggsave(filename = "../../01_images/04_sankey.png", width = 8, height = 7)

# ----- Selección de características: 20 genes más importantes -----

# Método mRMR (mínima redundancia, máxima relevancia)
mrmrRanking <- featureSelection(particion.entrenamiento, labels_train, colnames(particion.entrenamiento),
                                mode = "mrmr")
mrmrRanking <- names(mrmrRanking)[1:20]

# Método random forest
rfRanking <- featureSelection(particion.entrenamiento, labels_train, colnames(particion.entrenamiento),
                              mode = "rf")
rfRanking <- rfRanking[1:20]

# Método Disease association ranking (en base a scores obtenidos en la literatura)
daRanking <- featureSelection(particion.entrenamiento, labels_train, colnames(particion.entrenamiento),
                              mode = "da", disease = "liver cancer")

# Selección de los 20 genes más relevantes. Si hay más de 20 genes con relevancia = 1, se seleccionan
# todos los que tengan relevancia = 1
ifelse(sum(daRanking == 1) <= 20,
       daRanking <- names(daRanking[1:20]),
       daRanking <- names(daRanking[1:sum(daRanking == 1)]))
      
mrmrRanking
rfRanking
daRanking

# ----- SVM: Resultados con validación cruzada para cada método de selección de características -----

numero_folds <- 5

# mRMR: 24 segundos
tic("svm_mrmr")
results_cv_svm_mrmr <- svm_CV(particion.entrenamiento, labels_train, mrmrRanking,
                              numFold = numero_folds)
toc()

# random forest: 22 segundos
tic("svm_rf")
results_cv_svm_rf <- svm_CV(particion.entrenamiento, labels_train, rfRanking,
                            numFold = numero_folds)
toc()

# disease association: 40 segundos
tic("svm_da")
results_cv_svm_da <- svm_CV(particion.entrenamiento, labels_train, daRanking,
                            numFold = numero_folds)
toc()

# Tabla con mejores parámetros
mejores_parametros_svm <- rbind("mrmr" = results_cv_svm_mrmr$bestParameters,
      "rf" = results_cv_svm_rf$bestParameters,
      "da" = results_cv_svm_da$bestParameters)

print(mejores_parametros_svm)

# ----- SVM: Resultados gráficos de validación cruzada -----

# Plotting the accuracy of all the folds evaluated in the CV process
png(filename = "../../01_images/05_svm_MRMR.png", width = 13, height = 8, units = "in", res = 300)
dataPlot(results_cv_svm_mrmr$accMatrix[, 1:12], colours = rainbow(numero_folds),
         mode = "classResults",
         main = "mRMR - Accuracy for each fold", xlab = "Genes", ylab = "Accuracy")
dev.off()

# Plotting the accuracy of all the folds evaluated in the CV process
png(filename = "../../01_images/06_svm_RF.png", width = 13, height = 8, units = "in", res = 300)
dataPlot(results_cv_svm_rf$accMatrix[, 1:12], colours = rainbow(numero_folds),
         mode = "classResults",
         main = "rf - Accuracy for each fold", xlab = "Genes", ylab = "Accuracy")
dev.off()

# Plotting the accuracy of all the folds evaluated in the CV process
png(filename = "../../01_images/07_svm_DA.png", width = 13, height = 8, units = "in", res = 300)
dataPlot(results_cv_svm_da$accMatrix[, 1:12], colours = rainbow(numero_folds),
         mode = "classResults",
         main = "da - Accuracy for each fold", xlab = "Genes", ylab = "Accuracy")
dev.off()

# ----- SVM: Mejor método en CV basado en precisión, especificidad, sensibilidad y F1-score -----
genes_a_usar <- c(1:10)

# Precisión
prec <- matrix(0, nrow = length(genes_a_usar), ncol = 3)
colnames(prec) <- c("MRMR", "RF", "DA")

for(i in 1:length(genes_a_usar)){
  prec[i, ] <- 100 * c(mean(results_cv_svm_mrmr$accMatrix[, genes_a_usar[i]]),
                       mean(results_cv_svm_rf$accMatrix[, genes_a_usar[i]]),
                       mean(results_cv_svm_da$accMatrix[, genes_a_usar[i]]))
}

which(prec == max(prec), arr.ind = TRUE)

#Especificidad 
spec <- matrix(0, nrow = length(genes_a_usar), ncol = 3)
colnames(spec) <- c("MRMR", "RF", "DA")

for(i in 1:length(genes_a_usar)){
  spec[i, ] <- 100 * c(mean(results_cv_svm_mrmr$specMatrix[, genes_a_usar[i]]),
                     mean(results_cv_svm_rf$specMatrix[, genes_a_usar[i]]),
                     mean(results_cv_svm_da$specMatrix[, genes_a_usar[i]]))
}

which(spec == max(spec), arr.ind = TRUE)

#Sensibilidad
sens <- matrix(0, nrow = length(genes_a_usar), ncol = 3)
colnames(sens) <- c("MRMR", "RF", "DA")

for(i in 1:length(genes_a_usar)){
  sens[i, ] <- 100 * c(mean(results_cv_svm_mrmr$sensMatrix[, genes_a_usar[i]]),
                       mean(results_cv_svm_rf$sensMatrix[, genes_a_usar[i]]),
                       mean(results_cv_svm_da$sensMatrix[, genes_a_usar[i]]))
}

which(sens == max(sens), arr.ind = TRUE)

#F1-Score
genes_a_usar <- c(1:10)
f1 <- matrix(0, nrow = length(genes_a_usar), ncol = 3)
colnames(f1) <- c("MRMR", "RF", "DA")

for(i in 1:length(genes_a_usar)){
  f1[i, ] <- 100 * c(mean(results_cv_svm_mrmr$f1Matrix[, genes_a_usar[i]]),
            mean(results_cv_svm_rf$f1Matrix[, genes_a_usar[i]]),
            mean(results_cv_svm_da$f1Matrix[, genes_a_usar[i]]))
}

which(f1 == max(f1), arr.ind = TRUE)
# El mejor F1-score se obtiene con 2, 3, 4 genes en RF -> Se elige el modelo con menos genes: 2

# Tablas y gráficos precisión, especificidad, sensibilidad y F1-score
which(f1 == max(f1), arr.ind = TRUE)
face_etiquetas <- rep("plain", 30)
face_etiquetas[c(5,8,11)] <- "bold"

ggplot(melt(t(f1)), aes(Var1, Var2, fill = value, label = sprintf(round(value, 2), fmt = '%#.2f'))) +
  geom_tile(colour = "black") + 
  scale_fill_gradient(low = "firebrick2", high = "chartreuse") +
  labs(x = "Método de selección de características",
       y = "Número de genes usados en el modelo",
       fill = "F1-Score",
       title = "F1-Score de SVM según método de selección de características\ny número de genes usados en el modelo",
       subtitle = "Resaltado en negrita el mejor F1-Score") +
  scale_y_continuous(breaks = 1:10) + 
  geom_label(fontface = face_etiquetas,
             label.size = 0) + 
  theme_minimal() + 
  theme(plot.title = element_text(hjust = .5),
        panel.grid = element_blank(),
        plot.subtitle = element_text(hjust = .5),
        panel.spacing = unit(c(0, 0, 0, 0), "null"),
        axis.text.x = element_text(size = 11, margin = margin(t = -15)),
        panel.grid.major.x = element_blank(),
        axis.text = element_text(color = "black"))

ggsave(filename = "../../01_images/08_svm_f1_score.png", width = 8, height = 8)

# Gráfico precisión
which(prec == max(prec), arr.ind = TRUE)
face_etiquetas <- rep("plain", 30)
face_etiquetas[c(5,8,11)] <- "bold"

ggplot(melt(t(prec)), aes(Var1, Var2, fill = value, label = sprintf(round(value, 2), fmt = '%#.2f'))) +
  geom_tile(colour = "black") + 
  scale_fill_gradient(low = "firebrick2", high = "chartreuse") +
  labs(x = "Método de selección de características",
       y = "Número de genes usados en el modelo",
       fill = "Precisión",
       title = "Precisión de SVM según método de selección de características\ny número de genes usados en el modelo",
       subtitle = "Resaltado en negrita la mejor precisión") +
  scale_y_continuous(breaks = 1:10) + 
  geom_label(fontface = face_etiquetas,
             label.size = 0) + 
  theme_minimal() + 
  theme(plot.title = element_text(hjust = .5),
        panel.grid = element_blank(),
        plot.subtitle = element_text(hjust = .5),
        panel.spacing = unit(c(0, 0, 0, 0), "null"),
        axis.text.x = element_text(size = 11, margin = margin(t = -15)),
        panel.grid.major.x = element_blank(),
        axis.text = element_text(color = "black"))

ggsave(filename = "../../01_images/09_svm_accuracy.png", width = 8, height = 8)

# ----- SVM: Resultados gráficos en train con el mejor método ------

# Imaginemos que el mejor ranking es rf con 2 genes
mejores_genes_svm <- rfRanking[1:2]

# Se grafican boxplots y mapas de calor con el mejor método
png(filename = "../../01_images/10_heatmap_svm_mejor_metodo.png", width = 13, height = 8, units = "in", res = 300)
DEGsMatrix_heatmap <- DEGsMatrix[mejores_genes_svm, c(which(labels == "Primary Tumor"),
                                                which(labels == "Solid Tissue Normal"))]
labels_heatmap <- labels[c(which(labels == "Primary Tumor"),
                           which(labels == "Solid Tissue Normal"))]
dataPlot(DEGsMatrix_heatmap[mejores_genes_svm, ], labels_heatmap, mode = "heatmap")
dev.off()

png(filename = "../../01_images/11_boxplots_svm_mejor_metodo.png", width = 13, height = 8, units = "in", res = 300)
dataPlot(DEGsMatrix[mejores_genes_svm, ], labels, mode = "genesBoxplot")
dev.off()

# ----- SVM: Resultados en test del mejor modelo -----

# mRMR
results_svm_mrmr <- svm_test(train = particion.entrenamiento, labels_train,
                             test = particion.test, labels_test, mrmrRanking,
                             bestParameters = results_cv_svm_mrmr$bestParameters)

# random forest
results_svm_rf <- svm_test(train = particion.entrenamiento, labels_train,
                           test = particion.test, labels_test, rfRanking,
                           bestParameters = results_cv_svm_rf$bestParameters)

# disease association
results_svm_da <- svm_test(train = particion.entrenamiento, labels_train,
                           test = particion.test, labels_test, daRanking,
                           bestParameters = results_cv_svm_da$bestParameters)

# Matriz de confusión mejor clasificador
results_svm_rf$cfMats[[2]]$table

# Gráficamente
png(filename = "../../01_images/12_matriz_confusion_svm_mejor_metodo.png", width = 13, height = 8, units = "in", res = 300)
dataPlot(results_svm_rf$cfMats[[4]]$table, labels_test,
         mode = "confusionMatrix")
dev.off()

# ----- kNN: Resultados con validación cruzada para cada método de selección de características -----
numero_folds <- 5

# mRMR
results_cv_knn_mrmr <- knn_CV(particion.entrenamiento, labels_train, mrmrRanking,
                              numFold = numero_folds)
results_cv_knn_mrmr$bestK

# random forest
results_cv_knn_rf <- knn_CV(particion.entrenamiento, labels_train, rfRanking,
                            numFold = numero_folds)
results_cv_knn_rf$bestK

# disease association
results_cv_knn_da <- knn_CV(particion.entrenamiento, labels_train, daRanking,
                            numFold = numero_folds)
results_cv_knn_da$bestK

# ----- kNN: Resultados gráficos de validación cruzada -----

# Plotting the accuracy of all the folds evaluated in the CV process
png(filename = "../../01_images/13_knn_MRMR.png", width = 13, height = 8, units = "in", res = 300)
dataPlot(results_cv_knn_mrmr$accMatrix[, 1:12], colours = rainbow(numero_folds), mode = "classResults",
         main = "mRMR - Accuracy for each fold", xlab = "Genes", ylab = "Accuracy")
dev.off()

# Plotting the accuracy of all the folds evaluated in the CV process
png(filename = "../../01_images/14_knn_RF.png", width = 13, height = 8, units = "in", res = 300)
dataPlot(results_cv_knn_rf$accMatrix[, 1:12], colours = rainbow(numero_folds), mode = "classResults",
         main = "rf - Accuracy for each fold", xlab = "Genes", ylab = "Accuracy")
dev.off()

# Plotting the accuracy of all the folds evaluated in the CV process
png(filename = "../../01_images/15_knn_DA.png", width = 13, height = 8, units = "in", res = 300)
dataPlot(results_cv_knn_da$accMatrix[, 1:12], colours = rainbow(numero_folds), mode = "classResults",
         main = "da - Accuracy for each fold", xlab = "Genes", ylab = "Accuracy")
dev.off()

# ----- kNN: Mejor método en CV basado en precisión, especificidad, sensibilidad y F1-score -----
genes_a_usar <- c(1:10)

# Precisión
prec <- matrix(0, nrow = length(genes_a_usar), ncol = 3)
colnames(prec) <- c("MRMR", "RF", "DA")

for(i in 1:length(genes_a_usar)){
  prec[i, ] <- 100 * c(mean(results_cv_knn_mrmr$accMatrix[, genes_a_usar[i]]),
                       mean(results_cv_knn_rf$accMatrix[, genes_a_usar[i]]),
                       mean(results_cv_knn_da$accMatrix[, genes_a_usar[i]]))
}

which(prec == max(prec), arr.ind = TRUE)

# Especificidad
spec <- matrix(0, nrow = length(genes_a_usar), ncol = 3)
colnames(spec) <- c("MRMR", "RF", "DA")

for(i in 1:length(genes_a_usar)){
  spec[i, ] <- 100 * c(mean(results_cv_knn_mrmr$specMatrix[, genes_a_usar[i]]),
                       mean(results_cv_knn_rf$specMatrix[, genes_a_usar[i]]),
                       mean(results_cv_knn_da$specMatrix[, genes_a_usar[i]]))
}

which(spec == max(spec), arr.ind = TRUE)

# Sensibilidad
sens <- matrix(0, nrow = length(genes_a_usar), ncol = 3)
colnames(sens) <- c("MRMR", "RF", "DA")

for(i in 1:length(genes_a_usar)){
  sens[i, ] <- 100 * c(mean(results_cv_knn_mrmr$sensMatrix[, genes_a_usar[i]]),
                       mean(results_cv_knn_rf$sensMatrix[, genes_a_usar[i]]),
                       mean(results_cv_knn_da$sensMatrix[, genes_a_usar[i]]))
}

which(sens == max(sens), arr.ind = TRUE)

# F1-Score
f1 <- matrix(0, nrow = length(genes_a_usar), ncol = 3)
colnames(f1) <- c("MRMR", "RF", "DA")

for(i in 1:length(genes_a_usar)){
  f1[i, ] <- 100 * c(mean(results_cv_knn_mrmr$f1Matrix[, genes_a_usar[i]]),
                     mean(results_cv_knn_rf$f1Matrix[, genes_a_usar[i]]),
                     mean(results_cv_knn_da$f1Matrix[, genes_a_usar[i]]))
}

which(f1 == max(f1), arr.ind = TRUE)
# El mejor F1-score se obtiene con 2, 3, 4 genes en RF -> Se elige el modelo con menos genes: 2

# Tablas y gráficos precisión, especificidad, sensibilidad y F1-score
which(prec == max(prec), arr.ind = TRUE)
which(spec == max(spec), arr.ind = TRUE)
which(sens == max(sens), arr.ind = TRUE)
which(f1 == max(f1), arr.ind = TRUE)

face_etiquetas <- rep("plain", 30)
face_etiquetas[c(5,8,11)] <- "bold"

ggplot(melt(t(f1)), aes(Var1, Var2, fill = value, label = sprintf(round(value, 2), fmt = '%#.2f'))) +
  geom_tile(colour = "black") + 
  scale_fill_gradient(low = "firebrick2", high = "chartreuse") +
  labs(x = "Método de selección de características",
       y = "Número de genes usados en el modelo",
       fill = "F1-Score",
       title = "F1-Score según método de selección de características\ny número de genes usados en el modelo",
       subtitle = "Resaltado en negrita el mejor F1-Score") +
  scale_y_continuous(breaks = 1:10) + 
  geom_label(fontface = face_etiquetas,
             label.size = 0) + 
  theme_minimal() + 
  theme(plot.title = element_text(hjust = .5),
        panel.grid = element_blank(),
        plot.subtitle = element_text(hjust = .5),
        panel.spacing = unit(c(0, 0, 0, 0), "null"),
        axis.text.x = element_text(size = 11, margin = margin(t = -15)),
        panel.grid.major.x = element_blank(),
        axis.text = element_text(color = "black"))

ggsave(filename = "../../01_images/16_knn_f1_score.png", width = 8, height = 8)

# Gráfico precisión
which(prec == max(prec), arr.ind = TRUE)

face_etiquetas <- rep("plain", 30)
face_etiquetas[c(5,8,11)] <- "bold"

ggplot(melt(t(prec)), aes(Var1, Var2, fill = value, label = sprintf(round(value, 2), fmt = '%#.2f'))) +
  geom_tile(colour = "black") + 
  scale_fill_gradient(low = "firebrick2", high = "chartreuse") +
  labs(x = "Método de selección de características",
       y = "Número de genes usados en el modelo",
       fill = "Precisión",
       title = "Precisión de kNN según método de selección de características\ny número de genes usados en el modelo",
       subtitle = "Resaltado en negrita la mejor precisión") +
  scale_y_continuous(breaks = 1:10) + 
  geom_label(fontface = face_etiquetas,
             label.size = 0) + 
  theme_minimal() + 
  theme(plot.title = element_text(hjust = .5),
        panel.grid = element_blank(),
        plot.subtitle = element_text(hjust = .5),
        panel.spacing = unit(c(0, 0, 0, 0), "null"),
        axis.text.x = element_text(size = 11, margin = margin(t = -15)),
        panel.grid.major.x = element_blank(),
        axis.text = element_text(color = "black"))

ggsave(filename = "../../01_images/17_knn_accuracy.png", width = 8, height = 8)

# --- Mejor método en CV basado en F1-score

# El mejor ranking es RF con 20 genes
mejores_genes_knn <- rfRanking[1:2]

# ----- kNN: Resultados gráficos en train con el mejor método ------

# Se grafican boxplots y mapas de calor con el mejor método
png(filename = "../../01_images/18_heatmap_knn_mejor_metodo.png", width = 13, height = 8, units = "in", res = 300)
DEGsMatrix_heatmap <- DEGsMatrix[mejores_genes_knn, c(which(labels == "Primary Tumor"),
                                                      which(labels == "Solid Tissue Normal"))]
labels_heatmap <- labels[c(which(labels == "Primary Tumor"),
                           which(labels == "Solid Tissue Normal"))]
dataPlot(DEGsMatrix_heatmap[mejores_genes_knn, ], labels_heatmap, mode = "heatmap")
dev.off()

png(filename = "../../01_images/19_boxplots_knn_mejor_metodo.png", width = 13, height = 8, units = "in", res = 300)
dataPlot(DEGsMatrix[mejores_genes_knn, ], labels, mode = "genesBoxplot")
dev.off()

# ----- kNN: Resultados en test del mejor modelo -----

# mRMR
results_knn_mrmr <- knn_test(train = particion.entrenamiento, labels_train,
                             test = particion.test, labels_test, mrmrRanking,
                             bestK = results_cv_knn_mrmr$bestK)

# random forest
results_knn_rf <- knn_test(train = particion.entrenamiento, labels_train,
                           test = particion.test, labels_test, rfRanking,
                           bestK = results_cv_knn_rf$bestK)

# disease association
results_knn_da <- knn_test(train = particion.entrenamiento, labels_train,
                           test = particion.test, labels_test, daRanking,
                           bestK = results_cv_knn_da$bestK)

# Matriz de confusión mejor clasificador
results_knn_rf$cfMats[[2]]$table

# Gráficamente
png(filename = "../../01_images/20_matriz_confusion_knn_mejor_metodo.png", width = 13, height = 8, units = "in", res = 300)
dataPlot(results_knn_rf$cfMats[[2]]$table, labels_test,
         mode = "confusionMatrix")
dev.off()


# ----- SVM: Descarga de información sobre enfermedades relacionadas con DEGs -----


# Se usan los genes más relevantes usando el mejor método encontrado anteriormente
# Se recuperan todas las enfermedades más importantes vinculadas a cada gen

tic("DEGsToDiseases:") # 90 segundos
diseases <- DEGsToDiseases(mejores_genes_svm, getEvidences = TRUE, size = 1000)
toc()

length(diseases)

# Extracción de todas las enfermedades relacionadas con los 10 genes
enfermedades <- c()
for(i in 1:length(diseases)){
  # Se extraen las enfermedades relacionadas
  enfermedades <- c(enfermedades, diseases[[i]]$summary[, 1])
  # Se eliminan duplicados
  enfermedades <- unique(enfermedades)
}

# Todas las enfermedades
enfermedades

# Buscar cáncer de hígado y enfermedades relacionadas
unique(c(enfermedades[grep(pattern = "cancer", ignore.case = T, x = enfermedades)],
         enfermedades[grep(pattern = "neopl",  ignore.case = T, x = enfermedades)],
         enfermedades[grep(pattern = "liver",  ignore.case = T, x = enfermedades)],
         enfermedades[grep(pattern = "alcoh",  ignore.case = T, x = enfermedades)]
         )
       )

# ----- kNN: Descarga de información sobre enfermedades relacionadas con DEGs -----

# Se usan los genes más relevantes usando el mejor método encontrado anteriormente
# Se recuperan todas las enfermedades más importantes vinculadas a cada gen

tic("DEGsToDiseases:") # 90 segundos
diseases <- DEGsToDiseases(mejores_genes_knn, getEvidences = TRUE, size = 1000)
toc()

length(diseases)

# Extracción de todas las enfermedades relacionadas con los 10 genes
enfermedades <- c()
for(i in 1:length(diseases)){
  # Se extraen las enfermedades relacionadas
  enfermedades <- c(enfermedades, diseases[[i]]$summary[, 1])
  # Se eliminan duplicados
  enfermedades <- unique(enfermedades)
}

# Todas las enfermedades
enfermedades

# Buscar cáncer de hígado y enfermedades relacionadas
unique(c(enfermedades[grep(pattern = "cancer", ignore.case = T, x = enfermedades)],
         enfermedades[grep(pattern = "neopl",  ignore.case = T, x = enfermedades)],
         enfermedades[grep(pattern = "liver",  ignore.case = T, x = enfermedades)],
         enfermedades[grep(pattern = "alcoh",  ignore.case = T, x = enfermedades)]
         )
       )

# ----- Guardar imagen para Shiny -----
matriz <- rbind(labels, DEGsMatrix)
save(matriz, file =  "../saved_files/matriz.RData")
toc()


# Guardar imagen
#save.image(file =  "../saved_files/01_workspace_completo.RData")
#load("../saved_files/01_workspace_completo.RData")
