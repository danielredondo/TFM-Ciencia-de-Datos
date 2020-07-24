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
# Si hay problemas del tipo "Error: (converted from warning)" : Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")library(KnowSeq)     # Para trabajar con datos de transcriptómica de GDC Portal
library(KnowSeq)     # Para trabajar con datos de transcriptómica de GDC Portal
library(dplyr)       # Para select, filter, pipes, ...
library(tictoc)      # Para medir tiempos con tic() y toc() a lo MATLAB 
library(beepr)       # Para avisar con beeps cuando acaba un proceso
library(caret)       # Para machine learning 
library(e1071)       # Para SVM
library(gplots)      # Para heatmaps
library(reshape2)    # Para melt
library(ggalluvial)  # Para diagrama de Sankey
library(readr)       # Para leer el fichero de datos clínicos

# -----  Sobreescribir la función dataPlot con la nueva función que pinta líneas discontinuas -----

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

# Añadir morfología del tumor: leer datos clínicos
clinical <- read_delim("../data/clinical.cart.2020-06-11/clinical.tsv", "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
  select(case_submitter_id, primary_diagnosis) %>% 
  unique

# Anexar morfología
tumores <- samplesInfo %>% filter(Sample.Type == "Primary Tumor")
nrow(tumores) == nrow(clinical) # Comprobación: todos los tumores tienen datos clínicos

tumores <- left_join(x = tumores, y = clinical, by = c("Case.ID" = "case_submitter_id"))

samplesInfo <- samplesInfo %>% filter(Sample.Type != "Primary Tumor") %>%
  mutate(primary_diagnosis = "Solid Tissue Normal") %>%
  rbind(tumores)

# Recodificación morfologías
table(samplesInfo$primary_diagnosis)
samplesInfo[samplesInfo$primary_diagnosis %in% c("Combined hepatocellular carcinoma and cholangiocarcinoma",
                                           "Clear cell adenocarcinoma, NOS"), "primary_diagnosis"] <- "otros"
samplesInfo[samplesInfo$primary_diagnosis == "Solid Tissue Normal", "primary_diagnosis"] <- "tejido_normal"
samplesInfo[samplesInfo$primary_diagnosis %in% c("Hepatocellular carcinoma, NOS", "Hepatocellular carcinoma, clear cell type",
                                                 "Hepatocellular carcinoma, fibrolamellar", "Hepatocellular carcinoma, spindle cell variant"),
            "primary_diagnosis"] <- "carc_hepatocelular"
samplesInfo[samplesInfo$primary_diagnosis == "Cholangiocarcinoma", "primary_diagnosis"] <- "colangiocarcinoma"

table(samplesInfo$primary_diagnosis)

# Quitar morfología de clases minoritarias
samplesInfo <- samplesInfo %>% filter(primary_diagnosis != "otros")

# Definición de parámetros
Run <- samplesInfo$File.Name
Path <- samplesInfo$File.ID
Class <- samplesInfo$primary_diagnosis

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
write.csv(SamplesDataFrame, file = "04_SamplesDataFrame.csv")

# ----- Unificación de ficheros en formato matriz -----

tic("countsToMatrix") # 85 segundos
countsInformation <- countsToMatrix("04_SamplesDataFrame.csv", extension = "counts")
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

png(filename = "../../04_analisis_multiclase/01_boxplot_all.png", width = 13, height = 8, units = "in", res = 300)
dataPlot(expressionMatrix, labels, mode = "boxplot", colours = c("red", "orange", "green"))
legend(1, 29, c("colangiocarcinoma", "carc_hepatocelular", "tejido_normal"),
       fill = c("red", "orange", "green"), xpd = TRUE)
title("Boxplot de matriz de expresión para todos los genes")
dev.off()

png(filename = "../../04_analisis_multiclase/02_boxplot_all_ordered.png", width = 13, height = 8, units = "in", res = 300)
dataPlot(expressionMatrix, labels, mode = "orderedBoxplot", colours = c("red", "orange", "green"))
legend(1, 29, c("colangiocarcinoma", "carc_hepatocelular", "tejido_normal"),
       fill = c("red", "orange", "green"), xpd = TRUE)
title("Boxplot de matriz de expresión para todos los genes")
dev.off()

# ----- Controlando por el efecto batch -----

tic("batchEffectRemoval") # 207 segundos
svaMod <- batchEffectRemoval(expressionMatrix, as.factor(labels), method = "sva")
toc()

# ----- Extracción de DEG (Expresión Diferencial de Genes) -----

tic("DEGsExtraction") # 1 segundo
DEGsInformation <- DEGsExtraction(expressionMatrix, as.factor(labels),
                                  # p-valor
                                  pvalue = 0.001,
                                  # Número de genes
                                  #number = 12,
                                  # Ajuste por batchEffect
                                  svaCorrection = TRUE, svaMod = svaMod)
toc()

# Número de genes extraídos: 8533
print(nrow(DEGsInformation$DEGsMatrix))

topTable <- DEGsInformation$Table
DEGsMatrix <- DEGsInformation$DEGsMatrix

# Se traspone la matriz
DEGsMatrixML <- t(DEGsMatrix)

# ----- Representación de DEG -----

# Boxplots para todas las muestras de los primeros genes
png(filename = "../../04_analisis_multiclase/03_boxplot_primeros_genes.png", width = 13, height = 8, units = "in", res = 300)
dataPlot(DEGsMatrix[1:numero_de_genes, ], labels, mode = "genesBoxplot", colours = c("green", "orange", "red"), columnasGenesBoxplot = 2)
dev.off()

# Mapa de calor para todas las muestras de los primeros genes
png(filename = "../../04_analisis_multiclase/04_heatmap_primeros_genes.png", width = 13, height = 8, units = "in", res = 300)
DEGsMatrix_heatmap <- DEGsMatrix[, c(which(labels == "colangiocarcinoma"),
                                     which(labels == "carc_hepatocelular"),
                                     which(labels == "tejido_normal"))]
labels_heatmap <- labels[c(which(labels == "colangiocarcinoma"),
                           which(labels == "carc_hepatocelular"),
                           which(labels == "tejido_normal"))]
dataPlot(DEGsMatrix_heatmap[1:numero_de_genes, ], labels_heatmap, mode = "heatmap", colours = c("red", "orange", "green"))
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
entr_car <- table(labels_train)[1]
entr_col <- table(labels_train)[2]
entr_nor <- table(labels_train)[3]

# Test
table(labels_test)
test_car <- table(labels_test)[1]
test_col <- table(labels_test)[2]
test_nor <- table(labels_test)[3]

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
datos_sankey <- data.frame(tipo = c(paste0("Carcinoma hepatocelular\n", entr_car + test_car, " casos"),
                                    paste0("Carcinoma hepatocelular\n", entr_car + test_car, " casos"),
                                    paste0("Colangiocarcinoma\n", entr_col + test_col, " casos"),
                                    paste0("Colangiocarcinoma\n", entr_col + test_col, " casos"),
                                    paste0("Tejido normal\n", entr_nor + test_nor, " casos"),
                                    paste0("Tejido normal\n", entr_nor + test_nor, " casos")),
                           traintest = c(paste0("Entrenamiento\n", entr_col, " colangiocarcinoma\n", entr_car, " carcinoma hepatocelular\n", entr_nor, " tejido normal"),
                                         paste0("Test\n", test_col, " colangiocarcinoma\n", test_car, " carcinoma hepatocelular\n", test_nor, " tejido normal"),
                                         paste0("Entrenamiento\n", entr_col, " colangiocarcinoma\n", entr_car, " carcinoma hepatocelular\n", entr_nor, " tejido normal"),
                                         paste0("Test\n", test_col, " colangiocarcinoma\n", test_car, " carcinoma hepatocelular\n", test_nor, " tejido normal"),
                                         paste0("Entrenamiento\n", entr_col, " colangiocarcinoma\n", entr_car, " carcinoma hepatocelular\n", entr_nor, " tejido normal"),
                                         paste0("Test\n", test_col, " colangiocarcinoma\n", test_car, " carcinoma hepatocelular\n", test_nor, " tejido normal")),
                           value = c(entr_car, test_car, entr_col, test_col, entr_nor, test_nor))

# Pequeño reorden para que mejorar la presentación de los datos
datos_sankey$tipo <- factor(datos_sankey$tipo,
                            levels = c(paste0("Colangiocarcinoma\n", entr_col + test_col, " casos"),
                                       paste0("Carcinoma hepatocelular\n", entr_car + test_car, " casos"),
                                       paste0("Tejido normal\n", entr_nor + test_nor, " casos")),
                            ordered = T)

ggplot(data = datos_sankey,
       aes(axis1 = tipo, axis2 = traintest, y = value)) +
  scale_x_discrete(limits = c("Tipo de muestra", "Partición\nentrenamiento-test"),
                   expand = c(.1, .05)) +
  ylab("") +
  geom_alluvium(col = "black", alpha = 1) +
  geom_alluvium(fill = rep(c("red", "red", "orange", "orange",
                             "green", "green"), 2), alpha = .6, show.legend = FALSE) +
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

ggsave(filename = "../../04_analisis_multiclase/05_sankey.png", width = 8, height = 7)

# ----- Selección de características: numero_de_genes genes más importantes -----

# Método mRMR (mínima redundancia, máxima relevancia)
mrmrRanking <- featureSelection(particion.entrenamiento, labels_train, colnames(particion.entrenamiento),
                                mode = "mrmr")
mrmrRanking <- names(mrmrRanking)[1:numero_de_genes]

# Método random forest
rfRanking <- featureSelection(particion.entrenamiento, labels_train, colnames(particion.entrenamiento),
                              mode = "rf")
rfRanking <- rfRanking[1:numero_de_genes]

# Método Disease association ranking (en base a scores obtenidos en la literatura)
daRanking <- NULL

daRanking <- featureSelection(particion.entrenamiento, labels_train, colnames(particion.entrenamiento),
                              mode = "da", disease = "liver cancer")

# Si ha habido algún problema en la llamada a la API, se repite tras un descanso de 5 segundos
while(is.null(daRanking)){
  Sys.sleep(5)
  daRanking <- featureSelection(particion.entrenamiento, labels_train, colnames(particion.entrenamiento),
                                mode = "da", disease = "liver cancer")
}
daRanking <- names(daRanking[1:numero_de_genes])

genes <- cbind(mrmrRanking, rfRanking, daRanking)
write.csv2(genes, file = "../../04_analisis_multiclase/genes.csv")

# ----- SVM: Resultados con validación cruzada para cada método de selección de características -----

numero_folds <- 5

# mRMR: 24 segundos
tic("svm_mrmr")
results_cv_svm_mrmr <- svm_trn(particion.entrenamiento, labels_train, mrmrRanking,
                               numFold = numero_folds)
toc()

# random forest: 22 segundos
tic("svm_rf")
results_cv_svm_rf <- svm_trn(particion.entrenamiento, labels_train, rfRanking,
                             numFold = numero_folds)
toc()

# disease association: 40 segundos
tic("svm_da")
results_cv_svm_da <- svm_trn(particion.entrenamiento, labels_train, daRanking,
                             numFold = numero_folds)
toc()

# Tabla con mejores parámetros
mejores_parametros_svm <- rbind("mrmr" = results_cv_svm_mrmr$bestParameters,
                                "rf" = results_cv_svm_rf$bestParameters,
                                "da" = results_cv_svm_da$bestParameters)

write.csv2(mejores_parametros_svm, file = "../../04_analisis_multiclase/mejores_parametros_svm.csv")

# ----- SVM: Resultados gráficos de validación cruzada -----

# Plotting the accuracy of all the folds evaluated in the CV process
png(filename = "../../04_analisis_multiclase/06_svm_MRMR_f1.png", width = 13, height = 6, units = "in", res = 300)
dataPlot(results_cv_svm_mrmr$F1Info$meanF1, colours = c("black", "darkblue"),
         mode = "classResults",
         main = "SVM + mRMR - Accuracy media", xlab = "Genes", ylab = "F1-Score")
dev.off()

# Plotting the accuracy of all the folds evaluated in the CV process
png(filename = "../../04_analisis_multiclase/07_svm_RF_f1.png", width = 13, height = 6, units = "in", res = 300)
dataPlot(results_cv_svm_rf$F1Info$meanF1, colours = c("black", "darkblue"),
         mode = "classResults",
         main = "SVM + RF - Accuracy media", xlab = "Genes", ylab = "F1-Score")
dev.off()

# Plotting the accuracy of all the folds evaluated in the CV process
png(filename = "../../04_analisis_multiclase/08_svm_DA_f1.png", width = 13, height = 6, units = "in", res = 300)
dataPlot(results_cv_svm_da$F1Info$meanF1, colours = c("black", "darkblue"),
         mode = "classResults",
         main = "SVM + DA - Accuracy media", xlab = "Genes", ylab = "F1-Score")
dev.off()

# ----- SVM: Mejor método en CV basado en precisión, especificidad, sensibilidad y F1-score -----
genes_a_usar <- c(1:numero_de_genes)

# Precisión
prec <- matrix(0, nrow = length(genes_a_usar), ncol = 3)
l_prec <- matrix("0", nrow = length(genes_a_usar), ncol = 3)
colnames(prec) <- c("MRMR", "RF", "DA")

for(i in 1:length(genes_a_usar)){
  prec[i, ] <- 100 * c(results_cv_svm_mrmr$accuracyInfo$meanAccuracy[genes_a_usar[i]],
                       results_cv_svm_rf$accuracyInfo$meanAccuracy[genes_a_usar[i]],
                       results_cv_svm_da$accuracyInfo$meanAccuracy[genes_a_usar[i]])
  l_prec[i, ] <- paste0(sprintf(prec[i, ], fmt = '%#.2f'), 
                        " (",
                        sprintf(100 * c(results_cv_svm_mrmr$accuracyInfo$standardDeviation[i],
                                        results_cv_svm_rf$accuracyInfo$standardDeviation[i],
                                        results_cv_svm_da$accuracyInfo$standardDeviation[i]
                        ), fmt = '%#.1f'), 
                        ")")
}

#F1-Score
f1 <- matrix(0, nrow = length(genes_a_usar), ncol = 3)
l_f1 <- matrix("0", nrow = length(genes_a_usar), ncol = 3)
colnames(f1) <- c("MRMR", "RF", "DA")

for(i in 1:length(genes_a_usar)){
  f1[i, ] <- 100 * c(results_cv_svm_mrmr$F1Info$meanF1[genes_a_usar[i]],
                     results_cv_svm_rf$F1Info$meanF1[genes_a_usar[i]],
                     results_cv_svm_da$F1Info$meanF1[genes_a_usar[i]])
  l_f1[i, ] <- paste0(sprintf(f1[i, ], fmt = '%#.2f'), 
                      " (",
                      sprintf(100 * c(results_cv_svm_mrmr$F1Info$standardDeviation[i],
                                      results_cv_svm_rf$F1Info$standardDeviation[i],
                                      results_cv_svm_da$F1Info$standardDeviation[i]
                      ), fmt = '%#.1f'), 
                      ")")
}

which(f1 == max(f1), arr.ind = TRUE)
# El mejor F1-score se obtiene con 7 genes en mRMR

# Tablas y gráficos precisión y F1-score
mejores <- which(f1 == max(f1), arr.ind = TRUE)
face_etiquetas <- matrix("plain", nrow = length(genes_a_usar), ncol = 3)
for(i in 1:nrow(mejores)) face_etiquetas[mejores[i, 1], mejores[i, 2]] <- "bold"
face_etiquetas <- as.vector(t(face_etiquetas))

data <- cbind(melt(t(f1)), melt(t(l_f1))$value)
names(data)[4] <- "l_f1"

ggplot(data, aes(Var1, Var2, fill = value, label = l_f1)) +
  geom_tile(colour = "black") + 
  scale_fill_gradient(low = "firebrick2", high = "chartreuse") +
  labs(x = "Método de selección de características",
       y = "Número de genes usados en el modelo",
       fill = "F1-Score",
       title = "F1-Score de SVM según método de selección de características\ny número de genes usados en el modelo",
       subtitle = "F1-Score medio en los 5 fold (desviación típica)") +
  scale_y_continuous(breaks = 1:length(genes_a_usar)) + 
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

ggsave(filename = "../../04_analisis_multiclase/09_svm_f1_score.png", width = 8, height = 8)

# Gráfico precisión
mejores <- which(prec == max(prec), arr.ind = TRUE)
face_etiquetas <- matrix("plain", nrow = length(genes_a_usar), ncol = 3)
for(i in 1:nrow(mejores)) face_etiquetas[mejores[i, 1], mejores[i, 2]] <- "bold"
face_etiquetas <- as.vector(t(face_etiquetas))

data <- cbind(melt(t(prec)), melt(t(l_prec))$value)
names(data)[4] <- "l_prec"

ggplot(data, aes(Var1, Var2, fill = value, label = l_prec)) +
  geom_tile(colour = "black") + 
  scale_fill_gradient(low = "firebrick2", high = "chartreuse") +
  labs(x = "Método de selección de características",
       y = "Número de genes usados en el modelo",
       fill = "Precisión",
       title = "Precisión de SVM según método de selección de características\ny número de genes usados en el modelo",
       subtitle = "Precisión media en los 5 fold (desviación típica)") +
  scale_y_continuous(breaks = 1:length(genes_a_usar)) + 
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

ggsave(filename = "../../04_analisis_multiclase/10_svm_accuracy.png", width = 8, height = 8)

# ----- SVM: Resultados gráficos en train con el mejor método ------

mejores_genes_svm <- mrmrRanking[1:7]

# Se grafican boxplots y mapas de calor con el mejor método
png(filename = "../../04_analisis_multiclase/11_svm_heatmap_mejor_metodo.png", width = 13, height = 8, units = "in", res = 300)
DEGsMatrix_heatmap <- DEGsMatrix[mejores_genes_svm, c(which(labels == "colangiocarcinoma"),
                                                      which(labels == "carc_hepatocelular"),
                                                      which(labels == "tejido_normal"))]
labels_heatmap <- labels[c(which(labels == "colangiocarcinoma"),
                           which(labels == "carc_hepatocelular"),
                           which(labels == "tejido_normal"))]
dataPlot(DEGsMatrix_heatmap[mejores_genes_svm, ], labels_heatmap, mode = "heatmap", colours = c("red", "orange", "green"))
dev.off()

png(filename = "../../04_analisis_multiclase/12_svm_boxplots_mejor_metodo.png", width = 13, height = 8, units = "in", res = 300)
dataPlot(DEGsMatrix[mejores_genes_svm, ], labels, mode = "genesBoxplot", colours = c("green", "red", "orange"))
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
tabla <- results_svm_mrmr$cfMats[[7]]$table

# Gráficamente
png(filename = "../../04_analisis_multiclase/13_svm_matriz_confusion_mejor_metodo.png", width = 13, height = 8, units = "in", res = 300)
dataPlot(tabla, labels_test,
         mode = "confusionMatrix")
dev.off()

# ----- RF: Resultados con validación cruzada para cada método de selección de características -----
numero_folds <- 5

# mRMR
results_cv_rf_mrmr <- rf_trn(particion.entrenamiento, labels_train, mrmrRanking,
                             numFold = numero_folds)

# random forest
results_cv_rf_rf <- rf_trn(particion.entrenamiento, labels_train, rfRanking,
                           numFold = numero_folds)

# disease association
results_cv_rf_da <- rf_trn(particion.entrenamiento, labels_train, daRanking,
                           numFold = numero_folds)

# ----- RF: Resultados gráficos de validación cruzada -----

# Plotting the accuracy of all the folds evaluated in the CV process
png(filename = "../../04_analisis_multiclase/14_rf_MRMR_f1.png", width = 13, height = 6, units = "in", res = 300)
dataPlot(results_cv_rf_mrmr$F1Info$meanF1, colours = c("black", "darkblue"),
         mode = "classResults",
         main = "RF + mRMR - F1-Score medio", xlab = "Genes", ylab = "F1-Score")
dev.off()

# Plotting the accuracy of all the folds evaluated in the CV process
png(filename = "../../04_analisis_multiclase/15_rf_RF_f1.png", width = 13, height = 6, units = "in", res = 300)
dataPlot(results_cv_rf_rf$F1Info$meanF1, colours = c("black", "darkblue"),
         mode = "classResults",
         main = "RF + RF - F1-Score medio", xlab = "Genes", ylab = "F1-Score")
dev.off()

# Plotting the accuracy of all the folds evaluated in the CV process
png(filename = "../../04_analisis_multiclase/16_rf_DA_f1.png", width = 13, height = 6, units = "in", res = 300)
dataPlot(results_cv_rf_da$F1Info$meanF1, colours = c("black", "darkblue"),
         mode = "classResults",
         main = "RF + DA - F1-Score medio", xlab = "Genes", ylab = "F1-Score")
dev.off()

# ----- RF: Mejor método en CV basado en precisión, especificidad, sensibilidad y F1-score -----
genes_a_usar <- c(1:numero_de_genes)

# Precisión
prec <- matrix(0, nrow = length(genes_a_usar), ncol = 3)
l_prec <- matrix("0", nrow = length(genes_a_usar), ncol = 3)
colnames(prec) <- c("MRMR", "RF", "DA")

for(i in 1:length(genes_a_usar)){
  prec[i, ] <- 100 * c(results_cv_rf_mrmr$accuracyInfo$meanAccuracy[genes_a_usar[i]],
                       results_cv_rf_rf$accuracyInfo$meanAccuracy[genes_a_usar[i]],
                       results_cv_rf_da$accuracyInfo$meanAccuracy[genes_a_usar[i]])
  l_prec[i, ] <- paste0(sprintf(prec[i, ], fmt = '%#.2f'), 
                        " (",
                        sprintf(100 * c(results_cv_rf_mrmr$accuracyInfo$standardDeviation[i],
                                        results_cv_rf_rf$accuracyInfo$standardDeviation[i],
                                        results_cv_rf_da$accuracyInfo$standardDeviation[i]
                        ), fmt = '%#.1f'), 
                        ")")
}

which(prec == max(prec), arr.ind = TRUE)

# F1-Score
f1 <- matrix(0, nrow = length(genes_a_usar), ncol = 3)
l_f1 <- matrix("0", nrow = length(genes_a_usar), ncol = 3)
colnames(f1) <- c("MRMR", "RF", "DA")

for(i in 1:length(genes_a_usar)){
  f1[i, ] <- 100 * c(results_cv_rf_mrmr$F1Info$meanF1[genes_a_usar[i]],
                       results_cv_rf_rf$F1Info$meanF1[genes_a_usar[i]],
                       results_cv_rf_da$F1Info$meanF1[genes_a_usar[i]])
  l_f1[i, ] <- paste0(sprintf(f1[i, ], fmt = '%#.2f'), 
                        " (",
                        sprintf(100 * c(results_cv_rf_mrmr$F1Info$standardDeviation[i],
                                        results_cv_rf_rf$F1Info$standardDeviation[i],
                                        results_cv_rf_da$F1Info$standardDeviation[i]
                        ), fmt = '%#.1f'), 
                        ")")
}

which(f1 == max(f1), arr.ind = TRUE)
# El mejor F1-score se obtiene con 7 genes en mrmr

# Tablas y gráficos precisión, especificidad, sensibilidad y F1-score
which(prec == max(prec), arr.ind = TRUE)
which(f1 == max(f1), arr.ind = TRUE)

mejores <- which(f1 == max(f1), arr.ind = TRUE)
face_etiquetas <- matrix("plain", nrow = length(genes_a_usar), ncol = 3)
for(i in 1:nrow(mejores)) face_etiquetas[mejores[i, 1], mejores[i, 2]] <- "bold"
face_etiquetas <- as.vector(t(face_etiquetas))

data <- cbind(melt(t(f1)), melt(t(l_f1))$value)
names(data)[4] <- "l_f1"

ggplot(data, aes(Var1, Var2, fill = value, label = l_f1)) +
  geom_tile(colour = "black") + 
  scale_fill_gradient(low = "firebrick2", high = "chartreuse") +
  labs(x = "Método de selección de características",
       y = "Número de genes usados en el modelo",
       fill = "F1-Score",
       title = "F1-Score según método de selección de características\ny número de genes usados en el modelo",
       subtitle = "F1-Score medio en los 5 fold (desviación típica)") +
  scale_y_continuous(breaks = 1:length(genes_a_usar)) + 
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

ggsave(filename = "../../04_analisis_multiclase/17_rf_f1_score.png", width = 8, height = 8)

# Gráfico precisión
which(prec == max(prec), arr.ind = TRUE)

mejores <- which(prec == max(prec), arr.ind = TRUE)
face_etiquetas <- matrix("plain", nrow = length(genes_a_usar), ncol = 3)
for(i in 1:nrow(mejores)) face_etiquetas[mejores[i, 1], mejores[i, 2]] <- "bold"
face_etiquetas <- as.vector(t(face_etiquetas))

data <- cbind(melt(t(prec)), melt(t(l_prec))$value)
names(data)[4] <- "l_prec"

ggplot(data, aes(Var1, Var2, fill = value, label = l_prec)) +
  geom_tile(colour = "black") + 
  scale_fill_gradient(low = "firebrick2", high = "chartreuse") +
  labs(x = "Método de selección de características",
       y = "Número de genes usados en el modelo",
       fill = "Precisión",
       title = "Precisión de RF según método de selección de características\ny número de genes usados en el modelo",
       subtitle = "Precisión media en los 5 fold (desviación típica)") +
  scale_y_continuous(breaks = 1:length(genes_a_usar)) + 
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

ggsave(filename = "../../04_analisis_multiclase/18_rf_accuracy.png", width = 8, height = 8)

# El mejor ranking es mRMR con 6 genes
mejores_genes_rf <- mrmrRanking[1:6]

# ----- RF: Resultados gráficos en train con el mejor método ------

# Se grafican boxplots y mapas de calor con el mejor método
png(filename = "../../04_analisis_multiclase/19_rf_heatmap_mejor_metodo.png", width = 13, height = 8, units = "in", res = 300)
DEGsMatrix_heatmap <- DEGsMatrix[mejores_genes_rf, c(which(labels == "colangiocarcinoma"),
                                                     which(labels == "carc_hepatocelular"),
                                                     which(labels == "tejido_normal"))]
labels_heatmap <- labels[c(which(labels == "colangiocarcinoma"),
                           which(labels == "carc_hepatocelular"),
                           which(labels == "tejido_normal"))]
dataPlot(DEGsMatrix_heatmap[mejores_genes_rf, ], labels_heatmap, mode = "heatmap", colours = c("red", "orange", "green"))
dev.off()

png(filename = "../../04_analisis_multiclase/20_rf_boxplots_mejor_metodo.png", width = 13, height = 8, units = "in", res = 300)
dataPlot(DEGsMatrix[mejores_genes_rf, ], labels, mode = "genesBoxplot", colours = c("green", "red", "orange"))
dev.off()

# ----- RF: Resultados en test del mejor modelo -----

# mRMR
results_rf_mrmr <- rf_test(train = particion.entrenamiento, labels_train,
                           test = particion.test, labels_test, mrmrRanking)

# random forest
results_rf_rf <- rf_test(train = particion.entrenamiento, labels_train,
                         test = particion.test, labels_test, rfRanking)

# disease association
results_rf_da <- rf_test(train = particion.entrenamiento, labels_train,
                         test = particion.test, labels_test, daRanking)

# Matriz de confusión mejor clasificador
tabla <- results_rf_mrmr$cfMats[[6]]$table

# Gráficamente
png(filename = "../../04_analisis_multiclase/21_rf_matriz_confusion_mejor_metodo.png", width = 13, height = 8, units = "in", res = 300)
dataPlot(tabla, labels_test,
         mode = "confusionMatrix")
dev.off()

# ----- kNN: Resultados con validación cruzada para cada método de selección de características -----
numero_folds <- 5

# mRMR
results_cv_knn_mrmr <- knn_trn(particion.entrenamiento, labels_train, mrmrRanking,
                               numFold = numero_folds)
results_cv_knn_mrmr$bestK

# random forest
results_cv_knn_rf <- knn_trn(particion.entrenamiento, labels_train, rfRanking,
                             numFold = numero_folds)
results_cv_knn_rf$bestK

# disease association
results_cv_knn_da <- knn_trn(particion.entrenamiento, labels_train, daRanking,
                             numFold = numero_folds)
results_cv_knn_da$bestK

mejores_parametros_knn <- cbind(results_cv_knn_mrmr$bestK, results_cv_knn_rf$bestK, results_cv_knn_da$bestK)
write.csv2(mejores_parametros_knn, file = "../../04_analisis_multiclase/mejores_parametros_knn.csv")

# ----- kNN: Resultados gráficos de validación cruzada -----

# Plotting the accuracy of all the folds evaluated in the CV process
png(filename = "../../04_analisis_multiclase/22_knn_MRMR_f1.png", width = 13, height = 6, units = "in", res = 300)
dataPlot(results_cv_knn_mrmr$F1Info$meanF1, colours = c("black", "darkblue"),
         mode = "classResults",
         main = "kNN + mRMR - F1-Score medio", xlab = "Genes", ylab = "F1-Score")
dev.off()

# Plotting the accuracy of all the folds evaluated in the CV process
png(filename = "../../04_analisis_multiclase/23_knn_RF_f1.png", width = 13, height = 6, units = "in", res = 300)
dataPlot(results_cv_knn_rf$F1Info$meanF1,  colours = c("black", "darkblue"),
         mode = "classResults",
         main = "kNN + RF - F1-Score medio", xlab = "Genes", ylab = "F1-Score")
dev.off()

# Plotting the accuracy of all the folds evaluated in the CV process
png(filename = "../../04_analisis_multiclase/24_knn_DA_f1.png", width = 13, height = 6, units = "in", res = 300)
dataPlot(results_cv_knn_da$F1Info$meanF1,  colours = c("black", "darkblue"),
         mode = "classResults",
         main = "kNN + DA - F1-Score medio", xlab = "Genes", ylab = "F1-Score")
dev.off()

# ----- kNN: Mejor método en CV basado en precisión, especificidad, sensibilidad y F1-score -----
genes_a_usar <- c(1:numero_de_genes)

# Precisión
prec <- matrix(0, nrow = length(genes_a_usar), ncol = 3)
l_prec <- matrix("0", nrow = length(genes_a_usar), ncol = 3)
colnames(prec) <- c("MRMR", "RF", "DA")

for(i in 1:length(genes_a_usar)){
  prec[i, ] <- 100 * c(results_cv_knn_mrmr$accuracyInfo$meanAccuracy[genes_a_usar[i]],
                       results_cv_knn_rf$accuracyInfo$meanAccuracy[genes_a_usar[i]],
                       results_cv_knn_da$accuracyInfo$meanAccuracy[genes_a_usar[i]])
  l_prec[i, ] <- paste0(sprintf(prec[i, ], fmt = '%#.2f'), 
                        " (",
                        sprintf(100 * c(results_cv_knn_mrmr$accuracyInfo$standardDeviation[i],
                                        results_cv_knn_rf$accuracyInfo$standardDeviation[i],
                                        results_cv_knn_da$accuracyInfo$standardDeviation[i]
                        ), fmt = '%#.1f'), 
                        ")")
}

which(prec == max(prec), arr.ind = TRUE)

# F1-Score
f1 <- matrix(0, nrow = length(genes_a_usar), ncol = 3)
l_f1 <- matrix("0", nrow = length(genes_a_usar), ncol = 3)
colnames(f1) <- c("MRMR", "RF", "DA")

for(i in 1:length(genes_a_usar)){
  f1[i, ] <- 100 * c(results_cv_knn_mrmr$F1Info$meanF1[genes_a_usar[i]],
                     results_cv_knn_rf$F1Info$meanF1[genes_a_usar[i]],
                     results_cv_knn_da$F1Info$meanF1[genes_a_usar[i]])
  l_f1[i, ] <- paste0(sprintf(f1[i, ], fmt = '%#.2f'), 
                      " (",
                      sprintf(100 * c(results_cv_knn_mrmr$F1Info$standardDeviation[i],
                                      results_cv_knn_rf$F1Info$standardDeviation[i],
                                      results_cv_knn_da$F1Info$standardDeviation[i]
                      ), fmt = '%#.1f'), 
                      ")")
}

which(f1 == max(f1), arr.ind = TRUE)
# El mejor F1-score se obtiene con 1 gen en mrmr

# Tablas y gráficos precisión, especificidad, sensibilidad y F1-score
which(prec == max(prec), arr.ind = TRUE)
which(f1 == max(f1), arr.ind = TRUE)

mejores <- which(f1 == max(f1), arr.ind = TRUE)
face_etiquetas <- matrix("plain", nrow = length(genes_a_usar), ncol = 3)
for(i in 1:nrow(mejores)) face_etiquetas[mejores[i, 1], mejores[i, 2]] <- "bold"
face_etiquetas <- as.vector(t(face_etiquetas))

data <- cbind(melt(t(f1)), melt(t(l_f1))$value)
names(data)[4] <- "l_f1"

ggplot(data, aes(Var1, Var2, fill = value, label = l_f1)) +
  geom_tile(colour = "black") + 
  scale_fill_gradient(low = "firebrick2", high = "chartreuse") +
  labs(x = "Método de selección de características",
       y = "Número de genes usados en el modelo",
       fill = "F1-Score",
       title = "F1-Score según método de selección de características\ny número de genes usados en el modelo",
       subtitle = "F1-Score medio en los 5 fold (desviación típica)") +
  scale_y_continuous(breaks = 1:length(genes_a_usar)) + 
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

ggsave(filename = "../../04_analisis_multiclase/25_knn_f1_score.png", width = 8, height = 8)

# Gráfico precisión
which(prec == max(prec), arr.ind = TRUE)

mejores <- which(prec == max(prec), arr.ind = TRUE)
face_etiquetas <- matrix("plain", nrow = length(genes_a_usar), ncol = 3)
for(i in 1:nrow(mejores)) face_etiquetas[mejores[i, 1], mejores[i, 2]] <- "bold"
face_etiquetas <- as.vector(t(face_etiquetas))

data <- cbind(melt(t(prec)), melt(t(l_prec))$value)
names(data)[4] <- "l_prec"

ggplot(data, aes(Var1, Var2, fill = value, label = l_prec)) +
  geom_tile(colour = "black") + 
  scale_fill_gradient(low = "firebrick2", high = "chartreuse") +
  labs(x = "Método de selección de características",
       y = "Número de genes usados en el modelo",
       fill = "Precisión",
       title = "Precisión de kNN según método de selección de características\ny número de genes usados en el modelo",
       subtitle = "Precisión media en los 5 fold (desviación típica)") +
  scale_y_continuous(breaks = 1:length(genes_a_usar)) + 
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

ggsave(filename = "../../04_analisis_multiclase/26_knn_accuracy.png", width = 8, height = 8)

# --- Mejor método en CV basado en F1-score

# El mejor ranking es mrmr con 7 genes
mejores_genes_knn <- mrmrRanking[1:7]

# ----- kNN: Resultados gráficos en train con el mejor método ------

# Se grafican boxplots y mapas de calor con el mejor método
png(filename = "../../04_analisis_multiclase/27_knn_heatmap_mejor_metodo.png", width = 13, height = 8, units = "in", res = 300)
DEGsMatrix_heatmap <- DEGsMatrix[mejores_genes_knn, c(which(labels == "colangiocarcinoma"),
                                                      which(labels == "carc_hepatocelular"),
                                                      which(labels == "tejido_normal"))]
labels_heatmap <- labels[c(which(labels == "colangiocarcinoma"),
                           which(labels == "carc_hepatocelular"),
                           which(labels == "tejido_normal"))]
dataPlot(DEGsMatrix_heatmap[mejores_genes_knn, ], labels_heatmap, mode = "heatmap", colours = c("red", "orange", "green"))
title("Segundo gen añadido sólo para mostrar el gráfico. Recortar figura")
dev.off()

png(filename = "../../04_analisis_multiclase/28_knn_boxplots_mejor_metodo.png", width = 13, height = 8, units = "in", res = 300)
dataPlot(DEGsMatrix[mejores_genes_knn, ], labels, mode = "genesBoxplot", colours = c("green", "red", "orange"), 
         main = "Segundo gen añadido sólo para mostrar el gráfico. Recortar figura")
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
tabla <- results_knn_mrmr$cfMats[[7]]$table

# Gráficamente
png(filename = "../../04_analisis_multiclase/29_knn_matriz_confusion_mejor_metodo.png", width = 13, height = 8, units = "in", res = 300)
dataPlot(tabla, labels_test, mode = "confusionMatrix")
dev.off()

# ----- SVM: Descarga de información sobre enfermedades relacionadas con DEGs -----

# Se usan los genes más relevantes usando el mejor método encontrado anteriormente
# Se recuperan todas las enfermedades más importantes vinculadas a cada gen

tic("DEGsToDiseases:") # 90 segundos
diseases <- DEGsToDiseases(mejores_genes_svm, getEvidences = TRUE, size = 10)
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
         enfermedades[grep(pattern = "alcoh",  ignore.case = T, x = enfermedades)]))

# ----- RF: Descarga de información sobre enfermedades relacionadas con DEGs -----

# Se usan los genes más relevantes usando el mejor método encontrado anteriormente
# Se recuperan todas las enfermedades más importantes vinculadas a cada gen

tic("DEGsToDiseases:") # 90 segundos
diseases <- DEGsToDiseases(mejores_genes_rf, getEvidences = TRUE, size = 10)
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
         enfermedades[grep(pattern = "alcoh",  ignore.case = T, x = enfermedades)]))

# ----- kNN: Descarga de información sobre enfermedades relacionadas con DEGs -----

# Se usan los genes más relevantes usando el mejor método encontrado anteriormente
# Se recuperan todas las enfermedades más importantes vinculadas a cada gen

tic("DEGsToDiseases:") # 90 segundos
diseases <- DEGsToDiseases(mejores_genes_knn, getEvidences = TRUE, size = 10)
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
         enfermedades[grep(pattern = "alcoh",  ignore.case = T, x = enfermedades)]))

# Session_info para reproducibilidad
devtools::session_info()