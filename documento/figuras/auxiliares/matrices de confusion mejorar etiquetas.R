library(KnowSeq)

data <- table(as.factor(c(rep("carc. hep.", 87), rep("colang.", 7), rep("tejido_normal", 14), rep("colang.", 3), "carc. hep.")),
              as.factor(c(rep("carc. hep.", 87), rep("colang.", 7), rep("tejido_normal", 14), rep("carc. hep.", 3), "colang.")))
plotConfMatrix(data)


png(filename = "28_higado_multiclase_15_svm_matriz_confusion_mejor_metodo_mejorado.png", width = 13, height = 8, units = "in", res = 300)
plotConfMatrix(data)
dev.off()


data <- table(as.factor(c(rep("carc. hep.", 88), rep("colang.", 7), rep("tejido_normal", 12), rep("colang.", 2), "carc. hep.", rep("carc. hep.", 2))),
              as.factor(c(rep("carc. hep.", 88), rep("colang.", 7), rep("tejido_normal", 12), rep("carc. hep.", 2), "colang.", rep("tejido_normal", 2))))
plotConfMatrix(data)


png(filename = "29_higado_multiclase_25_rf_matriz_confusion_mejor_metodo_mejorado.png", width = 13, height = 8, units = "in", res = 300)
plotConfMatrix(data)
dev.off()
