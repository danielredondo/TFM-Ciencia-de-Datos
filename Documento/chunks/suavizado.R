# Lectura de ficheros
incidencia <- read.csv2("incidencia.csv", fileEncoding="latin1")
mortalidad <- read.csv2("mortalidad.csv", fileEncoding="latin1")
poblacion <- read.csv2("poblacion.csv", fileEncoding="latin1")
parametros.registro <- read.csv2("parametros.registro.csv")

# Casos observados 2004-2013 hombres ----
o114 <- c(138, 203, 156, 180, 172, 227, 208, 214, 243, 222)
o115 <- c(90, 97, 93, 95, 115, 136, 112, 110, 115, 121)
o116 <- c(307, 350, 332, 344, 377, 317, 348, 312, 353, 351) 
o117 <- c(327, 351, 405, 442, 478, 498, 463, 519, 494, 447) 
o118 <- c(257, 245, 245, 253, 251, 275, 264, 265, 312, 268) 
o119 <- c(87, 70, 89, 76, 94, 76, 83, 67, 86, 73) 
o199 <- c(784, 891, 824, 872, 876, 939, 892, 950, 924, 947) 

# Casos observados 2004-2013 mujeres ----
o614 <- c(117, 138, 135, 153, 144, 153, 152, 144, 168, 163)
o615 <- c(67, 61, 75, 60, 67, 70, 70, 77, 66, 56)
o616 <- c(45, 52, 51, 53, 48, 62, 59, 71, 60, 56)
o617 <- c(338, 345, 389, 410, 415, 494, 419, 466, 450, 494) 
o618 <- c(114, 106, 99, 107, 103, 123, 127, 124, 107, 124) 
o619 <- c(60, 77, 58, 67, 68, 69, 62, 50, 53, 51)
o699 <- c(718, 763, 703, 822, 790, 838, 850, 856, 791, 858)

# Selección de localización anatómica y número de alphas a incluir en el análisis
loc <- 699
numero_alfas <- 20 # 20 -> alpha = 0.05, ..., 0.95, 1.00
observado <- eval(as.name(paste0("o", loc)))


for(i in 0:numero_alfas){
  
  # Cálculo de alpha
  alfa <- i / numero_alfas 
  
  if(i == 0){
    vector <- NA # Se crea el vector que almacenará los MAPE
    alfa <- NULL # El primer caso estima el coeficiente alpha
  }
  
  for(año.a.estimar in 2004:2013){
    
    # Parámetros para estimar la mortalidad para el año a estimar
    añoFrim <- año.a.estimar - 5
    añoIrim <- 1985
    crear.parametros.tumor.RCG()
    crear.parametros.registro.RCG()
    añoIP1 <- año.a.estimar - 22
    añoFP1 <- año.a.estimar - 18
    añoIP2 <- año.a.estimar - 17
    añoFP2 <- año.a.estimar - 13
    añoIP3 <- año.a.estimar - 12
    añoFP3 <- año.a.estimar - 8
    añoIP4 <- año.a.estimar - 7
    añoFP4 <- año.a.estimar - 3
    añoIP5 <- año.a.estimar - 2
    añoFP5 <- año.a.estimar + 2
    
    # Selección de incidencia y mortalidad, creación de RIM
    incidencia.seleccion <- subset(incidencia,
                             incidencia$AÑO %in% añoIrim:añoFrim &
                               !incidencia$ICDO %in% c(121, 621))[, -23]
    mortalidad.seleccion <- subset(mortalidad,
                             mortalidad$AÑO %in% añoIrim:añoFrim &
                               mortalidad$PROVINCIA == 18)
    tabla.rim <- merge(incidencia.seleccion, mortalidad.seleccion)
    tabla.rim$R <- (tabla.rim$I01 + tabla.rim$I02 +
          tabla.rim$I03 + tabla.rim$I04 + tabla.rim$I05 +
          tabla.rim$I06 + tabla.rim$I07 + tabla.rim$I08 +
          tabla.rim$I09 + tabla.rim$I10 + tabla.rim$I11 +
          tabla.rim$I12 + tabla.rim$I13 + tabla.rim$I14 +
          tabla.rim$I15 + tabla.rim$I16 + tabla.rim$I17 +
          tabla.rim$I18) /
          (tabla.rim$M01 + tabla.rim$M02 + tabla.rim$M03 +
          tabla.rim$M04 + tabla.rim$M05 + tabla.rim$M06 +
          tabla.rim$M07 + tabla.rim$M08 + tabla.rim$M09 +
          tabla.rim$M10 + tabla.rim$M11 + tabla.rim$M12 +
          tabla.rim$M13 + tabla.rim$M14 + tabla.rim$M15 +
          tabla.rim$M16 + tabla.rim$M17 + tabla.rim$M18)
    
    # Selección de RIM para la localización seleccionada
    a <- subset(tabla.rim$R, tabla.rim$ICDO == loc)
    
    # Modelo de alisamiento exponencial
    hw <- HoltWinters(a, alpha = alfa, gamma = FALSE, beta = FALSE)
    
    # Predicción de 5 años (para llegar hasta el año a estimar)
    hw.p <- predict(hw, n.ahead = 5)
    
    # NORDPRED para el año a estimar
    prep.nordpred.lista()
    prep.nordpred.sexo("M")
    prep.nordpred.sexo("F")
    escribir.prediccion.año.a.estimar()
    prediccion.mortalidad <- read.csv2("prediccion.año.a.estimar.csv")
    
    M <- subset(prediccion.mortalidad, prediccion.mortalidad$icdo == loc)
    M <- M$M01 + M$M02 + M$M03 + M$M04 + M$M05 + M$M06 +
         M$M07 + M$M08 + M$M09 + M$M10 + M$M11 + M$M12 +
         M$M13 + M$M14 + M$M15 + M$M16+ M$M17 + M$M18
    
    # Cálculo de casos incidentes estimados: mortalidad * RIM
    total <- M * hw.p[5, ]
    
    # Se anexan los casos estimados año a año
    ifelse(año.a.estimar == 2004,
           total_periodos <- total,
           total_periodos <- rbind(total_periodos, total))
  }
  
  # Se imprime el MAPE, el valor de alpha, y los números de casos estimados
  print(paste0("MAPE = ",
               round(10 * sum(abs(total_periodos - observado) / observado), 1)))
  print(paste0("alfa = ", alfa))
  print("__________________")
  
  # Se guarda el valor del MAPE
  vector <- rbind(vector, 10 * sum(abs(total_periodos - observado) / observado))
  
  if(i == numero_alfas){   # Al terminar el proceso:
    
    # Se quita el primer elemento, que es NA
    vector <- vector[-1, ]
    
    # Se extrae el primer valor, resultante de la estimación automática de alpha
    alpha_estimado <- vector[1]
    vector <- vector[-1]
    
    # El resto de valores se guardan en un vector. Se añaden los alphas.
    vector <- cbind(c(1:numero_alfas) / numero_alfas, vector)
    colnames(vector) <- c("al", "v")
    
    # Se representará en rojo el valor más bajo de MAPE
    colores <- rep("black", numero_alfas)
    for(i in 1:numero_alfas){
      if(vector[i, 2] == min(vector[, 2])) colores[i] <- "red"
    }
    
    # Se muestra la tabla con coeficientes de suavizado, MAPE,
    # y color para detectar el MAPE más bajo
    print(cbind(vector, colores))
    
    # Mostramos el gráfico
    plot(vector, main = paste0(loc), pch = 20, cex = 1.5,
         col = colores, ylab = "MAPE", xlab = "alpha")
    abline(h = alpha_estimado, lty = 3)
    
    # Guardamos el gráfico como PDF
    pdf(paste0(loc, "-", numero_alfas, ".pdf"))
    plot(vector, main = paste0(loc), pch = 20, cex = 1.5,
         col = colores, ylab = "MAPE", xlab = "alpha")
    abline(h = alpha_estimado, lty = 3)
    dev.off()
  }
}