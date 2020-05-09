cuartiles <- function(icdo)
{
  # Se selecciona Granada (18 = 101) y los datos de los cuatro primeros periodos de NORDPRED (años añoIP1 hasta año añoFP4)
  datos <- mortalidad[mortalidad$PROVINCIA==101 & mortalidad$ICDO==icdo & mortalidad$AÑO>=añoIP1 & mortalidad$AÑO<=añoFP4,]	
  datos <- aggregate(datos[,4:22],list(ICDO=datos$ICDO),sum)
  datos <- datos[,3:20]
  suma <- sum(datos)
  frecuencia <- (datos/suma)*100
  delta1 <- 100
  delta2 <- 100
  delta3 <- 100
  delta4 <- 100
  
  # Se calculan P10, T1, T2 y P90.
  # corteedad es el número del grupo de edad (de 1 a 18) que corresponde al P10.
  for (i in 1:18) {
    suma.sup <- sum(frecuencia[,1:i])
    if (abs(10-suma.sup) < delta1) {
      delta1 <- abs(10-suma.sup)
      corteedad <- i
      p10 <- i*5
    }
    if (abs(33.3333-suma.sup) < delta2) {
      delta2 <- abs(33.3333-suma.sup)
      t2 <- i*5
    }
    if (abs(66.6666-suma.sup) < delta3) {
      delta3 <- abs(66.6666-suma.sup)
      t3 <- i*5
    }
    if (abs(90-suma.sup) < delta4) {
      delta4 <- abs(90-suma.sup)
      p90 <- i*5
    }
  }
  aux <- c(icdo, corteedad,p10,t2,t3,p90)
  aux
}	