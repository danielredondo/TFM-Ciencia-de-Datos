crear.parametros.tumor.RCG<-function(){
  salida <- cuartiles(114)
  salida <- rbind(salida, cuartiles(115))
  salida <- rbind(salida, cuartiles(116))
  salida <- rbind(salida, cuartiles(117))
  salida <- rbind(salida, cuartiles(118))
  salida <- rbind(salida, cuartiles(119))
  salida <- rbind(salida, cuartiles(199))
  salida <- rbind(salida, cuartiles(614))
  salida <- rbind(salida, cuartiles(615))
  salida <- rbind(salida, cuartiles(616))
  salida <- rbind(salida, cuartiles(617))
  salida <- rbind(salida, cuartiles(618))
  salida <- rbind(salida, cuartiles(619))
  salida <- rbind(salida, cuartiles(699))
  # para respetar la estructura de parametros.tumor
  subsalida1 <- salida[c(1:6), c(2:6)]
  subsalida2 <- salida[7, c(2:6)]
  subsalida3 <- salida[c(8:13), c(2:6)]
  subsalida4 <- salida[14, c(2:6)]
  cero5 <- rep(0, 5)
  salida <- rbind(subsalida1, cero5, subsalida2, cero5, subsalida3, cero5, subsalida4, cero5)
  icdo <- c(114, 115, 116, 117, 118, 119, 121, 199, 100, 614, 615, 616, 617, 618, 619, 621, 699, 600)
  sexo <- c("M", "M", "M", "M", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F", "F", "F", "F", "F")
  loc <- c("Colon", "Recto", "Pulmón", "Próstata", "Vejiga", "Estómago", "Colon-recto", "Otros", "Total excepto piel no melanoma", "Colon", "Recto", "Pulmón", "Mama", "Cuerpo uterino", "Ovario", "Colon-recto", "Otros", "Total excepto piel no melanoma")
  tipo <- c("SP", "SP", "SP", "SP", "SP", "SP", "0", "SP", "0", "SP", "SP", "SP", "SP", "SP", "SP", "0", "SP", "0")
  escenario <- c("A", "A", "A", "A", "A", "A", "Especial", "A", "Especial", "A", "A", "A", "A", "A", "A", "Especial", "A", "Especial")
  salida <- cbind(icdo, sexo, loc,salida, tipo,escenario)
  cabecera <- c("CODIGO", "SEXO", "VALOR", "CORTEEDAD", "P10", "T2", "T3", "P90", "TIPO", "ESCENARIO")
  colnames(salida) <- cabecera
  setwd(paste0(ruta, "\\datos"))
  write.csv2(salida, "parametros.tumor.csv", row.names = FALSE, quote = FALSE)
  print(paste("Se ha creado parametros.tumor.csv"), quote = FALSE)
}