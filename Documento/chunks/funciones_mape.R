mejor.escenario.ajuste <- function(localizacion, sims = 100, nivelconfianza, imprimir = "FALSE"){
  esc.e1 <- funcion.ajuste(localizacion, "E1", sims, nivelconfianza, "FALSE")[[4]]
  esc.e3 <- funcion.ajuste(localizacion, "E3", sims, nivelconfianza, "FALSE")[[4]]
  esc.e5 <- funcion.ajuste(localizacion, "E5", sims, nivelconfianza, "FALSE")[[4]]
  esc.l  <- funcion.ajuste(localizacion, "L" , sims, nivelconfianza, "FALSE")[[4]]
  esc.c  <- funcion.ajuste(localizacion, "C" , sims, nivelconfianza, "FALSE")[[4]]
  
  #Vemos cuál es el mejor escenario
  if(esc.e1<esc.e3 && esc.e1<esc.e5 && esc.e1<esc.l && esc.e1<esc.c) esc.mejor="E1" else{
    if(esc.e3<esc.e1 && esc.e3<esc.e5 && esc.e3<esc.l && esc.e3<esc.c) esc.mejor="E3" else{
      if(esc.e5<esc.e1 && esc.e5<esc.e3 && esc.e5<esc.l && esc.e5<esc.c) esc.mejor="E5" else{
        if(esc.l<esc.e1 && esc.l<esc.e5 && esc.l<esc.e3 && esc.l<esc.c) esc.mejor="L" else{
          esc.mejor="C"
        }
      }
    }
  }
  
  ifelse(imprimir == "FALSE",
         salida <- paste0(funcion.ajuste(localizacion, esc.mejor, sims, nivelconfianza, "FALSE")[[3]]),
         salida <- funcion.ajuste(localizacion, esc.mejor, sims, nivelconfianza, "FALSE"))
  salida
}

funcion.ajuste.total <- function(sims = 100, nivelconfianza){
  
  # Hombres
  salida <- funcion.ajuste(114,"TODOS",sims,nivelconfianza)
  salida <- rbind(salida, funcion.ajuste(115, "TODOS", sims,nivelconfianza))
  salida <- rbind(salida, funcion.ajuste(116, "TODOS", sims,nivelconfianza))
  salida <- rbind(salida, funcion.ajuste(117, "TODOS", sims,nivelconfianza))
  salida <- rbind(salida, funcion.ajuste(118, "TODOS", sims,nivelconfianza))
  salida <- rbind(salida, funcion.ajuste(119, "TODOS", sims,nivelconfianza))
  salida <- rbind(salida, funcion.ajuste(199, "TODOS", sims,nivelconfianza))
  
  # Mujeres
  salida <- rbind(salida, funcion.ajuste(614, "TODOS", sims, nivelconfianza))
  salida <- rbind(salida, funcion.ajuste(615, "TODOS", sims, nivelconfianza))
  salida <- rbind(salida, funcion.ajuste(616, "TODOS", sims, nivelconfianza))
  salida <- rbind(salida, funcion.ajuste(617, "TODOS", sims, nivelconfianza))
  salida <- rbind(salida, funcion.ajuste(618, "TODOS", sims, nivelconfianza))
  salida <- rbind(salida, funcion.ajuste(619, "TODOS", sims, nivelconfianza))
  salida <- rbind(salida, funcion.ajuste(699, "TODOS", sims, nivelconfianza))
  
  salida
  
}

tabla.mejores.escenarios<-function(sims=100,nivelconfianza){
  v<-funcion.ajuste.total(sims,nivelconfianza)
  v<-v[order(v$f),]
  
  # Hombres
  salida <- t(subset(v, v$localizacion == "114"))
  salida <- rbind(salida, t(subset(v, v$localizacion == "115")))
  salida <- rbind(salida, t(subset(v, v$localizacion == "116")))
  salida <- rbind(salida, t(subset(v, v$localizacion == "117")))
  salida <- rbind(salida, t(subset(v, v$localizacion == "118")))
  salida <- rbind(salida, t(subset(v, v$localizacion == "119")))
  salida <- rbind(salida, t(subset(v, v$localizacion == "199")))
  
  # Mujeres
  salida <- rbind(salida, t(subset(v, v$localizacion == "614")))
  salida <- rbind(salida, t(subset(v, v$localizacion == "615")))
  salida <- rbind(salida, t(subset(v, v$localizacion == "616")))
  salida <- rbind(salida, t(subset(v, v$localizacion == "617")))
  salida <- rbind(salida, t(subset(v, v$localizacion == "618")))
  salida <- rbind(salida, t(subset(v, v$localizacion == "619")))
  salida <- rbind(salida, t(subset(v, v$localizacion == "699")))
  
  print(salida, quote = "FALSE")
}


