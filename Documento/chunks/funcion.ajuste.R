funcion.ajuste<-function(localizacion, escenario, sims = 100, nivelconfianza, imprimir = "FALSE"){
  
  if(escenario=="TODOS"){
    salida <- funcion.ajuste(localizacion,"E3",sims,nivelconfianza)
    salida <- rbind(salida, funcion.ajuste(localizacion, "E1", sims, nivelconfianza))
    salida <- rbind(salida, funcion.ajuste(localizacion, "E5", sims, nivelconfianza))
    salida <- rbind(salida, funcion.ajuste(localizacion, "L", sims, nivelconfianza))
    salida <- rbind(salida, funcion.ajuste(localizacion, "C", sims, nivelconfianza))
    
    salida
  }
  
  else{	
    # 121 = Colon-recto en hombres
    # 621 = Colon-recto en mujeres
    if(localizacion != 121 && localizacion != 621){
      observado <- estimado.observado(localizacion, escenario, sims, nivelconfianza)$observado
      estimado  <- estimado.observado(localizacion, escenario, sims, nivelconfianza)$estimado
    }
    else{
      if(localizacion==121){
        observado.colon <- estimado.observado(114, escenario, sims, nivelconfianza)$observado
        observado.recto <- estimado.observado(115, escenario, sims, nivelconfianza)$observado
        observado <- observado.colon+observado.recto
        estimado.colon <- estimado.observado(114, escenario, sims, nivelconfianza)$estimado
        estimado.recto <- estimado.observado(115, escenario, sims, nivelconfianza)$estimado
        estimado <- estimado.colon+estimado.recto
      }
      else{
        observado.colon <- estimado.observado(614, escenario, sims, nivelconfianza)$observado
        observado.recto <- estimado.observado(615, escenario, sims, nivelconfianza)$observado
        observado <- observado.colon+observado.recto
        estimado.colon <- estimado.observado(614, escenario, sims, nivelconfianza)$estimado
        estimado.recto <- estimado.observado(615, escenario, sims, nivelconfianza)$estimado
        estimado <- estimado.colon+estimado.recto		
      }
    }
    
    if(imprimir == "TRUE"){
      print("Observado:", quote = "FALSE")
      print(observado)
      print("Estimado:", quote = "FALSE")
      print(estimado)
    }
    
    # MAPE
    f <- 100 * (1 / 10) * sum(abs((estimado - observado) / observado))
    icdo <- parametros.tumor[parametros.tumor$CODIGO == localizacion, "VALOR"]
    data.frame(localizacion, icdo, escenario, f)
    
  }	
}
}