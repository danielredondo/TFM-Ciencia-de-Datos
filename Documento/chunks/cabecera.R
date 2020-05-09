# 240 para 500 simulaciones, 24000 para 50000 simulaciones
ciclos.seleccionados <- 24000 

ruta <- "RUTA"

año.a.estimar <- 2013 

# Años para evolución RIM
total.años.evolucion.rim <- 15
añoFrim <- año.a.estimar - 5
añoIrim <- añoFrim - total.años.evolucion.rim + 1
# Por ejemplo:
# Si año.a.estimar es 2017, y se consideran 15 años de evolución que acaban en 2012, el año inicial será 1998

# Años NORDPRED (I = Año inicial, F = Año final, P = Periodo)
añoIP1 <- año.a.estimar - 22
añoFP1 <- año.a.estimar - 18
añoIP2 <- año.a.estimar - 17
añoFP2 <- año.a.estimar - 13
añoIP3 <- año.a.estimar - 12
añoFP3 <- año.a.estimar - 8
añoIP4 <- año.a.estimar - 7
añoFP4 <- año.a.estimar - 3

# P5 es el periodo del que obtenemos la mortalidad tras aplicar NORDPRED a los periodos del 1 al 4
añoIP5 <- año.a.estimar - 2
añoFP5 <- año.a.estimar + 2