# Carga de paquetes
library(readr)
library(dplyr)

# Lectura de fichero con datos clínicos
clinical <- read_delim("data/clinical.cart.2020-06-11/clinical.tsv",
                       "\t", escape_double = FALSE, trim_ws = TRUE)

# Diccionario de variables:
# https://docs.gdc.cancer.gov/Data_Dictionary/viewer/#?view=table-definition-view&id=diagnosis

clinical <- clinical %>% 
  # Selección de variables
  select(case_id, tumor_stage, age_at_diagnosis, gender, ethnicity, race, vital_status, year_of_diagnosis, primary_diagnosis) %>% 
  # Eliminar duplicados
  unique

# Verificación: todos los casos de "Primary Tumor" tienen un registro clínico
length(clinical$case_id)
length(unique(clinical$case_id))

# Tabla con años de diagnóstico
table(clinical$year_of_diagnosis)

# Análisis estadio
table(clinical$tumor_stage, useNA = "always")
# Recodificación a estadios simples
clinical$estadio <- 0
clinical[which(clinical$tumor_stage == "stage i"), "estadio"] <- 1
clinical[which(clinical$tumor_stage == "stage ia"), "estadio"] <- 1
estadio2 <- c("stage ii", "stage iia", "stage iib", "stage iic")
clinical[which(clinical$tumor_stage %in% estadio2), "estadio"] <- 2
estadio3 <- c("stage iii", "stage iiia", "stage iiib", "stage iiic")
estadio4 <- c("stage iv", "stage iva", "stage ivb")
clinical[which(clinical$tumor_stage %in% estadio3), "estadio"] <- 3
clinical[which(clinical$tumor_stage %in% estadio4), "estadio"] <- 4
clinical[which(clinical$tumor_stage == "not reported"), "estadio"] <- NA
clinical[which(clinical$estadio == 0), "estadio"] <- NA
# Tabla recodificada
table(clinical$estadio, useNA = "always")

# Cálculo y análisis de edad al diagnóstico
clinical$age_at_diagnosis <- as.numeric(clinical$age_at_diagnosis) / 365.25
summary(clinical$age_at_diagnosis)

# Cálculo y análisis de grupos de edad
clinical$gedad <- cut(clinical$age_at_diagnosis, breaks = c(0, 40, 50, 60, 70, 80, 100),
                      include.lowest = F)
levels(clinical$gedad)
table(clinical$gedad, useNA = "always")

# Análisis de género
table(clinical$gender, useNA = "always")

# Análisis de raza 
table(clinical$race, useNA = "always")

# Análisis de estado vital
table(clinical$vital_status, useNA = "always")

# Recodificación y análisis de diagnóstico primario
table(clinical$primary_diagnosis, useNA = "always") %>%  sort(decreasing = T)
clinical[! clinical$primary_diagnosis %in% c(
 "Adenocarcinoma, NOS",
 "Mucinous adenocarcinoma"
), "primary_diagnosis"] <- "otros"
# Tabla final de diagnóstico primario
table(clinical$primary_diagnosis, useNA = "always")

# Tablas de estado vital (análisis completo, quitando datos faltantes)
clinical <- clinical %>% filter(vital_status != "Not Reported")

# Tabla estado vital - sexo + chi cuadrado
round(table(clinical$gender, clinical$vital_status, useNA = "always"), 1)
round(prop.table(table(clinical$gender, clinical$vital_status, useNA = "always"), margin = 1) * 100, 1)

chisq.test(clinical$gender, clinical$vital_status)

# Tabla estado vital - grupo de edad + chi cuadrado
clinical2 <- clinical %>% filter(is.na(gedad) != T)
round(table(clinical2$gedad, clinical2$vital_status, useNA = "always"), 1)
round(prop.table(table(clinical2$gedad, clinical2$vital_status, useNA = "always"), margin = 1) * 100, 1)

chisq.test(clinical2$gedad, clinical2$vital_status)

# Tabla estado vital - diagnóstico primario + chi cuadrado
round(table(clinical$primary_diagnosis, clinical$vital_status, useNA = "always"), 1)
round(prop.table(table(clinical$primary_diagnosis, clinical$vital_status, useNA = "always"), margin = 1) * 100, 1)

chisq.test(clinical$primary_diagnosis, clinical$vital_status)

# Tabla estado vital - estadio + chi cuadrado
clinical3 <- clinical %>% filter(is.na(estadio) != T)
round(table(clinical3$estadio, clinical3$vital_status, useNA = "always"), 1)
round(prop.table(table(clinical3$estadio, clinical3$vital_status, useNA = "always"), margin = 1) * 100, 1)

chisq.test(clinical3$estadio, clinical3$vital_status)
