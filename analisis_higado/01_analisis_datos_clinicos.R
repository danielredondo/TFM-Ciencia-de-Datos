library(readr)
library(dplyr)

clinical <- read_delim("data/clinical.cart.2020-06-11/clinical.tsv",
                       "\t", escape_double = FALSE, trim_ws = TRUE)

# Diccionario de variables: https://docs.gdc.cancer.gov/Data_Dictionary/viewer/#?view=table-definition-view&id=diagnosis

clinical <- clinical %>% 
  select(case_id, tumor_stage, age_at_diagnosis, gender, ethnicity, race, vital_status, year_of_diagnosis) %>% 
  unique

# Años de diagnóstico
clinical$year_of_diagnosis %>% table

# Todos los casos de "Primary Tumor" tienen un registro clínico
length(clinical$case_id)
length(unique(clinical$case_id))

# Análisis estadio
table(clinical$tumor_stage, useNA = "always")
# Recodificación a estadios simples
clinical$estadio <- 0
clinical[which(clinical$tumor_stage == "stage i"), "estadio"] <- 1
clinical[which(clinical$tumor_stage == "stage ii"), "estadio"] <- 2
estadio3 <- c("stage iiia", "stage iiib", "stage iiic")
estadio4 <- c("stage iva", "stage ivb")
clinical[which(clinical$tumor_stage %in% estadio3), "estadio"] <- 3
clinical[which(clinical$tumor_stage %in% estadio4), "estadio"] <- 4
clinical[which(clinical$tumor_stage == "'--"), "estadio"] <- NA
clinical[which(clinical$estadio == 0), "estadio"] <- NA

table(clinical$estadio, useNA = "always")

# Análisis edad al diagnóstico
clinical$age_at_diagnosis <- as.numeric(clinical$age_at_diagnosis) / 365.25
summary(clinical$age_at_diagnosis)

# Grupos de edad
clinical$gedad <- cut(clinical$age_at_diagnosis, breaks = c(0, 40, 50, 60, 70, 80, 100),
                      include.lowest = F)
levels(clinical$gedad)
table(clinical$gedad, useNA = "always")

# Análisis género
table(clinical$gender, useNA = "always")

# Análisis etnia
table(clinical$ethnicity, useNA = "always")

# Análisis raza 
table(clinical$race, useNA = "always")

# Análisis estadio vital
table(clinical$vital_status, useNA = "always")

# Tablas estadio vital (análisis completo, quitando datos faltantes)

clinical <- clinical %>% filter(vital_status != "Not Reported")

# Estado vital - sexo
round(table(clinical$gender, clinical$vital_status, useNA = "always"), 1)
round(prop.table(table(clinical$gender, clinical$vital_status, useNA = "always"), margin = 1) * 100, 1)

chisq.test(clinical$gender, clinical$vital_status)

# Estado vital - gedad
clinical2 <- clinical %>% filter(is.na(gedad) != T)
round(table(clinical2$gedad, clinical2$vital_status, useNA = "always"), 1)
round(prop.table(table(clinical2$gedad, clinical2$vital_status, useNA = "always"), margin = 1) * 100, 1)

chisq.test(clinical2$gedad, clinical2$vital_status)

# Estado vital - estadio
clinical3 <- clinical %>% filter(is.na(estadio) != T)
round(table(clinical3$estadio, clinical3$vital_status, useNA = "always"), 1)
round(prop.table(table(clinical3$estadio, clinical3$vital_status, useNA = "always"), margin = 1) * 100, 1)

chisq.test(clinical3$estadio, clinical3$vital_status)
