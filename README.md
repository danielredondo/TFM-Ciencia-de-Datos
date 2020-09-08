# Epidemiología y detección de biomarcadores en cáncer

Trabajo Fin de Máster. Máster Universitario Oficial en Ciencia de Datos e Ingeniería de Computadores. Curso 2019/20.

Universidad de Granada.

**Daniel Redondo Sánchez.**

**Tutores: Daniel Castillo Secilla, Luis Javier Herrera Maldonado.**

## Descripción del contenido del repositorio

- `analisis_cr`: Código, resultados y figuras del análisis realizado para cáncer de colon-recto.

- `analisis_higado`: Código, resultados y figuras del análisis realizado para cáncer de hígado.

- `documento`: Documento del TFM escrito en *LaTeX*. La versión final del documento en PDF está [aquí](https://github.com/danielredondo/TFM_ciencia_de_datos/blob/master/documento/documento_entregado.pdf).

- `epidemiologia`: Datos sobre la epidemiología del cáncer.

- `presentacion`: Presentación para la defensa del TFM.

- `shiny`: Aplicación web `biomarkeRs`, disponible online en [https://dredondo.shinyapps.io/biomarkeRs/](https://dredondo.shinyapps.io/biomarkeRs/).

## Resumen

**Introducción:** El cáncer es uno de los mayores problemas de salud pública del mundo con más de 17 millones de casos nuevos y 9 millones de defunciones al año. 

**Métodos:** Este trabajo se centra en el cáncer de hígado y el cáncer de colon-recto, describiendo sus principales indicadores epidemiológicos y usando *machine learning* para analizar más de 1.100 muestras de RNA-Seq procedentes de pacientes de cáncer. Para clasificación biclase (tumor vs. tejido normal) y multiclase (varios tipos de tumor vs. tejido normal) se identifican los 10 genes más relevantes y se construyen modelos predictivos con SVM, *random forest* y kNN con validación cruzada 5-fold.

**Resultados:** Los mejores clasificadores biclase para cada algoritmo son validados con excelentes resultados para cáncer de hígado (F1-Score en test: 99,5%) y cáncer de colon-recto (F1-Score: 100%). En los mejores modelos para clasificación multiclase se obtienen medidas de evaluación inferiores tanto en hígado (F1-Score: 91,8%) como en colon-recto (F1-Score: 79,3%).

Se ha desarrollado una aplicación web, `biomarkeRs`, que implementa análisis de transcriptómica y puede resultar de utilidad para usuarios sin conocimientos previos de programación.

**Conclusiones:** SVM, random forest y kNN obtienen resultados muy similares, y consiguen distinguir correctamente entre tejidos tumorales y sanos, con algunos problemas para distinguir entre diferentes tipos de cáncer. Es necesaria una validación externa e interpretaciones clínicas para establecer de forma clara una asociación gen-enfermedad.

**Palabras clave:** epidemiología, transcriptómica, cáncer de hígado, cáncer de colon-recto, RNA-Seq, machine learning, SVM, random forest, kNN.

## Abstract
**Introduction:** Cancer is one of the world's largest public health problems with more than 17 million new cases and 9 million deaths every year. 

**Methods:** This work focuses on liver cancer and colon-rectum cancer, describing their main epidemiological indicators and using machine learning to analyze more than 1,100 RNA-Seq samples from cancer patients. For binary (tumor vs. normal tissue) and multiclass (various tumor types vs. normal tissue) classification, the 10 most relevant genes are identified and predictive models are constructed with SVM, random forest and kNN with 5-fold cross-validation. 

**Results:** The best binary classifiers are validated with excellent results for liver cancer (F1-Score in test: 99.5%) and colon-rectum cancer (F1-Score: 100%). Lower evaluation measures are obtained in the best models for multiclass classification, in both liver (F1-Score: 91.8%) and colon-rectum (F1-Score: 79.3%). 

A web application has been developed, `biomarkeRs`, that implements transcriptomic analysis and can be useful for users with no previous knowledge of programming. 

**Conclusions:** SVM, random forest and kNN obtained very similar results, and managed to correctly distinguish between tumoral and normal tissues with some troubles distinguishing between different types of cancer. External validation and clinical interpretations are necessary to clearly establish a gene-disease association.


**Keywords:** epidemiology, transcriptomics, liver cancer, colorectal cancer, RNA-Seq, machine learning, SVM, random forest, kNN.