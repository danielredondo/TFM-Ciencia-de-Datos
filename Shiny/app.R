# Mejoras:
# Añadir tumores y sanos en la parte derecha del gráfico de Sankey
# da debe tener en cuenta la enfermedad!!

library(shiny)
library(shinydashboard)
library(dplyr)
library(KnowSeq)
library(reshape2)
library(caret)
library(ggplot2)
library(ggalluvial)

source("www/dataPlot.R")

ui <- dashboardPage(
  ## Tema
  skin = "black",
  ## Cabecera
  dashboardHeader(title = "biomaRcadores"#,titleWidth = 350
                  ),
  ## Barra lateral
  dashboardSidebar(
    sidebarMenu(
      menuItem("Introducción", tabName = "intro", icon = icon("file-alt")),
      menuItem("Carga de datos", tabName = "datos", icon = icon("database")),
      menuItem("Selección de genes", tabName = "genes", icon = icon("dna")),
      menuItem("Entrenamiento de modelos", tabName = "entrenamiento", icon = icon("play")),
      menuItem("Validación de modelos", tabName = "validacion", icon = icon("check-circle")),
      menuItem("Autor", tabName = "autor", icon = icon("id-card")),
      menuItem("Código", tabName = "codigo", icon = icon("code"))
    )
  ),
  ## Cuerpo
  dashboardBody(
    tabItems(
      # Tab 1
      tabItem(tabName = "intro",
              h2("Introducción App + abstract TFM"),
              
              h4(tags$b("Introducción")),
              "Texto",
              h4(tags$b("Métodos")),
              "Texto",
              h4(tags$b("Resultados")),
              "Texto",
              h4(tags$b("Conclusiones")),
              "Texto",
              
              h4(tags$b("Sobre esta aplicación")),
              "Texto",
              # Parte final
              br(), br(), br(),
              fluidRow(column(6, tags$img(src = "ugr.png", height = "150px")),
                       column(6, tags$img(src = "knowseq.png", height = "150px")))
              
      ),
      
      # Tab 2
      tabItem(tabName = "datos",
              h2("Carga de datos"),
              fileInput(inputId = "archivo_rdata",
                        label = "Seleccione el archivo .RData a importar",
                        buttonLabel = "Examinar...",
                        accept = ".RData",
                        placeholder = "No se ha seleccionado ningún archivo",
                        width = "50%"
              ),
              
              actionButton(inputId = "boton_importar",
                           label = "Importar archivo",
                           icon = icon("fas fa-file-import", lib = "font-awesome"),
                           width = "50%"),
              br(),
              
              tableOutput("tabla1"),
              
              h2("Partición entrenamiento-test"),
              
              sliderInput("porcentaje_entrenamiento",
                          label = "Porcentaje entrenamiento (%)",
                          value = 75, min = 5, max = 95, step = 5,
                          width = "50%"
                          ),
              
              plotOutput("sankey", width = "75%")
      ),
      
      # Tab 3
      tabItem(tabName = "genes",
              h2("Selección de genes"),
              sliderInput(inputId = "numero_genes", label = "Selecciona el número de genes relevantes", value = 20, min = 0, max = 50, step = 1),
              "Se muestran a continuación los mejores genes seleccionados por cada método de selección de características.",
              br(),
              tags$i("Puede tardar unos segundos en actualizarse."),
              fluidRow(
                column(4, h4(tags$b("  MRMR")), tableOutput("genes_mrmr")),
                column(4, h4(tags$b("  RF")), tableOutput("genes_rf")),
                column(4, h4(tags$b("  DA")), tableOutput("genes_da")),
              ),
                

              br(),
              
              h2("Enfermedades relacionadas")
      ),
      
      # Tab 4
      tabItem(tabName = "entrenamiento",
              h2("Entrenamiento de modelos"),
              
              # Elegir MRMR/RF/DA, se muestra MRMR por ahora
              selectInput("tipo_entrenamiento",
                          label = "Tipo de selección de genes",
                          choices = c("mRMR", "rf", "da"),
                          selected = "mRMR",
                          width = "50%"),
              
              # Se puede añadir número de CV que se hacen
              selectInput("numero_folds",
                          label = "Número de folds",
                          choices = c(3, 5, 10),
                          selected = 5,
                          width = "50%"),
              
              tableOutput("mejores_parametros"),
              
              plotOutput("resultados_entrenamiento")
              
              # Falta añadir el F1-score y decir cuál es el mejor método

              
      ),
      # Tab 5
      tabItem(tabName = "validacion",
              h2("Validación de modelos")
      ),
      # Tab 6
      tabItem(tabName = "autor",
              h2("Autor")
      ),
      # Tab 7
      tabItem(tabName = "codigo",
              h2("Código")
      )
      ) # Final tabs
  ) # Final dashboard body
) # Final dashboard page

# Ampliar tamaño de .RData a 15MB en lugar de los 5MB por defecto
options(shiny.maxRequestSize = 15*1024^2)

server <- function(input, output){

  observeEvent(input$boton_importar, {
    
    # Si se ha seleccionado un fichero, se importa
    load(input$archivo_rdata$datapath)
    # Extraer labels
    labels <- matriz[1, ] %>% as.vector
    # Extraer matriz
    DEGsMatrix <- matriz[2:nrow(matriz), ]
    filas <- rownames(DEGsMatrix)
    DEGsMatrix <- apply(DEGsMatrix, 2, as.numeric)
    rownames(DEGsMatrix) <- filas
    # Crear DEGsMatrixML
    DEGsMatrixML <- t(DEGsMatrix)
    
    # Parámetros generales
    
    # Partición 75% / 25% con balanceo de clase
    set.seed(1991)
    indices <- reactive(createDataPartition(labels, p = input$porcentaje_entrenamiento / 100, list = FALSE))
    particion <- reactive(list(training = DEGsMatrixML[indices(), ], test = DEGsMatrixML[-indices(), ]))
    
    # Conjuntos
    particion.entrenamiento <- reactive(particion()$training)
    particion.test <- reactive(particion()$test)
    
    # Etiquetas
    labels_train <- reactive(labels[indices()])
    labels_test  <- reactive(labels[-indices()])
    
    # Se muestra la tabla
    output$tabla1 <- renderTable({
        if(is.null(input$archivo_rdata))
          return(NULL)
        
      # Mensaje de OK
      showModal(modalDialog(
        h3(icon("check-circle", lib = "font-awesome", class = "fa-1x"),
           " El archivo se ha importado correctamente"),
        easyClose = TRUE,
        footer = NULL
      ))
    
      tabla_aux <- as.data.frame(table(labels)) %>% rename(Etiqueta = labels, Frecuencia = Freq)
      return(tabla_aux)
    })
  
  output$sankey <- renderPlot({
    if(is.null(input$archivo_rdata))
      return(NULL)
    
    # Número de casos
    # Train
    #table(labels_train)
    entr_tum <- table(labels_train())[1]
    entr_san <- table(labels_train())[2]
    
    # Test
    #table(labels_test)
    test_tum <- table(labels_test())[1]
    test_san <- table(labels_test())[2]
    # Total
    #table(labels)
    
    # Verificar balanceo de clase en entrenamiento y test
    # Train
    #labels_train %>% table %>% prop.table %>% round(3) * 100
    # Test
    #labels_test %>% table %>% prop.table %>% round(3) * 100
    # Total
    #labels %>% table %>% prop.table %>% round(3) * 100
    
    # Diagrama de Sankey
    datos_sankey <- data.frame(tipo = c(paste0("Tumor\n", entr_tum + test_tum, " casos"), paste0("Tumor\n", entr_tum + test_tum, " casos"),
                                        paste0("Tejido normal\n", entr_san + test_san, " casos"), paste0("Tejido normal\n", entr_san + test_san, " casos")),
                               traintest = c("Entrenamiento", "Test", "Entrenamiento", "Test"),
                               value = c(entr_tum, test_tum, entr_san, test_san))
    
    # Pequeño reorden para que mejorar la presentación de los datos
    datos_sankey$tipo <- factor(datos_sankey$tipo,
                                levels = c(paste0("Tumor\n", entr_tum + test_tum, " casos"), paste0("Tejido normal\n", entr_san + test_san, " casos")),
                                ordered = T)
    
    print(datos_sankey)
    
    ggplot(data = datos_sankey,
           aes(axis1 = tipo, axis2 = traintest, y = value)) +
      scale_x_discrete(limits = c("Tipo de muestra", "Entrenamiento-test"),
                       expand = c(.1, .05)) +
      ylab("") +
      geom_alluvium(col = "black", alpha = 1) +
      geom_alluvium(aes(fill = tipo), alpha = .6, show.legend = FALSE) +
      geom_stratum() +
      geom_text(stat = "stratum", infer.label = TRUE, cex = 3) +
      theme_minimal() +
      ggtitle("Partición en conjuntos de entrenamiento y test",
              paste0("Reparto ", input$porcentaje_entrenamiento, "% - ", 100 - input$porcentaje_entrenamiento, "% con balanceo de clases")) +
      theme(plot.title = element_text(hjust = .5),
            plot.subtitle = element_text(hjust = .5),
            axis.text = element_text(color = "black", margin = margin(t = -30), size = 12),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank()) 
    })
  
  
  # Esto se puede mejorar para que calcule los mejores 50 genes y luego sólo tenga que hacer un subset de esa tabla
  # seleccionando sólo los primeros X elementos, para que reaccione antes la aplicación web.
  
  output$genes_mrmr <- renderTable({
    # Método mRMR (mínima redundancia, máxima relevancia)
    mrmrRanking <- featureSelection(particion.entrenamiento(), labels_train(), colnames(particion.entrenamiento()),
                                    mode = "mrmr")
    
    mrmrRanking <- names(mrmrRanking)[1:input$numero_genes]
    
    return(mrmrRanking)
  }, colnames = FALSE)
  
  output$genes_rf <- renderTable({
    # Método random forest
    rfRanking <- featureSelection(particion.entrenamiento(), labels_train(), colnames(particion.entrenamiento()),
                                  mode = "rf")
    rfRanking <- rfRanking[1:input$numero_genes]
    
    return(rfRanking)
  }, colnames = FALSE)
  
  output$genes_da <- renderTable({
    daRanking <- featureSelection(particion.entrenamiento(), labels_train(), colnames(particion.entrenamiento()),
                                  mode = "da", disease = "liver cancer")
    
    daRanking <- names(daRanking)[1:input$numero_genes]
    
    return(daRanking)
  }, colnames = FALSE)
  
  # Método mRMR (mínima redundancia, máxima relevancia)
  mrmrRanking <- reactive({
    aux <- featureSelection(particion.entrenamiento(), labels_train(), colnames(particion.entrenamiento()),
                                  mode = "mrmr")
    return(names(aux)[1:input$numero_genes])
  })
                          
  
  results_cv <- reactive({
    # Reestructurar para no calcular dos veces mrmrranking!

    svm_CV(particion.entrenamiento(), labels_train(), mrmrRanking(),
           numFold = as.numeric(input$numero_folds))
  })
  
  output$mejores_parametros <- renderTable({
    as.data.frame(results_cv()$bestParameters)
    }, rownames = TRUE)
  
  output$resultados_entrenamiento <- renderPlot({
    
    # Quizá mejor con F1
    dataPlot(results_cv()$accMatrix[, 1:12], colours = rainbow(as.numeric(input$numero_folds)),
             mode = "classResults",
             main = "mRMR - Accuracy for each fold",
             xlab = "Genes",
             ylab = "Accuracy")
    
  })
  
  }) # Cierre botón import
  
  set.seed(122)
  histdata <- rnorm(500)
  
  output$plot1 <- renderPlot({
    data <- histdata[seq_len(input$slider)]
    hist(data)
  })
  
}

shinyApp(ui, server)