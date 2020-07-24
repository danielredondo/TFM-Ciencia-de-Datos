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
library(waiter)

source("www/dataPlot.R")

spinner_abrir <- tagList(
  spin_folding_cube(),
  span(br(), h4("Loading application..."), style="color:white;")
)

spinner <- tagList(
  spin_chasing_dots(),
  span(br(), h4("Loading..."), style="color:white; display: inline-block;")
)

ui <- dashboardPage(title = "biomarkeRs", # Title in web browser
  ## Tema
  skin = "black",
  ## Cabecera
  dashboardHeader(title = span(
    "biomarkeRs",
    style = "font-family: Lucida Console; font-weight: bold"
  )),
  ## Barra lateral
  dashboardSidebar(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
    ),
    
    sidebarMenu(
      menuItem("Introduction", tabName = "intro", icon = icon("file-alt")),
      menuItem("Data loading", tabName = "datos", icon = icon("database")),
      menuItem("Genes selection", tabName = "genes", icon = icon("dna")),
      menuItem("Model training", tabName = "entrenamiento", icon = icon("play")),
      menuItem("Model validation", tabName = "validacion", icon = icon("check-circle")),
      menuItem("Related diseases", tabName = "enfermedades", icon = icon("disease")),
      menuItem("Authors", tabName = "autores", icon = icon("id-card")),
      menuItem("Code", tabName = "codigo", icon = icon("code"))
    )
  ),
  ## Cuerpo
  dashboardBody(
    use_waiter(),
    #waiter_show_on_load(spinner_abrir, color = "#027368"),
    #waiter_on_busy(spinner, color = "#027368"),
    tabItems(
      # Tab 1
      tabItem(tabName = "intro",
              h1("Epidemiology and biomarkers detection in cancer"),
              
              h3(tags$b("Introduction")),
              "[Text]",
              h3(tags$b("Methods")),
              "[Text]",
              h3(tags$b("Results")),
              "[Text]",
              h3(tags$b("Conclusions")),
              "[Text]",
              
              h2("About this web application"),
              "[Text]",
              # Parte final
              br(), br(), br(),
              fluidRow(column(6, tags$img(src = "ugr.png", height = "100px")),
                       column(6, tags$img(src = "knowseq.png", height = "120px")))
              
      ),
      
      # Tab 2
      tabItem(tabName = "datos",

              h1("Data loading"),
              fileInput(inputId = "archivo_rdata",
                        label = "Select .RData file. (see Example File)",
                        accept = ".RData",
                        width = "50%"
              ),
              
              actionButton(inputId = "boton_importar",
                           label = "Import file",
                           icon = icon("fas fa-file-import", lib = "font-awesome"),
                           width = "50%"),
              br(),
              
              conditionalPanel(condition = "input.boton_importar!=0",
                tableOutput("tabla1"),
                
                h2("Train-test partition"),
                
                sliderInput("porcentaje_entrenamiento",
                            label = "Train percentage (%)",
                            value = 75, min = 5, max = 95, step = 5,
                            width = "50%"
                            ),
                
                plotOutput("sankey", width = "50%")
              )
      ),
      
      # Tab 3
      tabItem(tabName = "genes",
              h1("Genes selection"),
              sliderInput(inputId = "numero_genes", label = "Select the number of genes to use", value = 20, min = 1, max = 51, step = 1),
              
              textInput(inputId = "disease_da", label = "Disease for DA algorithm", value = "liver cancer", width = "50%"),
              
              actionButton(inputId = "boton_genes",
                           label = "Select most relevant genes",
                           icon = icon("fas fa-calculator", lib = "font-awesome"),
                           width = "50%"),
              br(),
              br(),
              conditionalPanel(condition = "input.boton_genes!=0",
                                 
                h3("Table of more relevant genes by feature selection method:"),
                br(),
                fluidRow(
                  column(4, h4(tags$b("  MRMR")), tableOutput("genes_mrmr")),
                  column(4, h4(tags$b("  RF")), tableOutput("genes_rf")),
                  column(4, h4(tags$b("  DA")), tableOutput("genes_da")),
                )
              )
      ),
      
      # Tab 4
      tabItem(tabName = "entrenamiento",
              h1("Model training"),
              
              # Elegir MRMR/RF/DA, se muestra MRMR por ahora
              selectInput("tipo_entrenamiento",
                          label = "Feature selection algorithm",
                          choices = c("mRMR", "rf", "da"),
                          selected = "mRMR",
                          width = "50%"),
              
              # Se puede añadir número de CV que se hacen
              selectInput("numero_folds",
                          label = "Number of folds",
                          choices = c(3, 5, 10),
                          selected = 5,
                          width = "50%"),
              
              tableOutput("mejores_parametros"),
              
              plotOutput("resultados_entrenamiento")
              
              # Falta añadir el F1-score y decir cuál es el mejor método
      ),

      
      # Tab 6
      tabItem(tabName = "enfermedades",
            h1("Related diseases"),
            "Aquí irá un campo para poner un gen, y se devolverán los resultados de KnowSeq::DEGsToDiseases"
      ),
      
      # Tab 7
      tabItem(tabName = "autores",
              h1("Authors"),
              tags$h4(
                tags$li(tags$b("Daniel Redondo Sánchez")), br(),
                tags$li(tags$b("Ignacio Rojas")), br(),
                tags$li(tags$b("Luis Javier Herrera")), br(),
                tags$li(tags$b("Daniel Castillo"))
                )
              ),
      # Tab 8
      tabItem(tabName = "codigo",
              h1("Code"),
              tags$h4(
                "In ", tags$a(href = "https://github.com/danielredondo/TFM_ciencia_de_datos/blob/master/shiny/app.R", "this GitHub repository"),
                "you can find the code of the web application.")
              )
      ) # Final tabs
  ) # Final dashboard body
) # Final dashboard page

# Ampliar tamaño de .RData a 15MB en lugar de los 5MB por defecto
options(shiny.maxRequestSize = 15*1024^2)

server <- function(input, output){

  #Sys.sleep(2)
  # give time for wait screen to show
  #waiter_hide()
  
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
    set.seed(31415)
    indices <- reactive(createDataPartition(labels, p = input$porcentaje_entrenamiento / 100, list = FALSE))
    particion <- reactive(list(training = DEGsMatrixML[indices(), ], test = DEGsMatrixML[-indices(), ]))
    
    # Conjuntos
    particion.entrenamiento <- reactive(particion()$training)
    particion.test <- reactive(particion()$test)
    
    # Labels
    labels_train <- reactive(labels[indices()])
    labels_test  <- reactive(labels[-indices()])
    
    # Se muestra la tabla
    output$tabla1 <- renderTable({
        if(is.null(input$archivo_rdata))
          return(NULL)
        
      # Mensaje de OK
      showModal(modalDialog(
        h3(icon("check-circle", lib = "font-awesome", class = "fa-1x"),
           " File imported"),
        easyClose = TRUE,
        footer = NULL
      ))
    
      tabla_aux <- as.data.frame(table(labels)) %>% rename(Label = labels, Samples = Freq)
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

    # Diagrama de Sankey
    datos_sankey <- data.frame(tipo = c(paste0("Tumor\n", entr_tum + test_tum, " casos"), paste0("Tumor\n", entr_tum + test_tum, " casos"),
                                        paste0("Normal tissue\n", entr_san + test_san, " casos"), paste0("Normal tissue\n", entr_san + test_san, " casos")),
                               traintest = c("Train", "Test", "Train", "Test"),
                               value = c(entr_tum, test_tum, entr_san, test_san))
    
    # Pequeño reorden para que mejorar la presentación de los datos
    datos_sankey$tipo <- factor(datos_sankey$tipo,
                                levels = c(paste0("Tumor\n", entr_tum + test_tum, " casos"), paste0("Normal tissue\n", entr_san + test_san, " casos")),
                                ordered = T)
    
    ggplot(data = datos_sankey,
           aes(axis1 = tipo, axis2 = traintest, y = value, label = after_stat(stratum))) +
      scale_x_discrete(limits = c("Type of sample", "Train-test"),
                       expand = c(.1, .05)) +
      ylab("") +
      geom_alluvium(col = "black", alpha = 1) +
      geom_alluvium(aes(fill = tipo), alpha = .6, show.legend = FALSE) +
      geom_stratum() +
      geom_text(stat = "stratum", cex = 3) +
      theme_minimal() +
      ggtitle("Train-test partition") +
      theme(plot.title = element_text(hjust = .5),
            axis.text = element_text(color = "black", margin = margin(t = -30), size = 12),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            panel.grid = element_blank()) 
    })
  }) # Cierre botón import
  
  
  # Esto se puede mejorar para que calcule los mejores 50 genes y luego sólo tenga que hacer un subset de esa tabla
  # seleccionando sólo los primeros X elementos, para que reaccione antes la aplicación web.
  
  w <- Waiter$new(html = span(""))
  
  observeEvent(input$boton_genes, {
    
    w$show()
    
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
    
    # Partición 75% / 25% con balanceo de clase
    set.seed(31415)
    indices <- reactive(createDataPartition(labels, p = input$porcentaje_entrenamiento / 100, list = FALSE))
    particion <- reactive(list(training = DEGsMatrixML[indices(), ], test = DEGsMatrixML[-indices(), ]))
    
    # Conjuntos
    particion.entrenamiento <- reactive(particion()$training)
    particion.test <- reactive(particion()$test)
    
    # Labels
    labels_train <- reactive(labels[indices()])
    labels_test  <- reactive(labels[-indices()])
    w$hide()

    # Método mRMR (mínima redundancia, máxima relevancia)
    w <- Waiter$new(html = tagList(spin_loaders(39, color = "white", style = "scale: 4"),
                                   span(br(), br(), br(), h4("Running mRMR algorithm..."),
                                        style="color:white;")))
    w$show()
    mrmrRanking <- featureSelection(particion.entrenamiento(), labels_train(), colnames(particion.entrenamiento()),
                                    mode = "mrmr")
    w$hide()
    
    # Método random forest
    w <- Waiter$new(html = tagList(spin_loaders(39, color = "white", style = "scale: 4"),
                                   span(br(), br(), br(), h4("Running RF algorithm..."),
                                        style="color:white;")))
    w$show()
    rfRanking <- featureSelection(particion.entrenamiento(), labels_train(), colnames(particion.entrenamiento()),
                                  mode = "rf")
    w$hide()
    
    # Método DA
    w <- Waiter$new(html = tagList(spin_loaders(39, color = "white", style = "scale: 4"),
                                   span(br(), br(), br(), h4("Running DA algorithm..."),
                                        style="color:white;")))
    w$show()
    daRanking <- featureSelection(particion.entrenamiento(), labels_train(), colnames(particion.entrenamiento()),
                                  mode = "da", disease = input$disease_da)
    w$hide()
    
  output$genes_mrmr <- renderTable({
    mrmrRanking <- names(mrmrRanking)[1:input$numero_genes]
    return(mrmrRanking)
  }, colnames = FALSE)
  
  output$genes_rf <- renderTable({
    rfRanking <- rfRanking[1:input$numero_genes]
    return(rfRanking)
  }, colnames = FALSE)
  
  output$genes_da <- renderTable({
    daRanking <- names(daRanking)[1:input$numero_genes]
    return(daRanking)
  }, colnames = FALSE)

  }) # Cierre botón calcular genes

  
  # Leer
  # https://stackoverflow.com/questions/33671915/r-shiny-server-how-to-keep-variable-value-in-observeevent-function
  # https://groups.google.com/g/shiny-discuss/c/UVL3uENQ88k
  # https://shiny.rstudio.com/articles/action-buttons.html
  # para ejecutar sólo una vez particion.entrenamiento y demás funciones.
  
  
  # Método mRMR (mínima redundancia, máxima relevancia)
  mrmrRanking <- reactive({
    aux <- featureSelection(particion.entrenamiento(), labels_train(), colnames(particion.entrenamiento()),
                                  mode = "mrmr")
    return(names(aux)[1:input$numero_genes])
  })
                          
  
  results_trn <- reactive({
    # Reestructurar para no calcular dos veces mrmrranking!

    svm_trn(particion.entrenamiento(), labels_train(), mrmrRanking(),
           numFold = as.numeric(input$numero_folds))
  })
  
  output$mejores_parametros <- renderTable({
    as.data.frame(results_trn()$bestParameters)
    }, rownames = TRUE)
  
  output$resultados_entrenamiento <- renderPlot({
    
    # Quizá mejor con F1
    dataPlot(results_trn()$accMatrix[, 1:12], colours = rainbow(as.numeric(input$numero_folds)),
             mode = "classResults",
             main = "mRMR - Accuracy for each fold",
             xlab = "Genes",
             ylab = "Accuracy")
    
  })
  
}

shinyApp(ui, server)