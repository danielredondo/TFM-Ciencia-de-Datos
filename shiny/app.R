# Before the deploy to shinyapps.io, run:
# library(BiocManager)
# options(repos = BiocManager::repositories())
# See here:
# https://community.rstudio.com/t/failing-to-deploy-shinyapp-depending-on-bioconductor-packages/6970/5

# Packages
library(shiny)
library(shinydashboard)
library(dplyr)
library(KnowSeq)
library(reshape2)
library(caret)
library(ggplot2)
library(ggalluvial)
library(DT)
library(waiter)

# File to slightly modify dataPlot function
source("www/dataPlot.R")

# Define some spinners
spinner_abrir <- tagList(
  spin_folding_cube(),
  span(br(), h4("Loading application..."), style="color:white;")
)

spinner <- tagList(
  spin_chasing_dots(),
  span(br(), h4("Loading..."), style="color:white; display: inline-block;")
)

ui <- dashboardPage(title = "biomarkeRs", # Title in web browser
  ## Theme
  skin = "black",
  ## Header
  dashboardHeader(title = span(
    "biomarkeRs",
    style = "font-family: Lucida Console; font-weight: bold"
  )),
  ## Sidebar
  dashboardSidebar(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
    ),
    
    sidebarMenu(
      menuItem("Introduction", tabName = "intro", icon = icon("file-alt")),
      menuItem("Data loading", tabName = "datos", icon = icon("database")),
      menuItem("Genes selection", tabName = "genes", icon = icon("dna")),
      menuItem("Model training", tabName = "entrenamiento", icon = icon("play")),
      menuItem("Model validation", tabName = "validation", icon = icon("check-circle")),
      menuItem("Related diseases", tabName = "enfermedades", icon = icon("disease")),
      menuItem("Authors", tabName = "autores", icon = icon("id-card")),
      menuItem("Code", tabName = "codigo", icon = icon("code"))
    )
  ),
  ## Body
  dashboardBody(
    use_waiter(),
    # Spinners to show on load, or when the application is busy
    #waiter_show_on_load(spinner_abrir, color = "#027368"),
    #waiter_on_busy(spinner, color = "#027368"),
    tabItems(
      # Tab 1
      tabItem(tabName = "intro",
              
              h1("About this web application"),
              "This web application allows users with no previous knowledge of programming to analyze transcriptomics data using machine learning.",
              br(), br(),
              
              "The ", tags$i("biomarkeRs"), "application is part of the", 
              tags$a(
                "final master's project of Daniel Redondo Sánchez.",
                href = "https://github.com/danielredondo/TFM_ciencia_de_datos",
                target="_blank"
              ),
              "It´s developed in R-Shiny and the code is ",
              tags$a(
                "open source.",
                href = "https://github.com/danielredondo/TFM_ciencia_de_datos/blob/master/shiny/app.R",
                target="_blank"
              ),
              
              h2("Abstract of the final master's project: Epidemiology and biomarkers detection in cancer"),
              
              h3(tags$b("Introduction")),
              "Cancer is one of the world's largest public health problems with more than 17 million new cases and 9 million deaths every year.",
              h3(tags$b("Methods")),
              "This work focuses on liver cancer and colon-rectum cancer, describing their main epidemiological indicators and using machine learning
              to analyze more than 1,100 RNA-Seq samples from cancer patients. For binary (tumor vs. normal tissue) and multiclass (various tumor 
              types vs. normal tissue) classification, the 10 most relevant genes are identified and predictive models are constructed with SVM, 
              random forest and kNN with 5-fold cross-validation.",
              h3(tags$b("Results")),
              "The best binary classifiers are validated with excellent results for liver cancer (F1-Score in test: 99.5%) and colon-rectum cancer
              (F1-Score: 100%). Lower evaluation measures are obtained in the best models for multiclass classification, in both liver (F1-Score:
              91.8%) and colon-rectum (F1-Score: 79.3%). ",
              br(), br(),
              "A web application has been developed, biomarkeRs, that implements transcriptomic analysis and can be useful for users
              with no previous knowledge of programming.",
              h3(tags$b("Conclusions")),
              "SVM, random forest and kNN obtained very similar results, and managed to correctly distinguish between tumoral and normal tissues
              with some troubles distinguishing between different types of cancer. External validation and clinical interpretations are necessary
              to clearly establish a gene-disease association.",
              

              # Images
              br(), br(), br(),
              fluidRow(column(6, tags$img(src = "ugr.png", height = "100px")),
                       column(6, tags$img(src = "knowseq.png", height = "120px")))
              
      ),
      
      # Tab 2
      tabItem(tabName = "datos",

              # Left column
              fluidRow(column(6, 
                h1("Data loading"),
                fileInput(inputId = "file_labels",
                          label = span("Select CSV file with labels (see ",
                                       tags$a(
                                         "here",
                                         href = "https://raw.githubusercontent.com/danielredondo/TFM_ciencia_de_datos/master/shiny/datos/higado_200genes_labels.csv",
                                         target="_blank"
                                         ),
                                       "an example)"),
                          accept = ".csv",
                          width = "100%"
                ),
                fileInput(inputId = "file_DEGsMatrix",
                          label = span("Select CSV file with DEGsMatrix (see ",
                                       tags$a(
                                         "here",
                                         href = "https://github.com/danielredondo/TFM_ciencia_de_datos/raw/master/shiny/datos/higado_200genes_DEGsMatrix.csv",
                                         target="_blank"
                                       ),
                                       "an example)"),
                          accept = ".csv",
                          width = "100%"
                ),
                
                actionButton(inputId = "boton_importar",
                             label = "Import file",
                             icon = icon("fas fa-file-import", lib = "font-awesome"),
                             width = "100%"),
                br(),
                
                conditionalPanel(condition = "input.boton_importar!=0",
                                 
                                 h2("Distribution of classes"),
                                 
                                 tableOutput("tabla1")
                                 
                                 )),
                 # Right column
                 column(6, br(), br(), br(),
                          conditionalPanel(condition = "input.boton_importar!=0",
                             h2("Train-test partition"),
                             
                             sliderInput("porcentaje_entrenamiento",
                                         label = "Train percentage (%)",
                                         value = 75, min = 5, max = 95, step = 5,
                                         width = "100%"
                             ),
                            h2("Sankey plot"),
                            plotOutput("sankey", width = "100%")))
              )

      ),
      
      # Tab 3
      tabItem(tabName = "genes",
              h1("Genes selection"),
              sliderInput(inputId = "numero_genes", label = "Select the number of genes to use", value = 50, min = 1, max = 50, step = 1, width = "50%"),
              
              textInput(inputId = "disease_da", label = "Disease for DA algorithm", value = "liver cancer", width = "50%"),
              
              actionButton(inputId = "boton_genes",
                           label = "Select most relevant genes",
                           icon = icon("dna", lib = "font-awesome"),
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
              
              # Choose feature selection algorithm
              selectInput("fs_algorithm",
                          label = "Feature selection algorithm",
                          choices = c("mRMR", "RF", "DA"),
                          selected = "mRMR",
                          width = "50%"),
              
              # Choose classification algorithm
              selectInput("cl_algorithm",
                          label = "Classification algorithm",
                          choices = c("SVM", "RF", "kNN"),
                          selected = "SVM",
                          width = "50%"),
              
              # Choose number of folds
              selectInput("number_folds",
                          label = "Number of folds",
                          choices = c(3, 5, 10),
                          selected = 5,
                          width = "50%"),
              
              # Train model button
              actionButton(inputId = "boton_model_training",
                           label 
                           = "Train model",
                           icon = icon("play", lib = "font-awesome"),
                           width = "50%"),
              
              br(),
              br(),
              
              # Show optimal parameters (if the classification method is SVM or kNN)
              conditionalPanel(condition = "input.cl_algorithm == 'SVM'",
                               br(),
                               textOutput("optimal_svm"),
                               br()),
              conditionalPanel(condition = "input.cl_algorithm == 'kNN'",
                               br(),
                               textOutput("optimal_knn"),
                               br()),
                               
              dataTableOutput("results_cv")
              
      ),

      # Tab 5
      tabItem(tabName = "validation",
              h1("Model validation"),
              
              selectInput("fs_algorithm_validation",
                          label = "Feature selection algorithm:",
                          choices = c("mRMR", "RF", "DA"),
                          selected = "mRMR",
                          width = "50%"),
              
              selectInput("cl_algorithm_validation",
                          label = "Classification algorithm (for SVM and kNN it must be trained first to obtain optimal parameters):",
                          choices = c("SVM", "RF", "kNN"),
                          selected = "SVM",
                          width = "50%"),
              
              sliderInput(inputId = "numero_genes_validation", label = "Select the number of genes to use (must be equal or less than the number of genes selected at 'Genes selection'):",
                          value = 10, min = 1, max = 50, step = 1, width = "50%"),

              actionButton(inputId = "boton_model_validation",
                           label = "Validate model in test",
                           icon = icon("play", lib = "font-awesome"),
                           width = "50%"),
              
              br(),
              br(),
              
              plotOutput("results_validation",
                         width = "50%")
              
      ),
      
      # Tab 6
      tabItem(tabName = "enfermedades",
            h1("Related diseases"),
            textInput(inputId = "gene_for_disease", label = "Gene", value = "TERT", width = "50%"),
            
            dataTableOutput("gene_for_disease_table")
      ),
      
      # Tab 7
      tabItem(tabName = "autores",
              h1("Authors"),
              tags$h4(
                tags$li("Daniel Redondo Sánchez. Granada Cancer Registry, ibs.GRANADA."), br(),
                tags$li("Daniel Castillo. University of Granada."), br(),
                tags$li("Luis Javier Herrera. University of Granada.")),
                
                h2("Contact"),
                "Daniel Redondo Sánchez (daniel.redondo.easp at juntadeandalucia.es) is the main developer of this Shiny App.", br(),
                "Daniel Castillo (cased at ugr.es) is the main developer of the R package KnowSeq."
              
              ),
      # Tab 8
      tabItem(tabName = "codigo",
              h1("Code"),
              tags$h4(
                "In ", tags$a(href = "https://github.com/danielredondo/TFM_ciencia_de_datos/tree/master/shiny", "this GitHub repository"),
                "you can find the code of the web application.")
              )
      ) # Close tabs
  ) # Close dashboard body
) # Close dashboard page

# Extend size of accepted files (40MB instead of the 5MB - default)
options(shiny.maxRequestSize = 40*1024^2)

server <- function(input, output){

  values <- reactiveValues(ranking = NULL, optimalSVM_train = NULL, optimalkNN_train = NULL)
  
  # Server of tab: Data loading ------
  
  observeEvent(input$boton_importar, {
    
    # If files are selected, they are imported
    # Read labels
    labels <- as.vector(t(read.csv2(file = input$file_labels$datapath)))
    # Read DEGsMatrix
    DEGsMatrix <- as.data.frame(read.csv2(file = input$file_DEGsMatrix$datapath, row.names = 1))
    filas <- rownames(DEGsMatrix)
    DEGsMatrix <- apply(DEGsMatrix, 2, as.numeric)
    rownames(DEGsMatrix) <- filas
    # Create DEGsMatrixML (for machine learning purposes)
    DEGsMatrixML <- t(DEGsMatrix)
    
    # Train-test partition
    set.seed(31415)
    indices <- reactive(createDataPartition(labels, p = input$porcentaje_entrenamiento / 100, list = FALSE))
    particion <- reactive(list(training = DEGsMatrixML[indices(), ], test = DEGsMatrixML[-indices(), ]))
    
    particion.entrenamiento <- reactive(particion()$training)
    particion.test <- reactive(particion()$test)
    
    # Labels
    labels_train <- reactive(labels[indices()])
    labels_test  <- reactive(labels[-indices()])
    
    # Table
    output$tabla1 <- renderTable({
        if(is.null(input$file_labels)) return(NULL)
        
      # Message if file is correctly imported
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
    if(is.null(input$file_labels)) return(NULL)
    
    # Train
    #table(labels_train)
    entr_tum <- table(labels_train())[1]
    entr_san <- table(labels_train())[2]
    
    # Test
    #table(labels_test)
    test_tum <- table(labels_test())[1]
    test_san <- table(labels_test())[2]

    # Sankey diagram
    datos_sankey <- data.frame(tipo = c(paste0("Tumour\n", entr_tum + test_tum, " cases"), paste0("Tumour\n", entr_tum + test_tum, " cases"),
                                        paste0("Normal tissue\n", entr_san + test_san, " cases"), paste0("Normal tissue\n", entr_san + test_san, " cases")),
                               traintest = c(paste0("Train\n", entr_tum, " tumour\n", entr_san, " normal tissue"),
                                             paste0("Test\n", test_tum, " tumour\n", test_san, " normal tissue"),
                                             paste0("Train\n", entr_tum, " tumour\n", entr_san, " normal tissue"),
                                             paste0("Test\n", test_tum, " tumour\n", test_san, " normal tissue")),
                               value = c(entr_tum, test_tum, entr_san, test_san))
    
    # Reordering types
    datos_sankey$tipo <- factor(datos_sankey$tipo,
                                levels = c(paste0("Tumour\n", entr_tum + test_tum, " cases"), paste0("Normal tissue\n", entr_san + test_san, " cases")),
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
  }) # Close import button
  
  
  # Server of tab: Genes selection ------
  
  w <- Waiter$new(html = span(""))
  
  observeEvent(input$boton_genes, {
    
    w$show()
    
    labels <- as.vector(t(read.csv2(file = input$file_labels$datapath)))
    DEGsMatrix <- as.data.frame(read.csv2(file = input$file_DEGsMatrix$datapath, row.names = 1))
    filas <- rownames(DEGsMatrix)
    DEGsMatrix <- apply(DEGsMatrix, 2, as.numeric)
    rownames(DEGsMatrix) <- filas
    DEGsMatrixML <- t(DEGsMatrix)
    
    set.seed(31415)
    indices <- reactive(createDataPartition(labels, p = input$porcentaje_entrenamiento / 100, list = FALSE))
    particion <- reactive(list(training = DEGsMatrixML[indices(), ], test = DEGsMatrixML[-indices(), ]))
    
    particion.entrenamiento <- reactive(particion()$training)
    particion.test <- reactive(particion()$test)
    
    labels_train <- reactive(labels[indices()])
    labels_test  <- reactive(labels[-indices()])
    w$hide()

    # mRMR method
    w <- Waiter$new(html = tagList(spin_folding_cube(),
                                   span(br(), br(), br(), h4("Running mRMR algorithm..."),
                                        style="color:white;")))
    w$show()
    mrmrRanking <- featureSelection(particion.entrenamiento(), labels_train(), colnames(particion.entrenamiento()),
                                    mode = "mrmr")
    mrmrRanking <- names(mrmrRanking)
    w$hide()
    
    # RF method
    w <- Waiter$new(html = tagList(spin_folding_cube(),
                                   span(br(), br(), br(), h4("Running RF algorithm..."),
                                        style="color:white;")))
    w$show()
    rfRanking <- featureSelection(particion.entrenamiento(), labels_train(), colnames(particion.entrenamiento()),
                                  mode = "rf")
    w$hide()
    
    # DA method
    w <- Waiter$new(html = tagList(spin_folding_cube(),
                                   span(br(), br(), br(), h4("Running DA algorithm..."),
                                        style="color:white;")))
    w$show()
    daRanking <- NULL
    
    daRanking <- featureSelection(particion.entrenamiento(), labels_train(), colnames(particion.entrenamiento()),
                                  mode = "da", disease = input$disease_da)
    daRanking <- names(daRanking)
    
    # If there has been any problem in the API call, it's repeated after a 3 second break
    while(is.null(daRanking)){
      Sys.sleep(3)
      daRanking <- featureSelection(particion.entrenamiento(), labels_train(), colnames(particion.entrenamiento()),
                                    mode = "da", disease = input$disease_da)
      daRanking <- names(daRanking)
    }
  
    w$hide()
    
    
  values$ranking <- cbind(mrmrRanking, rfRanking, daRanking)
    
  # Ranking tables
    
  output$genes_mrmr <- renderTable({
    mrmrRanking <- mrmrRanking[1:input$numero_genes]
    return(mrmrRanking)
  }, colnames = FALSE)
  
  output$genes_rf <- renderTable({
    rfRanking <- rfRanking[1:input$numero_genes]
    return(rfRanking)
  }, colnames = FALSE)
  
  output$genes_da <- renderTable({
    daRanking <- daRanking[1:input$numero_genes]
    return(daRanking)
  }, colnames = FALSE)

  }) # Close button

  # Server of tab: Model training ------
  
  w2 <- Waiter$new(html = tagList(spin_folding_cube(),
                                 span(br(), br(), br(), h4(""),
                                 style="color:white;")))  
  
  observeEvent(input$boton_model_training, {
    
    w2$show()
    

    labels <- as.vector(t(read.csv2(file = input$file_labels$datapath)))
    DEGsMatrix <- as.data.frame(read.csv2(file = input$file_DEGsMatrix$datapath, row.names = 1))
    filas <- rownames(DEGsMatrix)
    DEGsMatrix <- apply(DEGsMatrix, 2, as.numeric)
    rownames(DEGsMatrix) <- filas
    DEGsMatrixML <- t(DEGsMatrix)
    
    set.seed(31415)
    indices <- reactive(createDataPartition(labels, p = input$porcentaje_entrenamiento / 100, list = FALSE))
    particion <- reactive(list(training = DEGsMatrixML[indices(), ], test = DEGsMatrixML[-indices(), ]))
    
    particion.entrenamiento <- reactive(particion()$training)
    particion.test <- reactive(particion()$test)
    
    labels_train <- reactive(labels[indices()])
    labels_test  <- reactive(labels[-indices()])
    w2$hide()
    
    if(input$fs_algorithm == "mRMR"){
      ranking <- values$ranking[1:input$numero_genes, 1]
    }
    
    if(input$fs_algorithm == "RF"){
      ranking <- values$ranking[1:input$numero_genes, 2]
    }
    
    if(input$fs_algorithm == "DA"){
      ranking <- values$ranking[1:input$numero_genes, 3]
    }
    
    if(input$cl_algorithm == "SVM"){
      w3 <- Waiter$new(html = tagList(spin_folding_cube(),
                                      span(br(), br(), br(), h4("Training SVM algorithm..."),
                                           style="color:white;")))  
      w3$show()
      results_cv <- svm_trn(particion.entrenamiento(), labels_train(), ranking,
                            numFold = as.numeric(input$number_folds))
      values$optimalSVM_train <- results_cv$bestParameters
      w3$hide()
    }
    
    if(input$cl_algorithm == "RF"){
      w3 <- Waiter$new(html = tagList(spin_folding_cube(),
                                      span(br(), br(), br(), h4("Training RF algorithm..."),
                                      style="color:white;")))  
      w3$show()
      results_cv <- rf_trn(particion.entrenamiento(), labels_train(), ranking,
                                numFold = as.numeric(input$number_folds))
      w3$hide()
    }
    
    if(input$cl_algorithm == "kNN"){
      w3 <- Waiter$new(html = tagList(spin_folding_cube(),
                                      span(br(), br(), br(), h4("Training kNN algorithm..."),
                                           style="color:white;")))  
      w3$show()
      results_cv <- knn_trn(particion.entrenamiento(), labels_train(), ranking,
                                numFold = as.numeric(input$number_folds))
      values$optimalkNN_train <- results_cv$bestK
      
      w3$hide()
    }
    
    output$optimal_svm <- renderText(paste0("\nOptimal coefficients for ", input$numero_genes, " genes : cost = ", results_cv$bestParameters[1], "; gamma = ", results_cv$bestParameters[2]))
    
    output$optimal_knn <- renderText(paste0("\nOptimal number of neighbours for ", input$numero_genes, " genes = ", results_cv$bestK))
    
    output$results_cv <- renderDataTable({
      df <- data.frame(`Number of genes` = 1:length(ranking),
               Accuracy = round(100 * results_cv$accuracyInfo$meanAccuracy, 2),
               `F1-Score` = round(100 * results_cv$F1Info$meanF1, 2), 
               Sensitivity = round(100 * results_cv$sensitivityInfo$meanSensitivity, 2),
               Specificity = round(100 * results_cv$specificityInfo$meanSpecificity, 2),
               check.names = F)
      dat <- datatable(df, rownames = F, options = list(pageLength = 10)) %>% formatStyle(names(df)[2:5],
                               background = styleColorBar(range(df[, 2:5]) - c(1, 0), "forestgreen"),
                               backgroundSize = "98% 88%",
                               backgroundRepeat = "no-repeat",
                               backgroundPosition = "center")
      return(dat)}
    )
    
  }) 
  
  
  # Server of tab: Model validation  ------
  
  w3 <- Waiter$new(html = tagList(spin_folding_cube(),
                                  span(br(), br(), br(), h4("Validating model..."),
                                  style="color:white;")))  
  
  observeEvent(input$boton_model_validation, {
    
    w3$show()
    
    labels <- as.vector(t(read.csv2(file = input$file_labels$datapath)))
    DEGsMatrix <- as.data.frame(read.csv2(file = input$file_DEGsMatrix$datapath, row.names = 1))
    filas <- rownames(DEGsMatrix)
    DEGsMatrix <- apply(DEGsMatrix, 2, as.numeric)
    rownames(DEGsMatrix) <- filas
    DEGsMatrixML <- t(DEGsMatrix)
    
    set.seed(31415)
    indices <- reactive(createDataPartition(labels, p = input$porcentaje_entrenamiento / 100, list = FALSE))
    particion <- reactive(list(training = DEGsMatrixML[indices(), ], test = DEGsMatrixML[-indices(), ]))
    
    particion.entrenamiento <- reactive(particion()$training)
    particion.test <- reactive(particion()$test)
    
    labels_train <- reactive(labels[indices()])
    labels_test  <- reactive(labels[-indices()])
    w3$hide()
    
    if(input$fs_algorithm_validation == "mRMR"){
      ranking <- values$ranking[1:input$numero_genes, 1]
    }
    
    if(input$fs_algorithm_validation == "RF"){
      ranking <- values$ranking[1:input$numero_genes, 2]
    }
    
    if(input$fs_algorithm_validation == "DA"){
      ranking <- values$ranking[1:input$numero_genes, 3]
    }
    
    if(input$cl_algorithm_validation == "SVM"){
      w3 <- Waiter$new(html = tagList(spin_folding_cube(),
                                      span(br(), br(), br(), h4("Validating SVM algorithm..."),
                                           style="color:white;")))  
      w3$show()
      results_validation <- svm_test(train = particion.entrenamiento(), labels_train(),
                             test = particion.test(), labels_test(),
                             ranking,
                             bestParameters = values$optimalSVM_train)
      w3$hide()
    }
    
    if(input$cl_algorithm_validation == "RF"){
      w3 <- Waiter$new(html = tagList(spin_folding_cube(),
                                      span(br(), br(), br(), h4("Validating RF algorithm..."),
                                           style="color:white;")))  
      w3$show()
      results_validation <- rf_test(train = particion.entrenamiento(), labels_train(),
                            test = particion.test(), labels_test(),
                            ranking)
      w3$hide()
    }
    
    if(input$cl_algorithm_validation == "kNN"){
      w3 <- Waiter$new(html = tagList(spin_folding_cube(),
                                      span(br(), br(), br(), h4("Validating kNN algorithm..."),
                                           style="color:white;")))  
      w3$show()
      results_validation <- knn_test(train = particion.entrenamiento(), labels_train(),
                                     test = particion.test(), labels_test(),
                                     ranking, bestK = values$optimalkNN_train)
      w3$hide()
    }
    
    output$results_validation <- renderPlot({
      tabla <- results_validation$cfMats[[input$numero_genes_validation]]$table
      plotConfMatrix(tabla)
    })
    
  })
  

  # Server of tab: Related diseases ------
  output$gene_for_disease_table <- renderDataTable(
    {dis <- as.data.frame(DEGsToDiseases(input$gene_for_disease, size = 10000))

    # Round coefficients
    for(i in 2:9){
      dis[, i] <- round(as.numeric(dis[, i]), 2)
    }
    
    names(dis) <- c("Disease", "Overall score", "Literature", "RNA Expr.", "Genetic", "Somatic Mut.", "Drug", "Animal", "Pathways")
    return(dis)}
    , filter = "top", options = list(pageLength = 10)
  )
  
}

shinyApp(ui, server)