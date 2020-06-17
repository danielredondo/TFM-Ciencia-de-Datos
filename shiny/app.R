library(shiny)
library(shinydashboard)
library(dplyr)
library(KnowSeq)
library(reshape2)

source("www/dataPlot.R")

ui <- dashboardPage(
  ## Tema
  skin = "black",
  ## Cabecera
  dashboardHeader(title = "TFM: Biomarcadores en cáncer",
                  titleWidth = 350),
  ## Barra lateral
  dashboardSidebar(
    sidebarMenu(
      menuItem("Introducción", tabName = "intro", icon = icon("file-alt")),
      menuItem("Carga de datos", tabName = "datos", icon = icon("database")),
      menuItem("Selección de genes", tabName = "genes", icon = icon("dna")),
      menuItem("Gráficos", tabName = "graficos", icon = icon("chart-bar")),
      menuItem("Autor", tabName = "autor", icon = icon("id-card"))
      
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
              tags$img(src = "ugr.png", width = "250px")
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
              
              tableOutput("tabla1")
      ),
      
      # Tab 3
      tabItem(tabName = "genes",
              h2("Selección de genes"),
              
              br(),
              
              h2("Enfermedades relacionadas")
      ),
      
      # Tab 4
      tabItem(tabName = "graficos",
              h2("Gráficos"),
              
              fluidRow(
                column(6, plotOutput("plot1")),
                column(6, 
                       title = "Controls",
                       sliderInput("slider", "Number of observations:",
                                   min = 0,
                                   max = 100,
                                   value = 50,
                                   step = 1)
                )
              )
      ),
      # Tab 5
      tabItem(tabName = "autor",
              h2("Autor")
        )
      ) # Final tabs
  ) # Final dashboard body
) # Final dashboard page

# Ampliar tamaño de .RData a 15MB en lugar de los 5MB por defecto
options(shiny.maxRequestSize = 15*1024^2)

server <- function(input, output){
  
  observeEvent(input$boton_importar, {
    
    # Se muestra la tabla
    output$tabla1 <- renderTable({
        if (is.null(input$archivo_rdata))
          return(NULL)
        
      # Mensaje de OK
      showModal(modalDialog(
        h3(icon("check-circle", lib = "font-awesome", class = "fa-1x"),
           " El archivo se ha importado correctamente"),
        easyClose = TRUE,
        footer = NULL
      ))
      
      # Si se ha seleccionado un fichero, se importa
      load(input$archivo_rdata$datapath)
      labels <- matriz[1, ] %>% as.vector
      tabla_aux <- as.data.frame(table(labels)) %>% rename(Etiqueta = labels, Frecuencia = Freq)
      return(tabla_aux)
    })
  })
  
  set.seed(122)
  histdata <- rnorm(500)
  
  output$plot1 <- renderPlot({
    data <- histdata[seq_len(input$slider)]
    hist(data)
  })
  
}

shinyApp(ui, server)