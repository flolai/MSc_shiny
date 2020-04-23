
library(shiny)
library(dplyr)
library(tidyr)
library(tidyverse)
library(plotly)
library(ropls)
library(ggiraph)
library(shinydashboard)

# ui object
ui <- fluidPage(
  titlePanel("Intergration of metabolomics and transcriptomics data"),
  
  tabsetPanel(
    tabPanel("Data Input",
             tags$style(type="text/css",
                        ".shiny-output-error { visibility: hidden; }",
                        ".shiny-output-error:before { visibility: hidden; }"
             ),  
             
             
             fileInput(inputId = "metadata",
                       label = "metabolomics data in csv format",
                       multiple = TRUE,
                       accept = c(".csv")),
             
             fileInput(inputId = "ngs",
                       label = "transcriptomics data in csv format",
                       multiple = TRUE,
                       accept = c(".csv")),
             
             fileInput(inputId = "comp_name",
                       label = "Compound name in csv format",
                       multiple = TRUE,
                       accept = c(".csv")),
             
             fileInput(inputId = "gene_name",
                       label = "Gene name in csv format",
                       multiple = TRUE,
                       accept = c(".csv")),
    ),
    
    tabPanel("Plots",
             fileInput(inputId = "samp_detail",
                       label = "Sample details in csv format",
                       multiple = TRUE,
                       accept = c(".csv")),
             #adding selection for sample discription classification for colouring
             
             selectizeInput('samp_colour', 'Variable selection for sample details',
                            choices = NULL,
                            options = list(placeholder = 'Select from below'),
                            multiple = FALSE),
             
             
             sidebarPanel(
               
               #DT::dataTableOutput("autoscale_no_norm"),
               #DT::dataTableOutput("autoscale_log2_full"),
               
               
               
               tags$hr(),
               actionButton("do", "Click to merge omics data - no log2 normalisation"),
               tags$hr(),
               actionButton("do_norm", "Click to merge omics data - log2 normalisation"),
               tags$hr(),
               actionButton("autoscale_no_norm", "Click to autoscaling merged data set - no normalisation"),
               tags$hr(),
               actionButton("autoscale_log2", "Click to autoscaling log2 normalised merged data set")
               
             ),
             mainPanel(
               actionButton("pca_full_scores", "PCA for autoscaled data"),
               actionButton("PCA_log2_scores", "PCA for log2, autoscaled data"),
               
               fluidRow(
                 splitLayout(cellWidths = c("50%", "50%"), plotOutput("pca_full_nice"), girafeOutput("pca_full_loadings"))
               ),
               
               fluidRow(
                 splitLayout(cellWidths = c("50%", "50%"), plotOutput("pca_log2_nice"), girafeOutput("pca_log2_loadings"))
               ),
               
               
               actionButton("oplsda_no_norm", "OPLS-DA model for not normalised data"),
               actionButton("oplsda_norm", "OPLS-DA model for normalised data"),
               
               fluidRow(
                 splitLayout(cellWidths = c("50%", "50%"),plotOutput("oplsda_no_norm_out"))
               ),
               
               fluidRow(
                 splitLayout(cellWidths = c("50%", "50%"),plotOutput("oplsda_norm_out"))
               ),
             ),
             actionButton("exp_neg_mk_no_norm", "negative loading markers no normalisation"),
             DT::dataTableOutput("oplsda_no_norm_load"),
    )
  )
)






# server()
server <- function(input, output,session){
  options(shiny.maxRequestSize=30*1024^2)  
  
  
  #reading metabolomics data, 1st column = sample name  
  meta <- reactive({
    if (is.null(input$metadata)){
      return(NULL)
    }
    read.csv(input$metadata$datapath, header = TRUE, row.names = 1)
  })
  
  #reading transcriptomics data, 1st column = sample name    
  ngs <- reactive({
    if (is.null(input$ngs)){
      return(NULL)
    }
    read.csv(input$ngs$datapath, header = TRUE, row.names = 1) 
  })
  
  #reading sample details for classification    
  samp_detail <- reactive({
    if (is.null(input$samp_detail)){
      return(NULL)
    }
    read.csv(input$samp_detail$datapath, header = TRUE, row.names = 1) 
  })
  
  
  observe({
    req(samp_detail())
    updateSelectizeInput(session, 
                         'samp_colour',
                         choices = as.character(names(samp_detail())),
                         selected = NULL)
  })
  
  samp_colour <- reactive({samp_detail()[,input$samp_colour]})
  
  
  #reading compounds and gene details for naming export  
  
  merge_full<- reactiveValues() 
  
  observeEvent(input$do, {
    merge_full$df <- transform(merge(meta(), ngs(), by.x = 0, by.y = 0), row.names = Row.names, Row.names = NULL)# merging metabolomics
    #and ngs data by row names = sample name
  })
  
  log2_full <- reactiveValues()
  
  observeEvent(input$do_norm, {
    log2_full$df <- transform(merge(log2(meta()), log2(ngs()), by.x = 0, by.y = 0, row.names = Row.names, Row.names = NULL)) #log2 normalisaion
    
  })
  
  
  autoscale_full <- reactiveValues()
  
  observeEvent(input$autoscale_no_norm, {
    merge_full$df[is.na(merge_full$df)] <- 0.001
    autoscale_full$df<-scale(merge_full$df, center = TRUE, scale = TRUE)
    
  })
  
  #output$autoscale_no_norm <-DT::renderDataTable({
  #  autoscale_full$df
  # })
  
  autoscale_log2_full <- reactiveValues()
  
  observeEvent(input$autoscale_log2, {
    rownames(log2_full$df)<- log2_full$df[,1]
    log2_full$df[is.na(log2_full$df)] <- 0.001
    autoscale_log2_full$df<-scale(log2_full$df[,-1], center = TRUE, scale = TRUE)
  })
  
  #    output$autoscale_log2_full <-DT::renderDataTable({
  #      autoscale_log2_full$df
  #    })
  
  
  
  #PCA for no log 2 normalised data
  pca_full <- reactive(opls(autoscale_full$df,scaleC= "none"))
  observeEvent(input$pca_full_scores,{  
    output$pca_full_nice <- renderPlot({
      pca_full_scores<- as.data.frame(getScoreMN(pca_full()))
      pca_full_nice<- ggplot(pca_full_scores,aes(x=p1,y=p2, label = row.names(pca_full_scores), colour = as.character(samp_colour())))
      pca_full_nice<- pca_full_nice + geom_point() + geom_text(size = 3)+ theme(legend.title= element_text(colour = 'black',face = 'bold')) + labs(color = input$samp_colour)
      pca_full_nice
    })
    
    
    output$pca_full_loadings <- renderGirafe({
      pca_full_load<-as.data.frame(getLoadingMN(pca_full()))
      pca_full_load$label<- row.names(pca_full_load)  
      pca_full_loadings<- ggplot(pca_full_load,aes(x=p1,y=p2,show.legend = FALSE, tooltip = label,
                                                   data_id = label)) + geom_point_interactive(size=3)
      pca_full_loadings<-girafe(ggobj = pca_full_loadings, options = list(opts_selection(type = "single", only_shiny = FALSE)))
      pca_full_loadings
    }) 
  })  
  
  #PCA for log2 normalalised data
  
  PCA_log2 <- reactive(opls(autoscale_log2_full$df,scaleC= "none"))
  observeEvent(input$PCA_log2_scores,{
    output$pca_log2_nice <- renderPlot({
      PCA_log2_scores<- as.data.frame(getScoreMN(PCA_log2()))
      pca_log2_nice<- ggplot(PCA_log2_scores,aes(x=p1,y=p2, label = row.names(PCA_log2_scores), colour = as.character(samp_colour())))
      pca_log2_nice<- pca_log2_nice + geom_point() + geom_text(size = 3)+ theme(legend.title= element_text(colour = 'black',face = 'bold')) + labs(color = input$samp_colour)
      pca_log2_nice
    })
    
    output$pca_log2_loadings <- renderGirafe({
      PCA_log2_load<-as.data.frame(getLoadingMN(PCA_log2()))
      PCA_log2_load$label<- row.names(PCA_log2_load)  
      pca_log2_loadings<- ggplot(PCA_log2_load,aes(x=p1,y=p2,show.legend = FALSE, tooltip = label,
                                                   data_id = label)) + geom_point_interactive(size=3)
      pca_log2_loadings<-girafe(ggobj = pca_log2_loadings, options = list(opts_selection(type = "single", only_shiny = FALSE)))
      pca_log2_loadings
    }) 
  })
  
  
  
  oplsda_no_norm <- reactive(opls(autoscale_full$df,as.character(samp_colour()) , scaleC= "none", predI = 1, orthoI = NA, permI = 100))
  
  observeEvent(input$oplsda_no_norm, {
    output$oplsda_no_norm_out <- renderPlot({
      oplsda_no_norm_out<- as.data.frame(merge(getScoreMN(oplsda_no_norm()),(getScoreMN(oplsda_no_norm(), orthoL = TRUE)), by.x = 0, by.y = 0))
      oplsda_no_norm_label <- as.data.frame(getScoreMN(oplsda_no_norm()))
      oplsda_no_norm_out<-ggplot(oplsda_no_norm_out,aes(x=p1,y=o1, label = row.names(oplsda_no_norm_label),colour = as.character(samp_colour())))
      oplsda_no_norm_out<-oplsda_no_norm_out + geom_text(size = 3) + geom_point()+ theme(legend.title= element_text(colour = 'black',face = 'bold')) + labs(color = input$samp_colour)
      oplsda_no_norm_out
      
    })
  })
  
  #formatting potential biomarker data for potential export
  marker_no_norm <- reactive(as.data.frame(merge(getLoadingMN(oplsda_no_norm()),(getVipVn(oplsda_no_norm())), by.x = 0, by.y = 0))) 
  order_vip_no_norm <- reactive({marker_no_norm()[order(marker_no_norm()$y, decreasing=TRUE), ]})
  top_500_marker_no_norm <- reactive(order_vip_no_norm()[1:500,])
  
  neg_load_no_norm <- reactiveValues()
  
  observeEvent(input$exp_neg_mk_no_norm, {
    neg_load_no_norm$df <- subset(top_500_marker_no_norm(), top_500_marker_no_norm()$p1 <= 0)
    names(neg_load_no_norm$df)[1] <- paste("Biomarkers for",head(samp_colour()))
    names(neg_load_no_norm$df)[2] <- 'P1 loading'
    names(neg_load_no_norm$df)[3] <- 'VIP'
  })
  
  pos_load_no_norm <- reactiveValues()
  
  observeEvent(input$exp_neg_mk_no_norm, {
    pos_load_no_norm$df <- subset(top_500_marker_no_norm(), top_500_marker_no_norm()$p1 >= 0)
    names(pos_load_no_norm$df)[1] <- paste("Biomarkers for",tail(samp_colour()))
    names(pos_load_no_norm$df)[2] <- 'P1 loading'
    names(pos_load_no_norm$df)[3] <- 'VIP'
  })
  
  
  
  #test printing dummy code
  #output$oplsda_no_norm_load <-DT::renderDataTable(pos_load_no_norm$df, rownames= FALSE)
  
  
  
  oplsda_norm <- reactive(opls(autoscale_log2_full$df,as.character(samp_colour()) , scaleC= "none", predI = 1, orthoI = NA, permI = 100))
  observeEvent(input$oplsda_norm, {
    output$oplsda_norm_out <- renderPlot({
      oplsda_norm_out<- as.data.frame(merge(getScoreMN(oplsda_norm()),(getScoreMN(oplsda_norm(), orthoL = TRUE)), by.x = 0, by.y = 0))
      oplsda_norm_label <- as.data.frame(getScoreMN(oplsda_norm()))
      oplsda_norm_out<-ggplot(oplsda_norm_out,aes(x=p1,y=o1,label = row.names(oplsda_norm_label), colour = as.character(samp_colour())))
      oplsda_norm_out<-oplsda_norm_out +geom_point()+ geom_text(size = 3) + theme(legend.title= element_text(colour = 'black',face = 'bold')) + labs(color = input$samp_colour)
      oplsda_norm_out
      
    })
  })
  
  
  
  
}

shinyApp(ui = ui, server = server)


