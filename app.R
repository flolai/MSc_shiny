#This is a R shiny app for data integration of transcriptomics and metabolomics data set
#Author: Florence Lai 2020

if(!require('pacman'))install.packages('pacman')
pacman::p_load(shiny,dplyr,reshape,plotly,tidyverse,ggiraph, gridExtra, DT)

library(ropls)
citation("ropls")

# ui object
ui <- fluidPage(
  titlePanel("TMMS-GUI: Transcriptmics, metabolomics multivatiate statistical analysis GUI"),
  
  tabsetPanel(
    
    #First tab to read multiple user CSV input of data, which include
    #metabolomics data, transcriptomics data, compound name and gene name all in csv format    
    
    tabPanel("Data Input",
             tags$style(type="text/css",
                        #".shiny-output-error { visibility: hidden; }",
                        ".shiny-output-error:after { visibility: hidden; }"
                        
             ),  
             fluidRow(column(width = 2,
                             tags$hr(),
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
             
             column(width = 6,
                    (includeMarkdown("My First Shiny App.md"))
             )
             )
             
    ),
    
    #Second tab to read sample information (e.g treated/control) and provide user with choice of normalisation or not
    #PCA scores plot and loadings plot will be generated. Loading plot have options to select clustered genes and metabolites
    #selected genes and metabolites will be displayed in GUI
    #OPLS-DA
    
    tabPanel("Plots",
             tags$hr(),
             sidebarPanel(
               width = 3,
             fileInput(inputId = "samp_detail",
                       label = "Sample details in csv format",
                       multiple = TRUE,
                       accept = c(".csv")),
             #adding selection for sample discription classification for colouring
             
             selectizeInput('samp_colour', 'Variable selection for sample details',
                            choices = NULL,
                            options = list(placeholder = 'Select from below'),
                            multiple = FALSE),
               tags$hr(),
               actionButton("pca_full_scores", "PCA for pre-normalised data"),
               tags$hr(),
               actionButton("oplsda_no_norm", "OPLS-DA model for pre-normalised data"),
               
               tags$hr(),
               actionButton("PCA_log2_scores", "PCA for log2 normalised data"),
               tags$hr(),
               actionButton("oplsda_norm", "OPLS-DA model for log2 normalised data"),
               tags$hr(),
             ),
             mainPanel(
               conditionalPanel(
                 condition = ("input.pca_full_scores != 0"),
                 fluidRow(
                   column(width = 5, 
                          h4("PCA scores plot for non normalised data"),
                          plotOutput("pca_full_nice")), 
                   column(width = 5, 
                          h4("PCA loadings plot for non normalised data"),
                          girafeOutput("pca_full_loadings")),     
                   
                   column(width =2 ,
                          h4("Selected genes and metabolites"),
                          DT::dataTableOutput("datatab")),
                 ),
               h6("PCA model P1 and P2 varance explained percentage"),
               tableOutput("no_norm_PCA_R2X"),
               ),
               
               conditionalPanel(
                 condition = ("input.PCA_log2_scores != 0"),
                 fluidRow(
                   column(width = 5, 
                          h4("PCA scores plot for normalised data"),
                          plotOutput("pca_log2_nice")), 
                   column(width = 5, 
                          h4("PCA loadings plot for normalised data"),
                          girafeOutput("pca_log2_loadings")),     
                   
                   column(width =2 ,
                          h4("Selected genes and metabolites"),
                          DT::dataTableOutput("datatab_log2")),
                 ),
                 h6("PCA model P1 and P2 varance explained percentage"),
                 tableOutput("norm_PCA_R2X"),
               ),
               
               
               conditionalPanel(
                 condition = ("input.oplsda_no_norm != 0"),
                 fluidRow(
                   column(width = 5, 
                          h4("OPLS-DA scores plot for non-normalised data"),
                          plotOutput("oplsda_no_norm_out")), 
                   column(width = 5, 
                          h4("Permutation model validation Plot"),
                          plotOutput("perm_no_norm_plot"))
                 )
               ),
               
               tableOutput("no_norm_sum"),
               
               conditionalPanel(
                 condition = ("input.oplsda_norm != 0"),
                 fluidRow(
                   column(width = 5, 
                          h4("OPLS-DA scores plot for normalised data"),
                          plotOutput("oplsda_norm_out")), 
                   column(width = 5, 
                          h4("Permutation model validation Plot"),
                          plotOutput("perm_norm_plot"))
                 )
               ),
               tableOutput("norm_sum"),
             ),
    ),
    tabPanel("Biomarkers",
             mainPanel( 
               tags$hr(),
               conditionalPanel(
                 condition = ("input.oplsda_no_norm != 0"),
                 downloadButton("download_gp1_no_norm", "Download group 1 (left) genes and metabolites"),
                 downloadButton("download_gp2_no_norm", "Download group 2 (right) genes and metabolites"),
               ),
               conditionalPanel(
                 condition = ("input.oplsda_norm != 0"),
                 downloadButton("download_gp1", "Download group 1 (left) genes and metabolites"),
                 downloadButton("download_gp2", "Download group 2 (right) genes and metabolites"),
               ),
               tags$hr(),             
               DT::dataTableOutput("biomarkers_gp1_no_norm"),
               DT::dataTableOutput("biomarkers_gp2_no_norm"),
               DT::dataTableOutput("biomarkers_gp1"),
               DT::dataTableOutput("biomarkers_gp2"),
             ),
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
  
  
  
  comp_name <- reactive({
    if (is.null(input$comp_name)){
      return(NULL)
    }
    read.csv(input$comp_name$datapath, header = FALSE)
    
  })
  
  gene_name <- reactive({
    if (is.null(input$gene_name)){
      return(NULL)
    }
    read.csv(input$gene_name$datapath, header = FALSE)
    
  })
  
  
  all_names <- reactive({
    temp_comp <- t(meta())
    temp_comp <- as.data.frame(row.names(temp_comp))
    temp_gene <- t(ngs())
    temp_gene <- as.data.frame(row.names(temp_gene))
    colnames(temp_gene)[1] <- 'ID'
    colnames(temp_comp)[1] <- 'ID'
    comp_gene <- rbind(temp_comp, temp_gene) # data part of gene and compound name, only taking the row name as no data is needed at this point
    gene_name <- cbind(as.data.frame(str_replace(gene_name()$V1, "-", ".")), gene_name()$V2) #in shiny, this is not needed
    colnames(gene_name)[1] <- 'V1' #ropls automatically change '-' in any name to '.' quick fix for now
    colnames(gene_name)[2] <- 'V2' #ropls automatically change '-' in any name to '.'
    all_comp_gene <- rbind(comp_name(),gene_name) #merging of gene name and compound names given
    all_names <- merge(comp_gene, all_comp_gene, by.x = 1, by.y = 1)
    colnames(all_names)[2] <- 'Name'
    all_names<- as.data.frame(all_names)
  })
  
  
  
  
  #reading sample details for classification    
  samp_detail <- reactive({
    if (is.null(input$samp_detail)){
      return(NULL)
    }
    read.csv(input$samp_detail$datapath, header = TRUE , row.names = 1) 
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
  
  observeEvent(input$pca_full_scores, {
    merge_full$df <- transform(merge(meta(), ngs(), by.x = 0, by.y = 0), row.names = Row.names, Row.names = NULL)# merging metabolomics
    #and ngs data by row names = sample name
  })
  
  log2_full <- reactiveValues()
  
  observeEvent(input$PCA_log2_scores, {
    log2_full$df <- transform(merge(log2(meta()), log2(ngs()), by.x = 0, by.y = 0, row.names = Row.names, Row.names = NULL)) #log2 normalisaion
  })
  
  
  autoscale_full <- reactiveValues()
  
  observeEvent(input$pca_full_scores, {
    merge_full$df[is.na(merge_full$df)] <- 0.001
    autoscale_full$df<-scale(merge_full$df, center = TRUE, scale = TRUE)
    autoscale_full$df[is.na(autoscale_full$df)] <- 0.001
  })
  
  
  
  autoscale_log2_full <- reactiveValues()
  
  observeEvent(input$PCA_log2_scores, {
    rownames(log2_full$df)<- log2_full$df[,1]
    log2_full$df <- log2_full$df[,-1]
    autoscale_log2_full$df<-scale(log2_full$df, center = TRUE, scale = TRUE)
    autoscale_log2_full$df[is.na(autoscale_log2_full$df)] <- 0.001
  })
  
  
  #PCA for no log 2 normalised data
  pca_full <- reactive(opls(autoscale_full$df,scaleC= "none"))
  observeEvent(input$pca_full_scores,{  
    output$pca_full_nice <- renderPlot({
      pca_full_scores<- as.data.frame(getScoreMN(pca_full()))
      pca_full_nice<- ggplot(pca_full_scores,aes(x=p1,y=p2, label = row.names(pca_full_scores), colour = as.character(samp_colour())))
      pca_full_nice<- pca_full_nice + geom_point(position = position_nudge(y = -1.8)) + geom_text(size = 3)+ theme(legend.title= element_text(colour = 'black',face = 'bold')) + labs(color = input$samp_colour)
      pca_full_nice
    })
    
    output$pca_full_loadings <- renderGirafe({
      pca_full_load<-as.data.frame(getLoadingMN(pca_full()))
      pca_full_load$label <- row.names(pca_full_load)  
      pca_full_loadings <- ggplot(pca_full_load) + geom_point_interactive(aes(x=p1,y=p2), size=3, tooltip=pca_full_load$label,
                                                                          data_id = pca_full_load$label)
      pca_full_loadings <- girafe(ggobj = pca_full_loadings, options = list(opts_zoom(max = 5)))
      
      pca_full_loadings <- girafe_options(pca_full_loadings,
                                          options = list(
                                            opts_selection(type = "multiple", css = "fill:red;stroke:gray;r:5pt;")
                                          )
      )
      pca_full_loadings
    }) 
    
  })  
  
  
  selected_marks <- reactive({
    if (is.null(input$pca_full_loadings_selected)){
      return(print('Please select genes and metabolite'))
  }
    input$pca_full_loadings_selected
  
  })
  
  
  
  output$datatab <- DT::renderDataTable({options = list(scrollX = TRUE)
  marks <- as.data.frame(row.names(getLoadingMN(pca_full())))
  out <- as.data.frame(selected_marks())
  names(out)[1] <- paste("Selected")
  out
  })
  
  
  observeEvent(input$pca_full_scores, {
    output$no_norm_PCA_R2X <- renderTable( rownames = FALSE, colnames = F,{
      no_norm_PCA_R2X <- t(as.data.frame(pca_full()@modelDF[["R2X"]]*100)[c(1:2),])
      no_norm_PCA_R2X <- cbind('Percentage of variance explained', no_norm_PCA_R2X)
      no_norm_PCA_R2X <- rbind(cbind('', 'P1','P2'),no_norm_PCA_R2X)
      no_norm_PCA_R2X
    })
  })
  
  #PCA for log2 normalalised data
  
  PCA_log2 <- reactive(opls(autoscale_log2_full$df,scaleC= "none"))
  observeEvent(input$PCA_log2_scores,{
    output$pca_log2_nice <- renderPlot({
      PCA_log2_scores<- as.data.frame(getScoreMN(PCA_log2()))
      pca_log2_nice<- ggplot(PCA_log2_scores,aes(x=p1,y=p2, label = row.names(PCA_log2_scores), colour = as.character(samp_colour())))
      pca_log2_nice<- pca_log2_nice + geom_point(position = position_nudge(y = -1.8)) + geom_text(size = 3)+ theme(legend.title= element_text(colour = 'black',face = 'bold')) + labs(color = input$samp_colour)
      pca_log2_nice
    })
    
    output$pca_log2_loadings <- renderGirafe({
      PCA_log2_load<-as.data.frame(getLoadingMN(PCA_log2()))
      PCA_log2_load$label<- row.names(PCA_log2_load)  
      pca_log2_loadings<- ggplot(PCA_log2_load,aes(x=p1,y=p2,show.legend = FALSE, tooltip = label,
                                                   data_id = label)) + geom_point_interactive(size=3)
      
      pca_log2_loadings <- girafe(ggobj = pca_log2_loadings, options = list(opts_zoom(max = 5)))
      
      pca_log2_loadings <- girafe_options(pca_log2_loadings,
                                          options = list(
                                            opts_selection(type = "multiple", css = "fill:red;stroke:gray;r:5pt;")
                                          )
      )
      pca_log2_loadings
    }) 
  })
  
  
  selected_marks_log2 <- reactive({
    if (is.null(input$pca_log2_loadings_selected)){
      return(print('Please select genes and metabolite'))
    }
    input$pca_log2_loadings_selected
    
  })
  
  
  
  output$datatab_log2 <- DT::renderDataTable({options = list(scrollX = TRUE)
  marks <- as.data.frame(row.names(getLoadingMN(PCA_log2())))
  out <- as.data.frame(selected_marks_log2())
  names(out)[1] <- paste("Selected")
  out
  })
  
  observeEvent(input$PCA_log2_scores, {
    output$norm_PCA_R2X <- renderTable(rownames = FALSE, colnames = F,{
      norm_PCA_R2X <- t(as.data.frame(PCA_log2()@modelDF[["R2X"]]*100)[c(1:2),])
      norm_PCA_R2X <- cbind('Percentage of variance explained', norm_PCA_R2X)
      norm_PCA_R2X <- rbind(cbind('', 'P1','P2'), norm_PCA_R2X)
      norm_PCA_R2X
    })
  })
  
  
  
  
  oplsda_no_norm <- reactive(opls(autoscale_full$df,as.character(samp_colour()) , scaleC= "none", predI = 1, orthoI = NA ,permI = 100))
  
  #adding OPLS-DA scores plot and permutation test
  observeEvent(input$oplsda_no_norm, {
    output$oplsda_no_norm_out <- renderPlot({
      oplsda_no_norm_out<- as.data.frame(merge(getScoreMN(oplsda_no_norm()),(getScoreMN(oplsda_no_norm(), orthoL = TRUE)), by.x = 0, by.y = 0))
      oplsda_no_norm_label <- as.data.frame(getScoreMN(oplsda_no_norm()))
      oplsda_no_norm_out<-ggplot(oplsda_no_norm_out,aes(x=p1,y=o1, label = row.names(oplsda_no_norm_label),colour = as.character(samp_colour())))
      oplsda_no_norm_out<-oplsda_no_norm_out + geom_text(size = 3) + geom_point(position = position_nudge(y = -1.8))+ theme(legend.title= element_text(colour = 'black',face = 'bold')) + labs(color = input$samp_colour)
      oplsda_no_norm_out
      
    })
    
    perm_no_norm <- reactiveValues()
    
    observeEvent(input$oplsda_no_norm, {
      perm_no_norm$df <- as.data.frame(oplsda_no_norm()@suppLs[["permMN"]])
      perm_no_norm$df <- subset(perm_no_norm$df, select= c(`R2Y(cum)`, `Q2(cum)`, sim))
      perm_no_norm$df <- melt(perm_no_norm$df, id = 'sim')
      
      
    })
    
    output$perm_no_norm_plot <- renderPlot({
      perm_no_norm_plot <- ggplot(perm_no_norm$df, aes(sim, value, col = variable)) + geom_point()
      no_norm_sum <- tableGrob(as.data.frame(getSummaryDF(oplsda_no_norm())))
      perm_no_norm_plot
    })
  })
  observeEvent(input$oplsda_no_norm, {
    output$no_norm_sum <- renderTable({
      no_norm_sum <- as.data.frame(getSummaryDF(oplsda_no_norm()))
    })
  })
  
  
  #formatting potential biomarker data for potential export
  marker_no_norm <- reactive(as.data.frame(merge(getLoadingMN(oplsda_no_norm()),(getVipVn(oplsda_no_norm())), by.x = 0, by.y = 0))) 
  order_vip_no_norm <- reactive({marker_no_norm()[order(marker_no_norm()$y, decreasing=TRUE), ]})
  top_marker_no_norm <- reactive(order_vip_no_norm()[1:1000,])
  
  biomarkers_gp1_no_norm <- reactiveValues()
  biomarkers_gp2_no_norm <- reactiveValues()
  neg_load_no_norm <- reactiveValues()
  
  observeEvent(input$oplsda_no_norm,{
    neg_load_no_norm$df <- subset(top_marker_no_norm(), top_marker_no_norm()$p1 <= 0)
    names(neg_load_no_norm$df)[1] <- paste("Biomarkers for left hand group from scores plot")
    names(neg_load_no_norm$df)[2] <- 'P1 loading'
    names(neg_load_no_norm$df)[3] <- 'VIP'
    neg_load_no_norm$df
  })
  
  
  #to view biomarkers
  observeEvent(input$oplsda_no_norm, {
    biomarkers_gp1_no_norm$df <- as.data.frame(merge(neg_load_no_norm$df,all_names(), by.x = 1, by.y = 1))
    biomarkers_gp1_no_norm$df <- biomarkers_gp1_no_norm$df[, -c(2)]
  })
  
  pos_load_no_norm <- reactiveValues()
  
  observeEvent(input$oplsda_no_norm, {
    pos_load_no_norm$df <- subset(top_marker_no_norm(), top_marker_no_norm()$p1 >= 0)
    names(pos_load_no_norm$df)[1] <- paste("Biomarkers for right hand group from scores plot")
    names(pos_load_no_norm$df)[2] <- 'P1 loading'
    names(pos_load_no_norm$df)[3] <- 'VIP'
  })
  
  observeEvent(input$oplsda_no_norm, {
    biomarkers_gp2_no_norm$df <- as.data.frame(merge(pos_load_no_norm$df,all_names(), by.x = 1, by.y = 1))
    biomarkers_gp2_no_norm$df <- biomarkers_gp2_no_norm$df[, -c(2)]
    
  })
  
  
  #normalised section
  
  oplsda_norm <- reactive(opls(autoscale_log2_full$df,as.character(samp_colour()) , scaleC= "none", predI = 1, orthoI = NA, permI = 100))
  
  observeEvent(input$oplsda_norm, {
    output$oplsda_norm_out <- renderPlot({
      oplsda_norm_out<- as.data.frame(merge(getScoreMN(oplsda_norm()),(getScoreMN(oplsda_norm(), orthoL = TRUE)), by.x = 0, by.y = 0))
      oplsda_norm_label <- as.data.frame(getScoreMN(oplsda_norm()))
      oplsda_norm_out<-ggplot(oplsda_norm_out,aes(x=p1,y=o1,label = row.names(oplsda_norm_label), colour = as.character(samp_colour())))
      oplsda_norm_out<-oplsda_norm_out +geom_point(position = position_nudge(y = -1.8))+ geom_text(size = 3) + theme(legend.title= element_text(colour = 'black',face = 'bold')) + labs(color = input$samp_colour)
      oplsda_norm_out
      
    })
  })
  
  perm_norm <- reactiveValues()
  
  observeEvent(input$oplsda_norm, {
    perm_norm$df <- as.data.frame(oplsda_norm()@suppLs[["permMN"]])
    perm_norm$df <- subset(perm_norm$df, select= c(`R2Y(cum)`, `Q2(cum)`, sim))
    perm_norm$df <- melt(perm_norm$df, id = 'sim')
    
    
  })
  
  output$perm_norm_plot <- renderPlot({
    perm_norm_plot <- ggplot(perm_norm$df, aes(sim, value, col = variable)) + geom_point()
    norm_sum <- tableGrob(as.data.frame(getSummaryDF(oplsda_norm())))
    perm_norm_plot
  })
  
  observeEvent(input$oplsda_norm, {
    output$norm_sum <- renderTable({
      norm_sum <- as.data.frame(getSummaryDF(oplsda_norm()))
    })
  })
  
  
  
  #  formatting potential biomarker data for potential export - normalised data
  marker_norm <- reactive(as.data.frame(merge(getLoadingMN(oplsda_norm()),(getVipVn(oplsda_norm())), by.x = 0, by.y = 0))) 
  order_vip_norm <- reactive({marker_norm()[order(marker_norm()$y, decreasing=TRUE), ]})
  top_marker_norm <- reactive(order_vip_norm()[1:1000,])
  
  biomarkers_gp1 <- reactiveValues()
  biomarkers_gp2 <- reactiveValues()
  neg_load_norm <- reactiveValues()
  
  observeEvent(input$oplsda_norm,{
    neg_load_norm$df <- subset(top_marker_norm(), top_marker_norm()$p1 <= 0)
    names(neg_load_norm$df)[1] <- paste("Biomarkers for left hand group from scores plot")
    names(neg_load_norm$df)[2] <- 'P1 loading'
    names(neg_load_norm$df)[3] <- 'VIP'
    neg_load_norm$df
  })
  
  
  #to view biomarkers
  observeEvent(input$oplsda_norm, {
    biomarkers_gp1$df <- as.data.frame(merge(neg_load_norm$df,all_names(), by.x = 1, by.y = 1))
    biomarkers_gp1$df <- biomarkers_gp1$df[, -c(2)]
  })
  
  pos_load_norm <- reactiveValues()
  
  observeEvent(input$oplsda_norm, {
    pos_load_norm$df <- subset(top_marker_norm(), top_marker_norm()$p1 >= 0)
    names(pos_load_norm$df)[1] <- paste("Biomarkers for right hand group from scores plot")
    names(pos_load_norm$df)[2] <- 'P1 loading'
    names(pos_load_norm$df)[3] <- 'VIP'
    
  })
  
  observeEvent(input$oplsda_norm, {
    biomarkers_gp2$df <- as.data.frame(merge(pos_load_norm$df,all_names(), by.x = 1, by.y = 1))
    biomarkers_gp2$df <- biomarkers_gp2$df[, -c(2)]
  })
  
  
  
  #Printing biomarkers on the page
  output$biomarkers_gp1_no_norm <-DT::renderDataTable(biomarkers_gp1_no_norm$df, rownames= FALSE)
  output$biomarkers_gp2_no_norm <-DT::renderDataTable(biomarkers_gp2_no_norm$df, rownames= FALSE)
  
  output$biomarkers_gp1 <-DT::renderDataTable(biomarkers_gp1$df, rownames= FALSE)
  output$biomarkers_gp2 <-DT::renderDataTable(biomarkers_gp2$df, rownames= FALSE)
  
  output$download_gp1_no_norm <- downloadHandler(
    filename = function() {
      paste(input$biomarkers_gp1_no_norm, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(biomarkers_gp1_no_norm$df, file, row.names = FALSE)
    }
  )
  
  output$download_gp2_no_norm <- downloadHandler(
    filename = function() {
      paste(input$biomarkers_gp2_no_norm, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(biomarkers_gp2_no_norm$df, file, row.names = FALSE)
    }
  )
  
  output$download_gp1 <- downloadHandler(
    filename = function() {
      paste(input$biomarkers_gp1, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(biomarkers_gp1$df, file, row.names = FALSE)
    }
  )
  
  output$download_gp2 <- downloadHandler(
    filename = function() {
      paste(input$biomarkers_gp2, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(biomarkers_gp2$df, file, row.names = FALSE)
    }
  )
  
}

shinyApp(ui = ui, server = server)


