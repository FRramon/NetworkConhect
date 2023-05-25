
#library(ggraph)

# Données d'exemple
#source("/home/imabrain/Documents/Brain_networks/vAllWeights/scripts/1_load_data.R",local=TRUE)
#source("/home/imabrain/Documents/Brain_networks/vAllWeights/scripts/2_graph_construction.R",local=TRUE)
#source("/home/imabrain/Documents/Brain_networks/vAllWeights/scripts/3_compute_metrics.R",local=TRUE)


###############################################################################################
# VISUALIZATION TOOL, THRESHOLDING, DISPLAY
###############################################################################################




#' Generate an interface to plot graph for a patient & visit & weighting scheme.
#'
#' @param metric a chr for the weighting scheme : "FA", "FBC", "GFA"...
#' @returns a shiny interface
#' @examples
#' interfaceConhect("FA")
#' @export
interfaceConhect <- function(metric){
  data_path <- getDataDir(metric)
  print(data_path)
  whole_data <-read_and_normalize_data(data_path,metric)
  max_val <- max(whole_data$weight)

    # Interface graphique
  ui <- fluidPage(

      # App title ----
    titlePanel(paste("Réseau",metric)),

      # Sidebar layout with input and output definitions ----
    sidebarLayout(

        # Sidebar panel for inputs ----
      sidebarPanel(

          # Input: Select the random distribution type ----
        selectInput("subject_id", "Patient : ", choices = unique(whole_data$subject_id)),
        selectInput("visit_id", "Numéro visite : ", choices = unique(whole_data$visit_id)),

        br(),
          # Input: Slider for the number of observations to generate ----
        sliderInput("n",
                    "Seuil :",
                    value = max_val/2,
                    min = 0,
                    max = max_val),
        br(),
        actionButton("Go", "Apply changes")),

        # Main panel for displaying outputs ----
        mainPanel(

          # Output: Tabset w/ plot, summary, and table ----
        tabsetPanel(type = "tabs",
                    tabPanel("Plot", plotOutput("plot")),
                    tabPanel("Strength distribution", plotOutput("plot4")),
                    tabPanel("Edge weight distribution", plotOutput("plot5")),
                    tabPanel("Plot Hubs", plotOutput("plothubs")),
                    tabPanel("Hubs labels", dataTableOutput("hubnames")),

                    tabPanel("Adjacency matrix", plotOutput("plot_adj")),
                    tabPanel("Summary", dataTableOutput("summary")),
                    tabPanel("LUT", dataTableOutput("lut_table")),
                    tabPanel("Table", dataTableOutput("summary_table"))
        )
      )
    )
  )

  server <- function(input, output) {
    #source("/home/imabrain/Documents/Brain_networks/vAllWeights/scripts/1_load_data.R",local=TRUE)
    #source("/home/imabrain/Documents/Brain_networks/vAllWeights/scripts/2_graph_construction.R",local=TRUE)
    #source("/home/imabrain/Documents/Brain_networks/vAllWeights/scripts/3_compute_metrics.R",local=TRUE)


    data <- eventReactive(input$Go,{
      edgelist <- makeEdgelist(whole_data,input$subject_id,input$visit_id,metric,input$n)
    })

    output$summary_table <- renderDataTable({
      data()
    })

    output$lut_table <- renderDataTable({
      rois <- c(data()$from,data()$to)
      rois_unique <- data.frame(No = unique(rois))
      lut<- getLUT()
      merge(rois_unique,lut,by="No")
    })

    output$plot <- renderPlot({
      g<- graph_from_data_frame(data(),directed=F)
      plot(g)
    })

    output$plot_adj <- renderPlot({
      g<- graph_from_data_frame(data(),directed=F)
      plot_adjacency(g,metric)
    })

    output$plot4 <- renderPlot({
      g<- graph_from_data_frame(data(),directed=F)
      degg<-strength(g)
      qplot(degg,geom='histogram',xlab="Strength")
    })

    output$plot5 <- renderPlot({
      g<- graph_from_data_frame(data(),directed=F)
      edge_density<-E(g)$weight
      qplot(edge_density,geom='histogram',xlab="Edge weights")
    })
    output$plothubs <- renderPlot({
      g<- graph_from_data_frame(data(),directed=F)
      par(mfrow = c(1,3))
      h1 <- hub_detection(g,"strength")
      h2 <- hub_detection(g,"betweenness")
      h3 <- hub_detection(g,"closeness")

      #print(unique(c(h1$labels[,2],h2$labels[,2],h3$labels[,2])))
    })
    output$hubnames <- renderDataTable({
      g <-  graph_from_data_frame(data(),directed=F)
      h1 <- hub_detection(g,"strength")
      h2 <- hub_detection(g,"betweenness")
      h3 <- hub_detection(g,"closeness")

      df <- data.frame(unique(c(h1$labels[,2],h2$labels[,2],h3$labels[,2])))
      colnames(df) <- "Hub names"
      df
    })

    output$summary <- renderDataTable({
      g <- graph_from_data_frame(data(),directed=F)
      gR <- sample_gnm(n = vcount(g),m=ecount(g))
      clustering_coeffR <- transitivity(gR)
      clustering_coeff<-transitivity(g)
      normalized_clustering_coeff <- clustering_coeff/clustering_coeffR
      normalized_shortest_path<-mean_distance(g)/mean_distance(gR)
      shortest_path<-mean_distance(g)
      S <-data.frame(clustering_coeff,shortest_path,normalized_shortest_path,normalized_clustering_coeff)
      S
    })
  }
  shinyApp(ui = ui, server = server)

}
# Lancement de l'application



