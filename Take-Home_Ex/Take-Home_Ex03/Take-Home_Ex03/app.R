# Load required packages
pacman::p_load(shiny, spdep, tmap, sf, ClustGeo, 
               ggpubr, cluster, factoextra, NbClust,
               heatmaply, corrplot, psych, tidyverse, GGally)

# Load data for use in the app
shan_sf <- st_read(dsn = "data/geospatial", 
                   layer = "myanmar_township_boundaries") %>%
  filter(ST %in% c("Shan (East)", "Shan (North)", "Shan (South)")) %>%
  select(c(2:7))

ict <- read_csv("data/aspatial/Shan-ICT.csv")

ict_derived <- ict %>%
  mutate(`RADIO_PR` = `Radio`/`Total households`*1000,
         `TV_PR` = `Television`/`Total households`*1000,
         `LLPHONE_PR` = `Land line phone`/`Total households`*1000,
         `MPHONE_PR` = `Mobile phone`/`Total households`*1000,
         `COMPUTER_PR` = `Computer`/`Total households`*1000,
         `INTERNET_PR` = `Internet at home`/`Total households`*1000) %>%
  rename(`DT_PCODE` = `District Pcode`,
         `DT` = `District Name`,
         `TS_PCODE` = `Township Pcode`, 
         `TS` = `Township Name`,
         `TT_HOUSEHOLDS` = `Total households`,
         `RADIO` = `Radio`, 
         `TV` = `Television`, 
         `LLPHONE` = `Land line phone`, 
         `MPHONE` = `Mobile phone`,
         `COMPUTER` = `Computer`, 
         `INTERNET` = `Internet at home`)

shan_sf <- left_join(shan_sf, ict_derived, by = c("TS_PCODE" = "TS_PCODE"))

# Hierarchical Cluster Analysis
cluster_vars <- shan_sf %>%
  st_set_geometry(NULL) %>%
  select("TS.x", "RADIO_PR", "TV_PR", "LLPHONE_PR", "MPHONE_PR", "COMPUTER_PR")

row.names(cluster_vars) <- cluster_vars$"TS.x"

shan_ict <- select(cluster_vars, c(2:6))

shan_ict.std <- normalize(shan_ict)
shan_ict.z <- scale(shan_ict)

# Define UI
ui <- fluidPage(
  titlePanel("Shan Region ICT Data Exploration"),
  
  sidebarLayout(
    sidebarPanel(
      conditionalPanel(
        condition = "input.tabSelected == 'Histogram' || input.tabSelected == 'Boxplot' || input.tabSelected == 'Map'",
        selectInput("variable", "Select Variable to Display", 
                    choices = c("RADIO_PR", "TV_PR", "LLPHONE_PR", "MPHONE_PR", "COMPUTER_PR", "INTERNET_PR"))
      ),
      conditionalPanel(
        condition = "input.tabSelected == 'Histogram'",
        sliderInput("bins", "Number of Bins for Histogram", min = 5, max = 50, value = 20)
      ),
      conditionalPanel(
        condition = "input.tabSelected == 'Standardization'",
        selectInput("standard_var", "Select Variable for Standardization Comparison", 
                    choices = c("RADIO_PR", "TV_PR", "LLPHONE_PR", "MPHONE_PR", "COMPUTER_PR")),
        radioButtons("plot_type", "Choose Plot Type", 
                     choices = c("Histogram" = "histogram", "Density" = "density"), 
                     selected = "density")
      ),
      conditionalPanel(
        condition = "input.tabSelected == 'Cluster'",
        selectInput("distance_metric", "Select Distance Metric", 
                    choices = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"),
                    selected = "euclidean")
      ),
      actionButton("refresh", "Refresh Data")
    ),
    
    mainPanel(
      tabsetPanel(
        id = "tabSelected",
        tabPanel("Histogram", value = "Histogram", plotOutput("histPlot")),
        tabPanel("Boxplot", value = "Boxplot", plotOutput("boxPlot")),
        tabPanel("Map", value = "Map", tmapOutput("mapPlot")),
        tabPanel("Correlation Analysis", value = "Correlation", plotOutput("corrPlot")),
        tabPanel("Hierarchy Cluster Analysis", value = "Cluster", 
                 plotOutput("clusterPlot"), 
                 plotOutput("heatmapPlot")),
        tabPanel("Standardization Comparison", value = "Standardization", plotOutput("standardizationPlot"))
      )
    )
  )
)

# Define the server logic
server <- function(input, output, session) {
  
  # Reactive expression to update plots based on selected variable
  reactive_data <- reactive({
    variable <- input$variable
    ict_derived %>%
      select(all_of(variable)) %>%
      pull()
  })
  
  # Render histogram based on selected variable
  output$histPlot <- renderPlot({
    data <- reactive_data()
    ggplot(ict_derived, aes(x = data)) +
      geom_histogram(bins = input$bins, color = "black", fill = "lightblue") +
      xlab(input$variable) +
      ggtitle(paste("Histogram of", input$variable))
  })
  
  # Render boxplot based on selected variable
  output$boxPlot <- renderPlot({
    data <- reactive_data()
    ggplot(ict_derived, aes(x = data)) +
      geom_boxplot(color = "black", fill = "lightblue") +
      xlab(input$variable) +
      ggtitle(paste("Boxplot of", input$variable))
  })
  
  # Render map based on selected variable
  output$mapPlot <- renderTmap({
    tm_shape(shan_sf) +
      tm_fill(col = input$variable, style = "jenks") +
      tm_borders(alpha = 0.5) +
      tm_layout(legend.position = c("right", "bottom"))
  })
  
  # Render correlation plot
  output$corrPlot <- renderPlot({
    cluster_vars.cor <- cor(ict_derived[,12:17])
    corrplot.mixed(cluster_vars.cor,
                   lower = "ellipse", 
                   upper = "number",
                   tl.pos = "lt",
                   diag = "l",
                   tl.col = "black")
  })
  
  # Render hierarchical cluster analysis dendrogram
  output$clusterPlot <- renderPlot({
    # Get the selected distance metric
    dist_method <- input$distance_metric
    
    # Calculate proximity matrix using the selected distance metric
    proxmat <- dist(shan_ict, method = dist_method)
    
    # Perform hierarchical clustering using Ward's method
    hclust_ward <- hclust(proxmat, method = 'ward.D')
    
    # Plot the dendrogram
    plot(hclust_ward, cex = 0.6, main = paste("Dendrogram: Ward's Method (", dist_method, ")", sep = ""))
    
    # Highlight clusters with colored borders
    rect.hclust(hclust_ward, k = 6, border = 2:5)
  })
  
  # Render heatmap for clusters
  output$heatmapPlot <- renderPlot({
    shan_ict_mat <- data.matrix(shan_ict)
    heatmaply(normalize(shan_ict_mat),
              Colv = NA,
              dist_method = input$distance_metric,
              hclust_method = "ward.D",
              seriate = "OLO",
              colors = Blues,
              k_row = 6,
              margins = c(NA, 200, 60, NA),
              fontsize_row = 4,
              fontsize_col = 5,
              main = "Geographic Segmentation of Shan State by ICT indicators",
              xlab = "ICT Indicators",
              ylab = "Townships of Shan State")
  })
  
  # Render Standardization Comparison plot
  output$standardizationPlot <- renderPlot({
    # Reactive selection for variable
    var_selected <- input$standard_var
    
    # Choose between histogram or density plot based on user selection
    if (input$plot_type == "histogram") {
      r <- ggplot(ict_derived, aes_string(x = var_selected)) +
        geom_histogram(bins = 20, color = "black", fill = "lightblue") +
        ggtitle("Raw values without standardisation")
      
      s <- ggplot(as.data.frame(shan_ict.std), aes_string(x = var_selected)) +
        geom_histogram(bins = 20, color = "black", fill = "lightblue") +
        ggtitle("Min-Max Standardisation")
      
      z <- ggplot(as.data.frame(shan_ict.z), aes_string(x = var_selected)) +
        geom_histogram(bins = 20, color = "black", fill = "lightblue") +
        ggtitle("Z-score Standardisation")
      
    } else {
      r <- ggplot(ict_derived, aes_string(x = var_selected)) +
        geom_density(color = "black", fill = "lightblue") +
        ggtitle("Raw values without standardisation")
      
      s <- ggplot(as.data.frame(shan_ict.std), aes_string(x = var_selected)) +
        geom_density(color = "black", fill = "lightblue") +
        ggtitle("Min-Max Standardisation")
      
      z <- ggplot(as.data.frame(shan_ict.z), aes_string(x = var_selected)) +
        geom_density(color = "black", fill = "lightblue") +
        ggtitle("Z-score Standardisation")
    }
    
    # Arrange plots side by side
    ggarrange(r, s, z, ncol = 3, nrow = 1)
  })
  
  observeEvent(input$refresh, {
    # Add any data refresh logic if needed
  })
}

# Run the Shiny app
shinyApp(ui = ui, server = server)
