library(shiny)
library(cluster)
library(Rtsne)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(factoextra)

#load datasets
filelist <- list.files(path = "./data2/", pattern = ".*.csv")
filenames <- gsub(".csv", "", filelist) 
filenames <- gsub("wasserstein_", "", filenames) # filtration

# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("Ivesheadiomorph Clustering"),
  sidebarLayout(
    
    sidebarPanel(
      "t-SNE visualisation of 2nd wasserstein distances of Persistence Diagrams",
      selectInput("dataset", label = "Filtration", choices = filenames),
      selectInput("method", label = "Method", choices = c("pam", "kmeans","hclust")),
      sliderInput("nclust",
                  "Number of clusters:",
                  min = 0,
                  max = 5,
                  value = 2),
      sliderInput("per",
                "t-SNE Perplexity:",
                min = 2,
                max = 8,
                value = 5)
    ),
    
    
    mainPanel(
      plotOutput(outputId = "clustPlot")
    )
  )

)

# Define server logic 
server <- function(input, output) {
  
  output$clustPlot <- renderPlot({
    
    x <- input$dataset
    df <- as.dist(m = read.csv(paste0("./data2/wasserstein_", x, ".csv"), header = F)[-1,-1])
    
    res <- eclust(df, FUNclust = input$method, k = input$nclust)
    clust <- res$cluster %>% as.factor()

    tsne <- Rtsne(df, perplexity = input$per, is_distance = T)
    tsne_plot <- data.frame(x = tsne$Y[,1], y = tsne$Y[,2])%>% as_tibble()
    colnames(tsne_plot) <- c("Dim.1", "Dim.2")
    tsne_plot <- tsne_plot  %>% mutate(groups = clust)
    
    ggscatter(tsne_plot, x = "Dim.1", y = "Dim.2",
              color = "groups",
              palette = 'pal5',
              label = as.character(seq(1,29,1)),
              size = 2)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
