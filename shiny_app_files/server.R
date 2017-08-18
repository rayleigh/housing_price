library(ggplot2)
library(tigris)
library(dplyr)
library(ggmap)
library(gridExtra)

options(tigris_use_cache = TRUE)

full_data <- read.csv(file = "Fulldataset_long.csv", header = T)
la_tract_info <- tracts(state = "06", county = "037", year = 2010)
la_fortify <- fortify(la_tract_info, region = "TRACTCE10")
la_fortify <- mutate(la_fortify, trtid10 = as.numeric(paste("6037", id, sep = "")))
plot_data <- left_join(la_fortify, full_data, by = "trtid10")

shinyServer(
  function(input, output) {
    output$year_slider <- renderUI({
      sliderInput("year",
                  "Year:",
                  min = 1990,
                  max = 2010,
                  value = 1990, 
                  sep = "",
                  animate=animationOptions(interval=input$wait_time, loop=T))
    })
    
    output$map <- renderPlot({
      tmp <- plot_data[plot_data$year == input$year, c("long", "lat", "group", input$var_list)]
      plots <- lapply(input$var_list, function(var) {
        ggplot(tmp) + geom_polygon(aes_string(x = "long", y = "lat", group = "group", fill = var), color = "black", size = 0.25) +
          scale_fill_continuous(na.value="white")
      })
      do.call("grid.arrange", c(plots, ncol = ceiling(sqrt(length(plots)))))
    }, height = 600, width = 800)
  })
