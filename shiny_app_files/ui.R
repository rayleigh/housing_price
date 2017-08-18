library(shiny)
library(shinycssloaders)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Housing price in LA"),
  
  # Sidebar with a slider input for the number of bins
  sidebarLayout(
    sidebarPanel(
    uiOutput("year_slider"),
    numericInput("wait_time", "Wait Time (ms)", value = 3000),
    checkboxGroupInput("var_list", label = "Variables of interest",
                       choices = list("Median Income" = "medianincome_", "White Proportion" = "white_p_",
                                      "Black Proportion" = "black_p_", "Hispanic Proportion" = "hispanic_p_",
                                      "Asian Proportion" = "asian_p_", "DataQuick loan" = "loanvalue2010usd",
                                      "HMDA loan" = "medianloanamount_", "DataQuick count" = "numloans",
                                      "HMDA count" = "loancount_c_", "Home price" = "pricepersqft2010usd"),
                       selected = "medianincome_")),
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("map")
    )
  )
))