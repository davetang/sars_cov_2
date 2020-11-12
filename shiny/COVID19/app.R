#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(plotly)
theme_set(theme_bw())
time_series_data <- read_csv(file = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv")

# Define UI for application that draws a histogram
ui <- fluidPage(

    titlePanel("COVID19 cases"),
    
    sidebarLayout(
        sidebarPanel(
            selectInput("country", "Country/Region:", choices = unique(time_series_data$`Country/Region`))
        ),

        mainPanel(
           plotlyOutput("case_plot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$case_plot <- renderPlotly({
        time_series_data %>%
            filter(`Country/Region` == input$country) %>%
            select_at(-(1:4)) %>%
            t() %>%
            as.data.frame() %>%
            rownames_to_column(var = "date") %>%
            mutate(confirmed = rowSums(.[-1])) %>%
            mutate(date = as.Date(date, format = "%m/%d/%y")) %>%
            mutate(confirmed = c(0, diff(confirmed))) -> df
        
        my_total <- sum(df$confirmed)
        my_title <- paste0(input$country, ": ", my_total)
        df$confirmed[df$confirmed < 0] <- 0
        
        ggplot(df, aes(date, confirmed)) +
            geom_col() +
            theme(axis.title.x = element_blank()) +
            labs(y = "Confirmed cases", title = my_title) -> p
        
        ggplotly(p)
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
