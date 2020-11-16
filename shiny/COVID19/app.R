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
library(DT)
theme_set(theme_bw())

confirmed_csv <- "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv"
read_csv(file = confirmed_csv) %>%
    mutate(`Country/Region` = sub(pattern = "US",
                                  replacement = "United States",
                                  x = `Country/Region`)
           ) -> time_series_confirmed


deaths_csv <- "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv"
read_csv(file = deaths_csv) %>%
    mutate(`Country/Region` = sub(pattern = "US",
                                  replacement = "United States",
                                  x = `Country/Region`)
           ) -> time_series_deaths

choices <- unique(time_series_confirmed$`Country/Region`)

prep_table <- function(country){
    
    time_series_deaths %>%
        filter(`Country/Region` == country) %>%
        select_at(-(1:4)) %>%
        t() %>%
        as.data.frame() %>%
        rownames_to_column(var = "date") %>%
        mutate(deaths = rowSums(.[-1])) %>%
        mutate(date = as.Date(date, format = "%m/%d/%y")) %>%
        mutate(deaths = c(0, diff(deaths))) %>%
        select(date, deaths) %>%
        arrange(desc(date)) -> deaths
    deaths$deaths[deaths$deaths < 0] <- 0
    
    time_series_confirmed %>%
        filter(`Country/Region` == country) %>%
        select_at(-(1:4)) %>%
        t() %>%
        as.data.frame() %>%
        rownames_to_column(var = "date") %>%
        mutate(confirmed = rowSums(.[-1])) %>%
        mutate(date = as.Date(date, format = "%m/%d/%y")) %>%
        mutate(confirmed = c(0, diff(confirmed))) %>%
        select(date, confirmed) %>%
        arrange(desc(date)) %>%
        inner_join(y = deaths, by = "date") %>%
        mutate(day = strftime(date, format = "%A")) %>%
        select(date, day, everything()) -> confirmed
    confirmed$confirmed[confirmed$confirmed < 0] <- 0
    return(confirmed)
}

# Define UI for application that draws a histogram
ui <- fluidPage(

    titlePanel("COVID-19 statistics"),
    
    sidebarLayout(
        sidebarPanel(
            selectInput("country", "Country/Region:", choices = choices, selected = "Sweden"),
            a("Data source", href="https://github.com/CSSEGISandData/COVID-19"),
            br(),
            a("Source code", href="https://github.com/davetang/sars_cov_2/tree/master/shiny/COVID19")
        ),

        mainPanel(
            tabsetPanel(type = "tabs",
                        tabPanel("Plots",
                                 plotlyOutput("plot_confirmed"),
                                 plotlyOutput("plot_deaths")
                        ),
                        tabPanel("Table",
                                 dataTableOutput("table"))
            )
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$plot_confirmed <- renderPlotly({
        time_series_confirmed %>%
            filter(`Country/Region` == input$country) %>%
            select_at(-(1:4)) %>%
            t() %>%
            as.data.frame() %>%
            rownames_to_column(var = "date") %>%
            mutate(confirmed = rowSums(.[-1])) %>%
            mutate(date = as.Date(date, format = "%m/%d/%y")) %>%
            mutate(confirmed = c(0, diff(confirmed))) -> df
        
        my_total <- sum(df$confirmed)
        my_title <- paste0(input$country, " total: ", my_total)
        df$confirmed[df$confirmed < 0] <- 0
        
        ggplot(df, aes(date, confirmed)) +
            geom_col() +
            theme(axis.title.x = element_blank()) +
            labs(y = "Confirmed cases", title = my_title) -> p
        
        ggplotly(p)
    })
    
    output$plot_deaths <- renderPlotly({
        time_series_deaths %>%
            filter(`Country/Region` == input$country) %>%
            select_at(-(1:4)) %>%
            t() %>%
            as.data.frame() %>%
            rownames_to_column(var = "date") %>%
            mutate(deaths = rowSums(.[-1])) %>%
            mutate(date = as.Date(date, format = "%m/%d/%y")) %>%
            mutate(deaths = c(0, diff(deaths))) -> df
        
        my_total <- sum(df$deaths)
        my_title <- paste0(input$country, " total: ", my_total)
        df$deaths[df$deaths < 0] <- 0
        
        ggplot(df, aes(date, deaths)) +
            geom_col() +
            theme(axis.title.x = element_blank()) +
            labs(y = "Confirmed deaths", title = my_title) -> p
        
        ggplotly(p)
    })
    
    output$table <- renderDataTable(
        datatable(prep_table(input$country))
    )
    
}

# Run the application 
shinyApp(ui = ui, server = server)
