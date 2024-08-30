library(shiny)
library(tidyr)
library(dplyr)
library(magrittr)
library(stringr)

ui <- fluidPage(
  titlePanel("RSV Scenario Projections"),
  sidebarLayout(
    sidebarPanel(
      fileInput(inputId = "file", label = "Upload Data", accept = c(".rds")),
      checkboxGroupInput(
        inputId = "checkGroup",
        label = tags$h3("Select Scenarios"),
        choices = NULL,  # Initially set to NULL, will be updated in server
        selected = "Counterfactual"
      )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Time Series",
                 plotOutput("timeseries", height = "800px")),
        tabPanel("Summary Table",
                 tableOutput("table"))
      )
    )
  )
)
server <- function(input, output, session) {

  # Define the reactive expression to read and return data
  dat_new <- reactive({
    req(input$file)  # Ensure a file is uploaded
    data1 <- readRDS(input$file$datapath)
    print(head(data1))  # Debug: Print the first few rows
    data1$scenario <- as.factor(data1$scenario)  # Ensure 'scenario' is a factor
    data1
  })

  # Update the checkbox group based on the 'scenario' column
  observe({
    data1 <- dat_new()
    print(levels(data1$scenario))  # Debug: Check levels of 'scenario'
    updateCheckboxGroupInput(session, "checkGroup",
                             label = "Select Scenarios",  # Simple string label
                             choices = levels(data1$scenario),
                             selected = "Counterfactual")
  })

  # Filter data based on selected scenarios
  filtered_data <- reactive({
    req(dat_new(), input$checkGroup)  # Ensure data and input are available
    print(input$checkGroup)  # Debug: Print selected scenarios
    filter(dat_new(), scenario %in% input$checkGroup)
  })

  # Render the timeseries plot
  output$timeseries <- renderPlot({
    req(filtered_data())  # Ensure filtered data is available
    ggplot(data = filtered_data()) +
      theme_bw() +
      geom_line(aes(x = date, y = hosp, group = scenario, color = scenario)) +
      geom_ribbon(aes(x = date, ymin = lower, ymax = upper, group = scenario, fill = scenario), alpha = 0.5) +
      facet_wrap(~Age, ncol = 3, scales = "free") +
      theme(legend.position = "top",
            axis.text = element_text(size = 15),
            strip.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 20),
            title = element_text(size = 20)) +
      guides(fill = FALSE) +
      labs(x = NULL, y = "Hospitalization rate per 10,000", color = "Scenario",
           title = "Weekly RSV hospitalization rate per 10,000 by scenario")
  })

  # Render the summary table
  output$table <- renderTable({
    req(filtered_data())  # Ensure filtered data is available
    print(head(filtered_data()))  # Debug: Print the filtered data
    filtered_data() %>%
      group_by(Age, scenario) %>%
      summarise(Total = round(sum(hosp), 1),
                lower = round(sum(lower), 1),
                upper = round(sum(upper), 1)) %>%
      ungroup() %>%
      mutate(combine = paste0(Total, " (", lower, "-", upper, ")"),
             combine = ifelse(grepl("NA", combine), str_sub(combine, end = -8), combine)) %>%
      select(-Total, -lower, -upper) %>%
      pivot_wider(names_from = scenario, values_from = combine)
  }, striped = TRUE, spacing = "s", align = "c", digits = 0, width = "auto")
}
shinyApp(ui = ui, server = server)
