library(shiny)
library(tidyr)
library(dplyr)
library(magrittr)
library(stringr)
library(scales)
library(ggplot2)


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
        tabPanel("Immunization Coverage",
                 plotOutput("coverage", height = "800px")),
        tabPanel("Hospitalization Time Series",
                 plotOutput("timeseries", height = "800px")),
        tabPanel("Total Hospitalizations",
                 tableOutput("table")),
        tabPanel("Difference from Counterfactual",
                 plotOutput("difference", height = "800px"))
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
    dat_new() %>%
      filter(scenario %in% input$checkGroup) %>%
      pivot_longer(cols = c("<6m":"All"), names_to = "Age", values_to = "value") %>%
      group_by(date, Age, scenario) %>%
      summarize(
        median = median(value),
        lower = quantile(value, probs = 0.025),
        upper = quantile(value, probs = 0.975),
        .groups = 'drop'  # Ensure ungrouping after summarization
      ) %>%
      mutate(Age = factor(Age, levels=c("<6m","6-11m","1-4yrs","5-64yrs","65-74yrs","75+yrs","All")))
  })

  diff_data <- reactive({
    req(dat_new(), input$checkGroup)  # Ensure data and input are available
    print(input$checkGroup)  # Debug: Print selected scenarios

    step1 = dat_new() %>%
      filter(scenario =="Counterfactual") %>%
      pivot_longer(cols = c("<6m":"All"), names_to = "Age", values_to = "value") %>%
      mutate(Counterfactual = value) %>%
      select(-scenario) %>%
      group_by(Age,sample) %>%
      summarize(Counterfactual = sum(Counterfactual))

    step2 = dat_new() %>%
      filter(scenario %in% input$checkGroup) %>%
      pivot_longer(cols = c("<6m":"All"), names_to = "Age", values_to = "value") %>%
      group_by(Age,scenario,sample) %>%
      summarize(total = sum(value)) %>%
      left_join(step1, by=c("Age","sample")) %>%
      mutate(diff = (Counterfactual - total)/Counterfactual*100) %>%
      group_by(Age, scenario) %>%
      summarize(
        median = median(diff),
        lower = quantile(diff, probs = 0.025),
        upper = quantile(diff, probs = 0.975),
        .groups = 'drop'  # Ensure ungrouping after summarization
      ) %>%
      mutate(Age = factor(Age, levels=c("<6m","6-11m","1-4yrs","5-64yrs","65-74yrs","75+yrs","All")))

    step2
  })


  # Filter data based on selected scenarios
  coverage_data <- reactive({
    req(dat_new(), input$checkGroup)  # Ensure data and input are available
    print(input$checkGroup)  # Debug: Print selected scenarios
    dat_new() %>%
      filter(scenario %in% input$checkGroup) %>%
      select(date:scenario) %>%
      pivot_longer(cols=c("monoclonal_birth":"adult_vax_75"),names_to = "immunization",values_to="doses") %>%
      group_by(date, immunization, scenario) %>%
      summarize(doses = unique(doses),
        .groups = 'drop'  # Ensure ungrouping after summarization
      ) %>%
      mutate(immunization = factor(immunization, levels=c("monoclonal_birth","monoclonal_catchup","maternal_vax","adult_vax_65to74","adult_vax_75"),
                                   labels=c("Monoclonal Birth Doses","Monoclonal Catchup Doses","Maternal Vaccinations","Adult Vaccinations (65-74yrs)","Adult Vaccinations (75+)")))
  })


  # Render the coverage plot
  output$coverage <- renderPlot({
    req(coverage_data())  # Ensure filtered data is available
    ggplot(data = coverage_data()) +
      theme_bw() +
      geom_line(aes(x = date, y = doses, group = scenario, color = scenario)) +
      scale_y_continuous(labels = label_comma(scale = 1))+
      facet_wrap(~immunization, ncol=3, scales = "free") +
      theme(legend.position = "top",
            axis.text = element_text(size = 15),
            strip.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 20),
            title = element_text(size = 20)) +
      labs(x = NULL, y = "Cumulative Immunization Doses", color = "Scenario")
  })


  # Render the timeseries plot
  output$timeseries <- renderPlot({
    req(filtered_data())  # Ensure filtered data is available
    ggplot(data = filtered_data()) +
      theme_bw() +
      geom_line(aes(x = date, y = median, group = scenario, color = scenario)) +
      geom_ribbon(aes(x = date, ymin = lower, ymax = upper, group = scenario, fill = scenario), alpha = 0.5) +
      facet_wrap(~Age, ncol = 2, scales = "free") +
      theme(legend.position = "top",
            axis.text = element_text(size = 15),
            strip.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 20),
            title = element_text(size = 20)) +
      guides(fill = "none") +
      labs(x = NULL, y = "RSV Hospitalizations", color = "Scenario")
  })



  # Render percent difference plot
  output$difference <- renderPlot({
    req(diff_data())  # Ensure filtered data is available
    ggplot(data = diff_data()) +
      theme_bw() +
      geom_point(aes(y = reorder(scenario,desc(scenario)), x = median, group = scenario, color = scenario),size=3) +
      geom_errorbar(aes(y=reorder(scenario,desc(scenario)), xmin = lower, xmax = upper, group = scenario, color = scenario), alpha = 0.5,width=.3,linewidth=2) +
      facet_wrap(~Age, ncol = 2, scales = "free") +
     # scale_y_reverse() +
      theme(legend.position = "top",
            axis.text = element_text(size = 15),
            strip.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 20),
            title = element_text(size = 20)) +
      labs(x = "% Difference", y = NULL, color = "Scenario")
  })


  # Render the summary table
  output$table <- renderTable({
    req(filtered_data())  # Ensure filtered data is available
    print(head(filtered_data()))  # Debug: Print the filtered data
    filtered_data() %>%
      group_by(Age, scenario) %>%
      summarise(Total = format(round(sum(median), 1),nsmall=1),
                lower = format(round(sum(lower), 1),nsmall=1),
                upper = format(round(sum(upper), 1),nsmall=1)) %>%
      ungroup() %>%
      mutate(combine = paste0(Total, " (", lower, "-", upper, ")"),
             combine = ifelse(grepl("NA", combine), str_sub(combine, end = -8), combine)) %>%
      select(-Total, -lower, -upper) %>%
      pivot_wider(names_from = scenario, values_from = combine)
  }, striped = TRUE, spacing = "s", align = "c", digits = 0, width = "auto")
}
shinyApp(ui = ui, server = server)
