# Load necessary libraries
library(shiny)
library(plotly)
library(shinyjs)
columns <- c("condition", "baseline", "HLA", "proportion_train", "model","feature_selection_method","microbiome_type","weighted","sample_number","features","horizon_time")

tab_to_column_df = data.frame(tabs=c("Condition","Baseline-Age", "HLA","proportion-train","model","feature_selection_method","microbiome_type","weighted","Sample-Number","features","horizon_time"),columns)

parameter_tabs <- tabsetPanel(
  id = "params",
  type = "hidden",
  tabPanel("Condition",
           selectInput("Condition","What condition would you like to highlight?",choices=NULL)
  ),
  tabPanel("Baseline-Age", 
           selectInput("Baseline-Age","What are the ages of subjects you would like to highlight?",choices=NULL)
  ),
  tabPanel("Sample-Number",
           sliderInput("Sample-Number", label="Choose the maximum samle size of specifications you want to see.",value = 0,min=0,max=0)
  ),
  tabPanel("HLA", 
           selectInput("HLA","Choose subjects to highlight with specific HLA haplotypes.",choices=NULL)
  ),
  tabPanel("proportion-train", 
           selectInput("proportion-train","Highlight the specifications based on the proportion of training data.",choices=NULL)
  ),
  tabPanel("model", 
           selectInput("model","Highlight the specifications with a specific type of machine learning method",choices=NULL)
  ),
  tabPanel("feature_selection_method", 
           selectInput("feature_selection_method","Highlight the specifications with a specific type of feature selection method",choices=NULL)
  ),
  tabPanel("microbiome_type", 
           selectInput("microbiome_type","Highlight the specifications based on the type of microbial features used.",choices=NULL)
  ),
  tabPanel("weighted", 
           selectInput("weighted","Highlight the specifications that did or did weight data to account for unbalanced data.",choices=NULL)
  ),
  tabPanel("features", 
           selectInput("features","Highlight the specifications based on the clincal and microbiome data used for prediction.",choices=NULL)
  ),
  tabPanel("horizon_time", 
           selectInput("horizon_time","Highlight the specifications based on the number of years in the future prediction was done for",choices=NULL)
  ),
  
)

# Define the UI
ui <- fluidPage(
  # Sidebar with controls to select the input file
  sidebarLayout(
    sidebarPanel(
      selectInput("Column", "What feature would you like to filter by?", c("Condition","Baseline-Age","Sample-Number","HLA","proportion-train","model","feature_selection_method","microbiome_type","weighted","features","horizon_time"),selected=""),
      parameter_tabs,
    ),
    mainPanel(
      plotlyOutput("plot"),
      textOutput("my_textbox")
    )
  )
)

# Define the server
server <- function(input, output,session) {
  # allow file uploads of up to 10,000 MBs
  options(shiny.maxRequestSize=10000*1024^2)
  
  output$my_textbox <- renderText({
    "Hi my name is Sam Zimmerman. Here you can interact with the results from my study The tenuous role of the metagenome in predicting future T1D: results from a specification curve analysis on the TEDDY cohort. Each dot is a different specification ranked by the predictive ability of the microbiome on T1D and related phenotypes by AUC. Put your mouse on a dot to see the analytical choices we made for the specification. Also, on the tab to the left, you can select specifications to highlight based on the analytical choices or number of samples for the given specification. Please go to our github page to view the code used to perform this analysis. https://github.com/b-tierney/teddy Thank you for visiting! If you have any questions, email me, Sam Zimmerman at samuel.e.zimmerman@gmail.com. Have a wonderful day!"
  })
  
  # Reactive expression to read the selected file
  losso_cox_res <- reactive({
    #req(input$file1)
    #read.csv(input$file1$datapath)
    req("~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/lasso_rf_regression_output_df_with_AUCs.csv")
    data <- read.csv("~/Dropbox (HMS)/Kostic_Lab/datasets/TEDDY/TEDDY_analysis_v2_march_2022/lasso_rf_regression_output_df_with_AUCs.csv")
  })
  
  observeEvent(input$Column, {
    #max_sample_number <- max(losso_cox_res()$sample_number)
    updateTabsetPanel(session,"params",selected=input$Column)
    #updateNumericInput(session,"sample_number", value = max_sample_number)
  }) 

  observeEvent(input$params, {
    if(input$params == "Sample-Number") {
      max_sample_number <- max(losso_cox_res()$sample_number,na.rm = TRUE)
      updateSliderInput(session,"Sample-Number", value = 0,max=max_sample_number)
    } else if(input$params == "Condition") {
      updateSelectizeInput(session, 'Condition', choices = c("all",unique(losso_cox_res()$condition)), server = TRUE)
    } else if(input$params == "Baseline-Age") {
      updateSelectizeInput(session, 'Baseline-Age', choices = c("all",unique(losso_cox_res()$baseline)), server = TRUE)
    } else if(input$params == "HLA") {
      updateSelectizeInput(session, 'HLA', choices = c("all",unique(losso_cox_res()$HLA)), server = TRUE)
    } else if(input$params == "proportion-train") {
      updateSelectizeInput(session, 'proportion-train', choices = c("all",unique(losso_cox_res()$proportion_train)), server = TRUE) 
    } else if(input$params == "model") {
      updateSelectizeInput(session, 'model', choices = c("all",unique(losso_cox_res()$model)), server = TRUE) 
    } else if(input$params == "feature_selection_method") {
      updateSelectizeInput(session, 'feature_selection_method', choices = c("all",unique(losso_cox_res()$feature_selection_method)), server = TRUE) 
    } else if(input$params == "microbiome_type") {
      updateSelectizeInput(session, 'microbiome_type', choices = c("all",unique(losso_cox_res()$microbiome_type)), server = TRUE) 
    } else if(input$params == "weighted") {
      updateSelectizeInput(session, 'weighted', choices = c("all",unique(losso_cox_res()$weighted)), server = TRUE) 
    } else if(input$params == "features") {
      updateSelectizeInput(session, 'features', choices = c("all",unique(losso_cox_res()$features)), server = TRUE) 
    } else if(input$params == "horizon_time") {
      updateSelectizeInput(session, 'horizon_time', choices = c("all",unique(losso_cox_res()$horizon_time)), server = TRUE) 
    }
  })

  #sample_number <- reactive({
  #    filter(losso_cox_res, sample_number <= input$sample_number)
  #  })
  #observeEvent(sample_number(), {
  #  choices <- max(territory()$sample_number)
  #  updateSelectInput(inputId = "sample_number", choices = choices) 
  #})
  
  all_condition_plot <- function(df,column_name) {
    plot_ly(data = df, x = ~rank, y = ~AUC, color = as.formula(paste("~",column_name,sep="")),type = "scatter",
            text = ~paste("HLA: ", HLA, '<br>cases_horizon:', cases_time_horizon,'<br>survivors_horizon:',survivors_time_horizon,'<br>condition:',condition,'<br>baseline',baseline,'<br>horizon',horizon_time,'<br>total_samples:',sample_number, '<br>CI:',paste(round(lower_CI_AUC,4),round(upper_CI_AUC,4),sep="-"))) %>% layout(legend = list(orientation = "h", xanchor = "center", x = 0.5,y=-0.3))
  }
  
  specific_condition_plot <- function(df,chosen_condition,column_name) {
    df$temp = df[,column_name] == chosen_condition
    df$temp[df$temp==TRUE] = chosen_condition
    df$temp[df$temp==FALSE] = "other"
    df$temp = factor(df$temp,levels=c("other",chosen_condition))
    plot_ly(data = df, x = ~rank, y = ~AUC, color = ~temp,colors=c("#cccccc","#000000"),type = "scatter",
            text = ~paste("HLA: ", HLA, '<br>cases_horizon:', cases_time_horizon,'<br>survivors_horizon:',survivors_time_horizon,'<br>condition:',condition,'<br>baseline',baseline,'<br>horizon',horizon_time,'<br>total_samples:',sample_number, '<br>CI:',paste(round(lower_CI_AUC,4),round(upper_CI_AUC,4),sep="-"))) %>% layout(legend = list(orientation = "h", xanchor = "center", x = 0.5,y=-0.3))
  }
  
  # Reactive expression to create the plot
  output$plot <- renderPlotly({
    # Get the data
    data <- losso_cox_res()
    # Remove AUC NAs and rank the data
    data <- data[!is.na(data$AUC),]
    data <- data[order(data$AUC),]
    data$rank <- 1:nrow(data)

    # Create the plot
    df_column_name = tab_to_column_df[match(input$params,tab_to_column_df$tabs),"columns"]

    if(input$params == "Sample-Number") {
      data = data[data$sample_number >= input$`Sample-Number`,]
      all_condition_plot(data,df_column_name)
    } else if(input$params == "Condition") {
      if(input$Condition == "all") {
        all_condition_plot(data,df_column_name)
      } else {
        chosen_condition = input$Condition
        specific_condition_plot(data,chosen_condition,df_column_name)
      }
    } else if(input$params == "Baseline-Age") {
      if(input$`Baseline-Age` == "all") {
        all_condition_plot(data,df_column_name)
      } else {
        chosen_condition = input$`Baseline-Age`
        specific_condition_plot(data,chosen_condition,df_column_name)
      }
    } else if(input$params == "HLA") {
      if(input$HLA == "all") {
        all_condition_plot(data,df_column_name)
      } else {
        chosen_condition = input$HLA
        specific_condition_plot(data,chosen_condition,df_column_name)
      }
    } else if(input$params == "proportion-train") {
      if(input$`proportion-train` == "all") {
        all_condition_plot(data,df_column_name)
      } else {
        chosen_condition = input$`proportion-train`
        specific_condition_plot(data,chosen_condition,df_column_name)
      }
    } else if(input$params == "model") {
      if(input$model == "all") {
        all_condition_plot(data,df_column_name)
      } else {
        chosen_condition = input$model
        specific_condition_plot(data,chosen_condition,df_column_name)
      }
    } else if(input$params == "feature_selection_method") {
      if(input$feature_selection_method == "all") {
      all_condition_plot(data,df_column_name)
      } else {
        chosen_condition = input$feature_selection_method
        specific_condition_plot(data,chosen_condition,df_column_name)
      }
    } else if(input$params == "microbiome_type") {
      if(input$microbiome_type == "all") {
        all_condition_plot(data,df_column_name)
      } else {
        chosen_condition = input$microbiome_type
        specific_condition_plot(data,chosen_condition,df_column_name)
      }
    } else if(input$params == "weighted") {
      if(input$weighted == "all") {
        all_condition_plot(data,df_column_name)
      } else {
        chosen_condition = input$weighted
        specific_condition_plot(data,chosen_condition,df_column_name)
      } 
    } else if(input$params == "features") {
      if(input$features == "all") {
        all_condition_plot(data,df_column_name)
      } else {
        chosen_condition = input$features
        specific_condition_plot(data,chosen_condition,df_column_name)
      }
    } else if(input$params == "horizon_time") {
      if(input$horizon_time == "all") {
        all_condition_plot(data,df_column_name)
      } else {
        chosen_condition = input$horizon_time
        specific_condition_plot(data,chosen_condition,df_column_name)
      }
    }
  })
}

# Run the application
shinyApp(ui = ui, server = server)