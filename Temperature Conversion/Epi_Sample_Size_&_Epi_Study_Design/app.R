
library(shiny)
library(bslib)
library(shinydashboard)

# Define UI for application that draws a histogram
ui <- page_fillable(
    # dashboardHeader(title = "Epidemiological Study Designer"),
    
    layout_column_wrap(
        width = 1/3,
        card(
            card_header("Questions"),
            
            
            selectInput("obj",
                        label = "What is your study objective?",
                        choices = c(
                            "Compare",
                            "Describe"
                        ),
                        
                        multiple = F
            ),
            
            
            selectInput(
                "time",
                label = "What is your time frame?",
                choices = c( 
                    "Prospective",
                    "Retrospective"
                ),
                multiple = F
            ),
            
            
            
            selectInput(
                "s_type",
                label = "What is your study type?",
                choices = c( 
                    "Observational",
                    "Experimental"
                ),
                multiple = F
            ),
            
            
            
            selectInput(
                "o_type",
                label = "What is your outcome type?",
                choices = c( 
                    "Continuous",
                    "Binary"
                ),
                multiple = F
            ),
            
            selectInput(
                "group",
                label = "How many groups?",
                choices = c( 
                    "Single",
                    "Multiple"
                ),
                multiple = F
            ),
            actionButton("run",
                         "Run")
            
        ),
        
        
        card(
            card_header("Recommended Study Design"),
            card_body(textOutput("epi_s_d"))
        ),
        
        
        card(
            card_header("Sample Size Formula "),
            card_body(textOutput("sample_si"))
        )
        
    )
    
    
    
)






# Define server logic required to draw a histogram
server <- function(input, output) {
    
    study_design <- reactive({
        
        req(
            input$obj,
            input$time
        )
        
        switch(input$obj,
               "Compare" = {
                   switch(input$time,
                          "Prospective" = "Randomised Controlled Trial",
                          "Retrospective" = "Case Control Study",
                          "Design not in database"
                   )
                   
               },
               
               
               "Describe" = {
                   switch(input$time,
                          "Prospective" = "Cohort Study",
                          "Retrospective" = "Cross Sectional Study",
                          "Design not in database"
                          
                   )
                   
               },
               
               "Design not in database"
               
        )
        
    }) |> bindEvent(input$run)
    
    output$epi_s_d <- renderText(
        paste(study_design())
        
    )
    
    
    
    
    sample_size <- reactive({
        
        req(
            input$s_type,
            input$o_type,
            input$group)
        
        switch(input$s_type,
               
               "Observational" = {
                   switch(input$o_type,
                          "Continuous" = {
                              switch(input$group,
                                     "Single" = "One sample mean formula",
                                     "Multiple" = "Two sample mean formula",
                                     "Not in database"
                                     
                              )
                          },
                          
                          
                          "Binary" = {
                              switch(input$group,
                                     "Single" = "One sample proportion formula",
                                     "Multiple" = "Two sample proportion formula",
                                     "Not in database"
                                     
                              )
                          }
                          
                   )
               },
               
               
               "Experimental" = {
                   switch(input$o_type,
                          "Continuous" = {
                              switch(input$group,
                                     "Single" = "One sample RCT formula",
                                     "Multiple" = "Two sample RCT formula",
                                     "Not in database"
                                     
                              )
                          },
                          
                          
                          "Binary" = {
                              switch(input$group,
                                     "Single" = "One sample RCT proportion",
                                     "Multiple" = "Two sample RCT proportion",
                                     "Not in database"
                                     
                              )
                          }
                          
                   )
               }
               
               
        )
        
    }) |> bindEvent(input$run)
    
    
    
    output$sample_si <- renderText(
        paste(sample_size())
        
    )
    
    
}

# Run the application 
shinyApp(ui = ui, server = server)
