

library(shiny)
library(bslib)
library(bsicons)


# Define UI for application that draws a histogram
ui <- page_fillable(
    
    div(
        style = "padding: 10px 16px;",
        h2("Temperature Conversion")
    ),
    
    input_dark_mode(),

    layout_column_wrap(
        width = 1/2, height = 300,
        
        card(
            card_header("Input Temperature"),
          numericInput(
              "input",
              label = "Input the value",
              value = 0
          ),
          
          radioButtons(
              "inputt",
              label = "Select the scale",
              choiceNames = c(
                  "Celsius (C)",
                  "Fahrenheit",
                  "Kelvin"
              ),
              choiceValues = c("C", "F", "K") # The values passed to input$opt
          ),
          
          br(),
          
          actionButton(
              "cal",
              "Calculate Temperature")
          
        ),
       
        
        card(
            card_header("Output Temperature"),
            # numericInput(
            #     "output",
            #     label = "Input the value",
            #     value = 0
            # ),
            
            value_box(
                title = "Temperature",
                value = textOutput("result"),
                showcase = bs_icon("thermometer")
                
            ),
            
            radioButtons(
                "out_opt",
                label = "Select the scale",
                choiceNames = c(
                    "Celsius (C)",
                    "Fahrenheit",
                    "Kelvin"
                ),
                choiceValues = c("C", "F", "K") # The values passed to input$opt
            )
            
        )
    )
)




# Define server logic required to draw a histogram
server <- function(input, output) {

    converted <- reactive({
        val <- input$input
        from <- input$inputt
        to <- input$out_opt
        
        
        celsius <- switch(from,
            
            "C" = val, 
            "F" = ((val - 32) * 5/9) , 
            "K" = val - 273.15
        )
        
        switch (to,
            "C" = celsius, 
            "F" = (celsius * (9 / 5)) + 32 , 
            "K" = celsius + 273.15
    ) 
        
        
    })|> bindEvent(input$cal)


    output$result <- renderText({
        paste(round(converted(), 2), input$out_opt)
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
