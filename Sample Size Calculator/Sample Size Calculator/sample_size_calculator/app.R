library(bslib)
library(bsicons)
library(shiny)
library(DT)


# ── Z value helpers ──────────────────────────────────────────────
z_alpha <- function(alpha) qnorm(1 - alpha / 2)
z_beta  <- function(power) qnorm(power)

# ── Sample Size Functions ────────────────────────────────────────

# 1. Prevalence
ss_prevalence <- function(p, precision, alpha) {
  za <- z_alpha(alpha)
  n  <- (za^2 * p * (1 - p)) / precision^2
  ceiling(n)
}

# 2. One Sample Mean
ss_one_mean <- function(mean0, mean1, sd, alpha, power) {
  za    <- z_alpha(alpha)
  zb    <- z_beta(power)
  delta <- abs(mean1 - mean0)
  n     <- ((za + zb) * sd / delta)^2
  ceiling(n)
}

# 3. Two Sample Means Equal SD
ss_two_means_equal <- function(mean1, mean2, sd, alpha, power) {
  za    <- z_alpha(alpha)
  zb    <- z_beta(power)
  delta <- abs(mean1 - mean2)
  n     <- 2 * ((za + zb) * sd / delta)^2
  ceiling(n)
}

# 4. Two Sample Means Unequal SD
ss_two_means_unequal <- function(mean1, mean2, sd1, sd2, alpha, power) {
  za    <- z_alpha(alpha)
  zb    <- z_beta(power)
  delta <- abs(mean1 - mean2)
  n     <- ((za + zb)^2 * (sd1^2 + sd2^2)) / delta^2
  ceiling(n)
}

# 5. Paired Samples
ss_paired <- function(mean1, mean2, sd_diff, alpha, power) {
  za    <- z_alpha(alpha)
  zb    <- z_beta(power)
  delta <- abs(mean1 - mean2)
  n     <- ((za + zb) * sd_diff / delta)^2
  ceiling(n)
}

# 6. Two Sample Proportions
ss_two_proportions <- function(p1, p2, alpha, power) {
  za    <- z_alpha(alpha)
  zb    <- z_beta(power)
  delta <- abs(p1 - p2)
  n     <- (za + zb)^2 * (p1*(1-p1) + p2*(1-p2)) / delta^2
  ceiling(n)
}

# 7. Logistic Regression
ss_logistic <- function(p1, OR, alpha, power) {
  za <- z_alpha(alpha)
  zb <- z_beta(power)
  n  <- (za + zb)^2 / (p1 * (1 - p1) * log(OR)^2)
  ceiling(n)
}

# 8. Relative Risk
ss_relative_risk <- function(p1, RR, alpha, power) {
  za    <- z_alpha(alpha)
  zb    <- z_beta(power)
  p2    <- p1 * RR
  delta <- abs(p1 - p2)
  n     <- (za + zb)^2 * (p1*(1-p1) + p2*(1-p2)) / delta^2
  ceiling(n)
}

# 9. Poisson Regression
ss_poisson <- function(rate1, rate2, alpha, power) {
  za <- z_alpha(alpha)
  zb <- z_beta(power)
  n  <- (za + zb)^2 * (1/rate1 + 1/rate2) / log(rate1/rate2)^2
  ceiling(n)
}

# 10. Survival Analysis
ss_survival <- function(S1, S2, alpha, power) {
  za <- z_alpha(alpha)
  zb <- z_beta(power)
  n  <- 2 * (za + zb)^2 / log(S1/S2)^2
  ceiling(n)
}


# Define UI for application that draws a histogram
library(shiny)
library(bslib)

ui <- page_sidebar(
    title = "Statistical Sample Size  Calculator",
    fillable = T,
  

    
    sidebar = sidebar(
        open = TRUE,
        width = 260,
        
        
        
        selectInput(
            "study_type",
            label    = "Select Study Type",
            choices  = c(
                "Prevalence Study"                  = "prevalence",
                "One Sample Mean"                   = "one_mean",
                "Two Sample Means (Equal SD)"       = "two_mean_equal",
                "Two Sample Means (Unequal SD)"     = "two_mean_unequal",
                "Paired Samples"                    = "paired",
                "Two Sample Proportions"            = "two_prop",
                "Logistic Regression"               = "logistic",
                "Relative Risk"                     = "rr",
                "Poisson Regression"                = "poisson",
                "Survival Analysis"                 = "survival"
            )
        ),
        
        
        # hr(),
        # h5("Parameters", style = "color:#2C7A7B; font-weight:bold;"),
        
        # ── Dynamic inputs ──
        uiOutput("dynamic_inputs"),
        
        hr(),
        h5("Test Parameters", style = "color:#2C7A7B; font-weight:bold;"),
        
        # Alpha slider (not shown for prevalence)
        conditionalPanel(
            condition = "input.study_type != 'prevalence'",
            sliderInput("power", "Power (1 - β)",
                        min = 0.70, 
                        max = 0.99,
                        value = 0.80, 
                        step = 0.01,
                        ticks = F)
        ),
        
        sliderInput("alpha", "Significance Level (α)",
                    ticks = F,
                    min = 0.01, 
                    max = 0.10,
                    value = 0.05, 
                    step = 0.01
                    ),
        
        hr(),
        actionButton("calculate", "Calculate Sample Size",
                     class = "btn-primary w-100",
                     icon  = icon("calculator")),
        
        
    ),
    
    
    
    # ── Main Panel ──
    layout_column_wrap(
        width = 1,

        # Results cards
        layout_column_wrap(
            width = 1/3,
            max_height = 120,
            
            value_box(
                title    = "Total Sample Size",
                value    = textOutput("total_n"),
                showcase = bsicons::bs_icon("universal-access")
            ),
            
            value_box(
                title    = "Sample Size Per Group",
                value    = textOutput("per_group_n"),
                showcase = bsicons::bs_icon("ui-radios-grid")
            ),
            
            value_box(
                title    = "Study Type",
                value    = textOutput("study_label"),
                showcase = bsicons::bs_icon("clock-history")
            )
        ),
        
        # Formula and interpretation
        layout_column_wrap(
            width = 1/2,
            min_height  = 100,
  
            
            card(
                card_header(
                    bsicons::bs_icon("file"), " Formula Used",
                    style = "background:#2C7A7B; color:white; font-weight:bold;"
                ),
                card_body(
                    uiOutput("formula_display")
                ),
                full_screen = T
            ),
            
            card(
                card_header(
                    bsicons::bs_icon("info-circle"), " Interpretation",
                    style = "background:#2C7A7B; color:white; font-weight:bold;"
                ),
                card_body(
                    uiOutput("interpretation")
                ),
                full_screen = T
            )
        ),
        
        # Parameters used
        card(
            card_header(
                bsicons::bs_icon("list-check"), " Parameters Used",
                style = "background:#2C7A7B; color:white; font-weight:bold;"
            ),
            card_body(
                tableOutput("params_table")
            ),
            full_screen = T
        )
    )
    
    
)

# ── Server ───────────────────────────────────────────────────────
server <- function(input, output, session) {
  
  # ── Dynamic UI inputs based on study type ──
  output$dynamic_inputs <- renderUI({
    switch(input$study_type,
           
           "prevalence" = tagList(
             sliderInput("p",         "Expected Prevalence (p)",  0.01, 0.99, 0.16, 0.01),
             sliderInput("precision", "Precision / Margin of Error (e)", 0.01, 0.20, 0.05, 0.01)
           ),
           
           "one_mean" = tagList(
             numericInput("mean0", "Null Mean (μ₀)",        value = 514),
             numericInput("mean1", "Alternative Mean (μ₁)", value = 534),
             numericInput("sd",    "Standard Deviation (σ)", value = 117)
           ),
           
           "two_mean_equal" = tagList(
             numericInput("mean1", "Mean Group 1 (μ₁)", value = 120),
             numericInput("mean2", "Mean Group 2 (μ₂)", value = 130),
             numericInput("sd",    "Common SD (σ)",      value = 8.5)
           ),
           
           "two_mean_unequal" = tagList(
             numericInput("mean1", "Mean Group 1 (μ₁)", value = 120),
             numericInput("mean2", "Mean Group 2 (μ₂)", value = 130),
             numericInput("sd1",   "SD Group 1 (σ₁)",   value = 8.5),
             numericInput("sd2",   "SD Group 2 (σ₂)",   value = 9.2)
           ),
           
           "paired" = tagList(
             numericInput("mean1",   "Mean Before (μ₁)",       value = 259),
             numericInput("mean2",   "Mean After (μ₂)",        value = 283),
             numericInput("sd_diff", "SD of Differences (σ_d)", value = 129)
           ),
           
           "two_prop" = tagList(
             sliderInput("p1", "Proportion Group 1 (p₁)", 0.01, 0.99, 0.19, 0.01),
             sliderInput("p2", "Proportion Group 2 (p₂)", 0.01, 0.99, 0.15, 0.01)
           ),
           
           "logistic" = tagList(
             sliderInput("p1", "Baseline Proportion (p₁)", 0.01, 0.99, 0.30, 0.01),
             numericInput("OR", "Odds Ratio (OR)",          value = 2.0)
           ),
           
           "rr" = tagList(
             sliderInput("p1", "Baseline Proportion (p₁)", 0.01, 0.99, 0.20, 0.01),
             numericInput("RR", "Relative Risk (RR)",       value = 1.5)
           ),
           
           "poisson" = tagList(
             numericInput("rate1", "Rate Group 1 (λ₁)", value = 0.10),
             numericInput("rate2", "Rate Group 2 (λ₂)", value = 0.20)
           ),
           
           "survival" = tagList(
             sliderInput("S1", "Survival Group 1 (S₁)", 0.01, 0.99, 0.60, 0.01),
             sliderInput("S2", "Survival Group 2 (S₂)", 0.01, 0.99, 0.40, 0.01)
           )
    )
  })
  
  # ── Reactive calculation ──
  results <- reactive({
    
    req(input$study_type, input$alpha)
    
    alpha <- input$alpha
    power <- if (input$study_type != "prevalence") input$power else NULL
    
    switch(input$study_type,
           
           "prevalence" = {
             req(input$p, input$precision)
             n     <- ss_prevalence(input$p, input$precision, alpha)
             label <- "Prevalence Study"
             formula <- "n = (Z²_α × p × (1-p)) / e²"
             list(total = n, per_group = n, label = label, formula = formula,
                  params = data.frame(
                    Parameter   = c("Prevalence (p)", "Precision (e)", "Alpha (α)", "Z_α"),
                    Value       = c(input$p, input$precision, alpha,
                                    round(z_alpha(alpha), 3))
                  ))
           },
           
           "one_mean" = {
             req(input$mean0, input$mean1, input$sd)
             n     <- ss_one_mean(input$mean0, input$mean1, input$sd, alpha, power)
             label <- "One Sample Mean"
             formula <- "n = ((Z_α + Z_β) × σ / δ)²"
             list(total = n, per_group = n, label = label, formula = formula,
                  params = data.frame(
                    Parameter = c("Null Mean (μ₀)", "Alt Mean (μ₁)", "SD (σ)",
                                  "Delta (δ)", "Alpha (α)", "Power (1-β)",
                                  "Z_α", "Z_β"),
                    Value     = c(input$mean0, input$mean1, input$sd,
                                  abs(input$mean1 - input$mean0),
                                  alpha, power,
                                  round(z_alpha(alpha), 3), round(z_beta(power), 3))
                  ))
           },
           
           "two_mean_equal" = {
             req(input$mean1, input$mean2, input$sd)
             n     <- ss_two_means_equal(input$mean1, input$mean2, input$sd, alpha, power)
             label <- "Two Sample Means (Equal SD)"
             formula <- "n = 2 × ((Z_α + Z_β) × σ / δ)²   per group"
             list(total = n * 2, per_group = n, label = label, formula = formula,
                  params = data.frame(
                    Parameter = c("Mean Group 1", "Mean Group 2", "Common SD (σ)",
                                  "Delta (δ)", "Alpha (α)", "Power (1-β)",
                                  "Z_α", "Z_β"),
                    Value     = c(input$mean1, input$mean2, input$sd,
                                  abs(input$mean1 - input$mean2),
                                  alpha, power,
                                  round(z_alpha(alpha), 3), round(z_beta(power), 3))
                  ))
           },
           
           "two_mean_unequal" = {
             req(input$mean1, input$mean2, input$sd1, input$sd2)
             n     <- ss_two_means_unequal(input$mean1, input$mean2,
                                           input$sd1,   input$sd2, alpha, power)
             label <- "Two Sample Means (Unequal SD)"
             formula <- "n = (Z_α + Z_β)² × (σ₁² + σ₂²) / δ²   per group"
             list(total = n * 2, per_group = n, label = label, formula = formula,
                  params = data.frame(
                    Parameter = c("Mean Group 1", "Mean Group 2",
                                  "SD Group 1 (σ₁)", "SD Group 2 (σ₂)",
                                  "Delta (δ)", "Alpha (α)", "Power (1-β)",
                                  "Z_α", "Z_β"),
                    Value     = c(input$mean1, input$mean2,
                                  input$sd1, input$sd2,
                                  abs(input$mean1 - input$mean2),
                                  alpha, power,
                                  round(z_alpha(alpha), 3), round(z_beta(power), 3))
                  ))
           },
           
           "paired" = {
             req(input$mean1, input$mean2, input$sd_diff)
             n     <- ss_paired(input$mean1, input$mean2, input$sd_diff, alpha, power)
             label <- "Paired Samples"
             formula <- "n = ((Z_α + Z_β) × σ_d / δ)²"
             list(total = n, per_group = n, label = label, formula = formula,
                  params = data.frame(
                    Parameter = c("Mean Before", "Mean After",
                                  "SD of Differences (σ_d)", "Delta (δ)",
                                  "Alpha (α)", "Power (1-β)", "Z_α", "Z_β"),
                    Value     = c(input$mean1, input$mean2,
                                  input$sd_diff, abs(input$mean1 - input$mean2),
                                  alpha, power,
                                  round(z_alpha(alpha), 3), round(z_beta(power), 3))
                  ))
           },
           
           "two_prop" = {
             req(input$p1, input$p2)
             n     <- ss_two_proportions(input$p1, input$p2, alpha, power)
             label <- "Two Sample Proportions"
             formula <- "n = (Z_α + Z_β)² × (p₁(1-p₁) + p₂(1-p₂)) / δ²   per group"
             list(total = n * 2, per_group = n, label = label, formula = formula,
                  params = data.frame(
                    Parameter = c("Proportion 1 (p₁)", "Proportion 2 (p₂)",
                                  "Delta (δ)", "Alpha (α)", "Power (1-β)",
                                  "Z_α", "Z_β"),
                    Value     = c(input$p1, input$p2,
                                  abs(input$p1 - input$p2),
                                  alpha, power,
                                  round(z_alpha(alpha), 3), round(z_beta(power), 3))
                  ))
           },
           
           "logistic" = {
             req(input$p1, input$OR)
             n     <- ss_logistic(input$p1, input$OR, alpha, power)
             label <- "Logistic Regression"
             formula <- "n = (Z_α + Z_β)² / (p₁ × (1-p₁) × ln(OR)²)"
             list(total = n, per_group = n, label = label, formula = formula,
                  params = data.frame(
                    Parameter = c("Baseline Proportion (p₁)", "Odds Ratio (OR)",
                                  "ln(OR)", "Alpha (α)", "Power (1-β)",
                                  "Z_α", "Z_β"),
                    Value     = c(input$p1, input$OR,
                                  round(log(input$OR), 3),
                                  alpha, power,
                                  round(z_alpha(alpha), 3), round(z_beta(power), 3))
                  ))
           },
           
           "rr" = {
             req(input$p1, input$RR)
             n     <- ss_relative_risk(input$p1, input$RR, alpha, power)
             p2    <- input$p1 * input$RR
             label <- "Relative Risk"
             formula <- "p₂ = p₁ × RR\nn = (Z_α + Z_β)² × (p₁(1-p₁) + p₂(1-p₂)) / δ²   per group"
             list(total = n * 2, per_group = n, label = label, formula = formula,
                  params = data.frame(
                    Parameter = c("Baseline Proportion (p₁)", "Relative Risk (RR)",
                                  "Derived p₂", "Delta (δ)",
                                  "Alpha (α)", "Power (1-β)", "Z_α", "Z_β"),
                    Value     = c(input$p1, input$RR,
                                  round(p2, 3), round(abs(input$p1 - p2), 3),
                                  alpha, power,
                                  round(z_alpha(alpha), 3), round(z_beta(power), 3))
                  ))
           },
           
           "poisson" = {
             req(input$rate1, input$rate2)
             n     <- ss_poisson(input$rate1, input$rate2, alpha, power)
             label <- "Poisson Regression"
             formula <- "n = (Z_α + Z_β)² × (1/λ₁ + 1/λ₂) / ln(λ₁/λ₂)²   per group"
             list(total = n * 2, per_group = n, label = label, formula = formula,
                  params = data.frame(
                    Parameter = c("Rate Group 1 (λ₁)", "Rate Group 2 (λ₂)",
                                  "ln(λ₁/λ₂)", "Alpha (α)", "Power (1-β)",
                                  "Z_α", "Z_β"),
                    Value     = c(input$rate1, input$rate2,
                                  round(log(input$rate1 / input$rate2), 3),
                                  alpha, power,
                                  round(z_alpha(alpha), 3), round(z_beta(power), 3))
                  ))
           },
           
           "survival" = {
             req(input$S1, input$S2)
             n     <- ss_survival(input$S1, input$S2, alpha, power)
             label <- "Survival Analysis"
             formula <- "n = 2 × (Z_α + Z_β)² / ln(S₁/S₂)²   per group"
             list(total = n * 2, per_group = n, label = label, formula = formula,
                  params = data.frame(
                    Parameter = c("Survival Group 1 (S₁)", "Survival Group 2 (S₂)",
                                  "ln(S₁/S₂)", "Alpha (α)", "Power (1-β)",
                                  "Z_α", "Z_β"),
                    Value     = c(input$S1, input$S2,
                                  round(log(input$S1 / input$S2), 3),
                                  alpha, power,
                                  round(z_alpha(alpha), 3), round(z_beta(power), 3))
                  ))
           }
    )
    
  }) |> bindEvent(input$calculate)
  
  # ── Outputs ──
  output$total_n <- renderText({
    req(results())
    format(results()$total, big.mark = ",")
  })
  
  output$per_group_n <- renderText({
    req(results())
    format(results()$per_group, big.mark = ",")
  })
  
  output$study_label <- renderText({
    req(results())
    results()$label
  })
  
  output$formula_display <- renderUI({
    req(results())
    tagList(
      tags$p(
        style = "font-family: 'Courier New', monospace;
                 font-size: 1.1em;
                 background: #f0f4f4;
                 padding: 12px;
                 border-left: 4px solid #2C7A7B;
                 border-radius: 4px;",
        results()$formula
      ),
      tags$p(
        style = "color: #666; font-size:0.9em;",
        "Based on Prof. Musonda's notes — University of Zambia"
      )
    )
  })
  
  output$interpretation <- renderUI({
    req(results())
    r <- results()
    
    # Interpretation text
    interp <- if (r$total == r$per_group) {
      paste0(
        "You need a total of ", format(r$total, big.mark = ","),
        " participants for this ", r$label, " study.",
        " This is a single group study so all ",
        format(r$total, big.mark = ","), " participants are in one group."
      )
    } else {
      paste0(
        "You need ", format(r$per_group, big.mark = ","),
        " participants per group, giving a total of ",
        format(r$total, big.mark = ","), " participants across both groups."
      )
    }
    
    tags$p(
      style = "font-size: 1em; line-height: 1.6;",
      interp
    )
  })
  
  output$params_table <- renderTable({
    req(results())
    results()$params
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
  
}

shinyApp(ui = ui, server = server)


# Run the application 
shinyApp(ui = ui, server = server)
