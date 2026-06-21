library(bslib)
library(bsicons)
library(shiny)

# ── Z value helpers ──────────────────────────────────────────────
z_alpha <- function(alpha) qnorm(1 - alpha / 2)
z_beta  <- function(power) qnorm(power)

# ── Sample Size Functions (ALL return n PER GROUP) ───────────────
ss_prevalence <- function(p, precision, alpha) {
    za <- z_alpha(alpha)
    ceiling((za^2 * p * (1 - p)) / precision^2)
}
ss_one_mean <- function(mean0, mean1, sd, alpha, power) {
    za <- z_alpha(alpha); zb <- z_beta(power)
    ceiling(((za + zb) * sd / abs(mean1 - mean0))^2)
}
ss_two_means_equal <- function(mean1, mean2, sd, alpha, power) {
    za <- z_alpha(alpha); zb <- z_beta(power)
    ceiling(((za + zb)^2 * (2 * sd^2)) / abs(mean1 - mean2)^2)
}
ss_two_means_unequal <- function(mean1, mean2, sd1, sd2, alpha, power) {
    za <- z_alpha(alpha); zb <- z_beta(power)
    ceiling(((za + zb)^2 * (sd1^2 + sd2^2)) / abs(mean1 - mean2)^2)
}
ss_paired <- function(mean1, mean2, sd_diff, alpha, power) {
    za <- z_alpha(alpha); zb <- z_beta(power)
    ceiling(((za + zb) * sd_diff / abs(mean1 - mean2))^2)
}
ss_two_proportions <- function(p1, p2, alpha, power) {
    za <- z_alpha(alpha); zb <- z_beta(power)
    ceiling(((za + zb)^2 * (p1*(1-p1) + p2*(1-p2))) / abs(p1 - p2)^2)
}
ss_logistic <- function(p1, OR, alpha, power) {
    za <- z_alpha(alpha); zb <- z_beta(power)
    ceiling((za + zb)^2 / (p1 * (1 - p1) * log(OR)^2))
}
ss_relative_risk <- function(p1, RR, alpha, power) {
    za <- z_alpha(alpha); zb <- z_beta(power)
    p2 <- p1 * RR
    ceiling(((za + zb)^2 * (p1*(1-p1) + p2*(1-p2))) / abs(p1 - p2)^2)
}
ss_poisson <- function(rate1, rate2, alpha, power) {
    za <- z_alpha(alpha); zb <- z_beta(power)
    ceiling((za + zb)^2 * (1/rate1 + 1/rate2) / log(rate1/rate2)^2)
}
ss_survival <- function(S1, S2, alpha, power) {
    za <- z_alpha(alpha); zb <- z_beta(power)
    ceiling(2 * (za + zb)^2 / log(S1/S2)^2)
}
adjust_nonresponse <- function(n, nonresponse) ceiling(n / (1 - nonresponse))

# ── Themes ────────────────────────────────────────────────────────
themes_list <- list(
    "Teal (Default)" = "teal",  "Ocean Blue" = "blue",
    "Forest Green"   = "green", "Royal Purple" = "purple",
    "Sunset Orange"  = "orange"
)
theme_colors <- list(
    teal   = list(primary = "#2C7A7B", light = "#e8f5f5", accent = "#1a5c5c"),
    blue   = list(primary = "#1a6fa8", light = "#e6f1fb", accent = "#0c447c"),
    green  = list(primary = "#3B6D11", light = "#eaf3de", accent = "#27500A"),
    purple = list(primary = "#534AB7", light = "#EEEDFE", accent = "#3C3489"),
    orange = list(primary = "#BA7517", light = "#FAEEDA", accent = "#854F0B")
)
make_bs_theme <- function(theme_name, dark_mode) {
    col <- theme_colors[[theme_name]]
    bs_theme(
        bootswatch  = if (dark_mode) "darkly" else "flatly",
        primary     = col$primary,
        "navbar-bg" = col$primary
    )
}

# ── Developer footer ──────────────────────────────────────────────
dev_footer <- function(primary, light, accent, dark_mode = FALSE) {
    bg    <- if (dark_mode) "#1e2430" else light
    text1 <- if (dark_mode) "#ffffff"  else accent
    text2 <- if (dark_mode) "#adb5bd"  else primary
    tags$div(
        style = paste0(
            "margin-top:32px; padding:18px 24px;",
            "background:", bg, ";",
            "border-top:3px solid ", primary, ";",
            "border-radius:8px;",
            "display:flex; align-items:center;",
            "justify-content:space-between; flex-wrap:wrap; gap:12px;"
        ),
        tags$div(
            style = "display:flex; align-items:center; gap:14px;",
            tags$div(
                style = paste0(
                    "width:46px; height:46px; border-radius:50%;",
                    "background:", primary, "; color:white;",
                    "display:flex; align-items:center; justify-content:center;",
                    "font-weight:bold; font-size:16px; flex-shrink:0;"
                ), "KC"
            ),
            tags$div(
                tags$p(style = paste0("margin:0; font-weight:600; font-size:1em; color:", text1, ";"),
                       "Kabinga Chanda"),
                tags$p(style = paste0("margin:2px 0 0; font-size:0.82em; color:", text2, ";"),
                       "MSc Field Epidemiology (2026) · School of Epidemiology & Biostatistics · University of Zambia"),
                tags$p(style = paste0("margin:2px 0 0; font-size:0.76em; font-style:italic; color:", text2, ";"),
                       "Ref: Lwanga & Lemeshow (WHO, 1991) · OpenEpi Standards · Prof. Musonda, UNZA")
            )
        ),
        tags$a(
            href  = "mailto:kabingachanda16@gmail.com",
            style = paste0("color:", primary, "; font-size:0.85em; text-decoration:none;"),
            icon("envelope"), " kabingachanda16@gmail.com"
        )
    )
}

# ── UI ────────────────────────────────────────────────────────────
ui <- page_sidebar(
    title    = tags$span(
        style = "font-size:1.1rem; font-weight:600; letter-spacing:0.3px;",
        icon("calculator"), " Statistical Sample Size Calculator"
    ),
    fillable = FALSE,          # ← KEY: lets cards grow naturally instead of squishing
    theme    = make_bs_theme("teal", FALSE),
    
    sidebar = sidebar(
        open  = TRUE,
        width = 270,
        padding = "16px",
        
        # ── Appearance ──
        tags$div(
            style = "margin-bottom:4px;",
            tags$p(style = "font-weight:600; font-size:0.88em; text-transform:uppercase;
                      letter-spacing:0.6px; color:#2C7A7B; margin-bottom:10px;",
                   icon("palette"), " Appearance"),
            selectInput("theme_choice", "Colour Theme",
                        choices = themes_list, selected = "teal",
                        width = "100%"),
            tags$div(
                style = "display:flex; align-items:center; justify-content:space-between;
                 padding:10px 12px; border-radius:8px;
                 border:1px solid #dee2e6; margin-top:-4px;",
                tags$span(style = "font-size:0.88em; font-weight:500;",
                          icon("moon"), "  Dark Mode"),
                tags$div(
                    class = "form-check form-switch mb-0",
                    tags$input(class = "form-check-input", type = "checkbox",
                               id = "dark_mode", role = "switch",
                               style = "cursor:pointer; width:2.4em; height:1.25em;")
                )
            )
        ),
        
        tags$hr(style = "margin:16px 0;"),
        
        # ── Study type ──
        tags$p(style = "font-weight:600; font-size:0.88em; text-transform:uppercase;
                    letter-spacing:0.6px; color:#2C7A7B; margin-bottom:10px;",
               icon("flask"), " Study Design"),
        selectInput("study_type", NULL,
                    choices = c(
                        "Prevalence Study"              = "prevalence",
                        "One Sample Mean"               = "one_mean",
                        "Two Sample Means (Equal SD)"   = "two_mean_equal",
                        "Two Sample Means (Unequal SD)" = "two_mean_unequal",
                        "Paired Samples"                = "paired",
                        "Two Sample Proportions"        = "two_prop",
                        "Logistic Regression"           = "logistic",
                        "Relative Risk"                 = "rr",
                        "Poisson Regression"            = "poisson",
                        "Survival Analysis"             = "survival"
                    ), width = "100%"
        ),
        
        uiOutput("dynamic_inputs"),
        
        tags$hr(style = "margin:16px 0;"),
        
        # ── Test parameters ──
        tags$p(style = "font-weight:600; font-size:0.88em; text-transform:uppercase;
                    letter-spacing:0.6px; color:#2C7A7B; margin-bottom:10px;",
               icon("sliders"), " Test Parameters"),
        
        conditionalPanel(
            condition = "input.study_type != 'prevalence'",
            tags$label("Power (1 - β)", class = "form-label",
                       style = "font-size:0.88em;"),
            sliderInput("power", NULL,
                        min = 0.70, max = 0.99, value = 0.80,
                        step = 0.01, ticks = FALSE, width = "100%")
        ),
        tags$label("Significance Level (α)", class = "form-label",
                   style = "font-size:0.88em;"),
        sliderInput("alpha", NULL,
                    min = 0.01, max = 0.10, value = 0.05,
                    step = 0.01, ticks = FALSE, width = "100%"),
        
        tags$hr(style = "margin:16px 0;"),
        
        # ── Non-response ──
        tags$p(style = "font-weight:600; font-size:0.88em; text-transform:uppercase;
                    letter-spacing:0.6px; color:#2C7A7B; margin-bottom:6px;",
               icon("person-x"), " Non-Response"),
        tags$label("Expected Non-Response Rate", class = "form-label",
                   style = "font-size:0.88em;"),
        sliderInput("nonresponse", NULL,
                    min = 0, max = 0.50, value = 0.10,
                    step = 0.01, ticks = FALSE, width = "100%"),
        tags$p(style = "font-size:0.78em; color:#6c757d; margin-top:-6px;",
               "Inflates n to account for dropouts / non-responders."),
        
        tags$hr(style = "margin:16px 0;"),
        
        actionButton("calculate", "  Calculate Sample Size",
                     class = "btn-primary w-100",
                     icon  = icon("calculator"),
                     style = "font-weight:600; padding:10px; font-size:0.95em;")
    ),
    
    # ── Main panel ──────────────────────────────────────────────────
    tags$div(
        style = "padding:8px 4px;",
        
        # Row 1 — three stat cards (fixed height, no text overflow)
        layout_column_wrap(
            width = 1/3,
            heights_equal = "row",
            style = "margin-bottom:20px;",
            
            value_box(
                title    = "Raw Sample Size",
                value    = textOutput("total_n"),
                showcase = bs_icon("people"),
                theme    = "primary",
                style    = "min-height:110px;"
            ),
            value_box(
                title    = "Adjusted for Non-Response",
                value    = textOutput("adj_n"),
                showcase = bs_icon("shield-check"),
                theme    = "primary",
                style    = "min-height:110px;"
            ),
            value_box(
                title    = "Study Design",
                value    = textOutput("study_label"),
                showcase = bs_icon("journal-medical"),
                theme    = "primary",
                style    = "min-height:110px; font-size:0.82rem;"
            )
        ),
        
        # Row 2 — formula | interpretation side by side
        layout_column_wrap(
            width = 1/2,
            style = "margin-bottom:20px;",
            
            card(
                height = 220,
                card_header(
                    style = "padding:10px 16px;",
                    tags$span(
                        style = "font-weight:600; font-size:0.9em;",
                        bs_icon("file-earmark-text"), "  Formula Used"
                    )
                ),
                card_body(
                    style  = "padding:14px 16px;",
                    uiOutput("formula_display")
                ),
                full_screen = TRUE
            ),
            
            card(
                height = 220,
                card_header(
                    style = "padding:10px 16px;",
                    tags$span(
                        style = "font-weight:600; font-size:0.9em;",
                        bs_icon("info-circle"), "  Interpretation"
                    )
                ),
                card_body(
                    style  = "padding:14px 16px;",
                    uiOutput("interpretation")
                ),
                full_screen = TRUE
            )
        ),
        
        # Row 3 — parameters table full width
        card(
            style = "margin-bottom:20px;",
            card_header(
                style = "padding:10px 16px;",
                tags$span(
                    style = "font-weight:600; font-size:0.9em;",
                    bs_icon("table"), "  Parameters Used"
                )
            ),
            card_body(
                style = "padding:0;",
                tableOutput("params_table")
            ),
            full_screen = TRUE
        ),
        
        # Footer
        uiOutput("footer_ui")
    )
)

# ── Server ────────────────────────────────────────────────────────
server <- function(input, output, session) {
    
    active_colors <- reactive({
        theme_colors[[ input$theme_choice %||% "teal" ]]
    })
    
    # Live theme + dark-mode switching
    observe({
        session$setCurrentTheme(
            make_bs_theme(input$theme_choice %||% "teal", isTRUE(input$dark_mode))
        )
    })
    
    output$footer_ui <- renderUI({
        col <- active_colors()
        dev_footer(col$primary, col$light, col$accent, isTRUE(input$dark_mode))
    })
    
    # Dynamic sidebar inputs
    output$dynamic_inputs <- renderUI({
        switch(input$study_type,
               "prevalence" = tagList(
                   sliderInput("p",         "Expected Prevalence (p)",         0.01, 0.99, 0.16, 0.01, ticks = F),
                   sliderInput("precision", "Precision / Margin of Error (e)", 0.01, 0.20, 0.05, 0.01, ticks = F)
               ),
               "one_mean" = tagList(
                   numericInput("mean0", "Null Mean (μ₀)",         value = 514),
                   numericInput("mean1", "Alternative Mean (μ₁)",  value = 534),
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
                   numericInput("mean1",   "Mean Before (μ₁)",        value = 259),
                   numericInput("mean2",   "Mean After (μ₂)",         value = 283),
                   numericInput("sd_diff", "SD of Differences (σ_d)", value = 129)
               ),
               "two_prop" = tagList(
                   sliderInput("p1", "Proportion Group 1 (p₁)", 0.01, 0.99, 0.19, 0.01, ticks = F),
                   sliderInput("p2", "Proportion Group 2 (p₂)", 0.01, 0.99, 0.15, 0.01, ticks = F)
               ),
               "logistic" = tagList(
                   sliderInput("p1", "Baseline Proportion (p₁)", 0.01, 0.99, 0.30, 0.01, ticks = F),
                   numericInput("OR", "Odds Ratio (OR)", value = 2.0)
               ),
               "rr" = tagList(
                   sliderInput("p1", "Baseline Proportion (p₁)", 0.01, 0.99, 0.20, 0.01, tick = F),
                   numericInput("RR", "Relative Risk (RR)", value = 1.5)
               ),
               "poisson" = tagList(
                   numericInput("rate1", "Rate Group 1 (λ₁)", value = 0.10),
                   numericInput("rate2", "Rate Group 2 (λ₂)", value = 0.20)
               ),
               "survival" = tagList(
                   sliderInput("S1", "Survival Group 1 (S₁)", 0.01, 0.99, 0.60, 0.01, ticks = F),
                   sliderInput("S2", "Survival Group 2 (S₂)", 0.01, 0.99, 0.40, 0.01, ticks = F)
               )
        )
    })
    
    # ── Main calculation ──
    results <- reactive({
        req(input$study_type, input$alpha)
        alpha <- input$alpha
        power <- if (input$study_type != "prevalence") input$power else NULL
        
        switch(input$study_type,
               
               "prevalence" = {
                   req(input$p, input$precision)
                   validate(
                       need(input$p > 0 & input$p < 1, "Prevalence must be between 0 and 1."),
                       need(input$precision > 0,       "Precision must be > 0.")
                   )
                   n <- ss_prevalence(input$p, input$precision, alpha)
                   list(total = n, per_group = NA, single_group = TRUE,
                        label   = "Prevalence Study",
                        formula = "n = (Z²_(α/2) × p × (1 - p)) / e²",
                        note    = NULL,
                        params  = data.frame(
                            Parameter = c("Prevalence (p)", "Precision (e)", "Alpha (α)", "Z_(α/2)"),
                            Value     = c(input$p, input$precision, alpha, round(z_alpha(alpha), 3))
                        ))
               },
               
               "one_mean" = {
                   req(input$mean0, input$mean1, input$sd)
                   validate(
                       need(input$sd > 0,               "SD must be > 0."),
                       need(input$mean0 != input$mean1, "Null and alternative means must differ.")
                   )
                   n <- ss_one_mean(input$mean0, input$mean1, input$sd, alpha, power)
                   list(total = n, per_group = NA, single_group = TRUE,
                        label   = "One Sample Mean",
                        formula = "n = ((Z_(α/2) + Z_β) × σ / δ)²",
                        note    = NULL,
                        params  = data.frame(
                            Parameter = c("Null Mean (μ₀)", "Alt Mean (μ₁)", "SD (σ)",
                                          "Delta (δ)", "Alpha (α)", "Power (1-β)", "Z_(α/2)", "Z_β"),
                            Value     = c(input$mean0, input$mean1, input$sd,
                                          abs(input$mean1 - input$mean0), alpha, power,
                                          round(z_alpha(alpha), 3), round(z_beta(power), 3))
                        ))
               },
               
               "two_mean_equal" = {
                   req(input$mean1, input$mean2, input$sd)
                   validate(
                       need(input$sd > 0,               "SD must be > 0."),
                       need(input$mean1 != input$mean2, "Group means must differ.")
                   )
                   n <- ss_two_means_equal(input$mean1, input$mean2, input$sd, alpha, power)
                   list(total = n * 2, per_group = n, single_group = FALSE,
                        label   = "Two Sample Means (Equal SD)",
                        formula = "n_per_group = ((Z_(α/2) + Z_β)² × 2σ²) / δ²",
                        note    = NULL,
                        params  = data.frame(
                            Parameter = c("Mean Group 1 (μ₁)", "Mean Group 2 (μ₂)", "Common SD (σ)",
                                          "Delta (δ)", "Alpha (α)", "Power (1-β)", "Z_(α/2)", "Z_β"),
                            Value     = c(input$mean1, input$mean2, input$sd,
                                          abs(input$mean1 - input$mean2), alpha, power,
                                          round(z_alpha(alpha), 3), round(z_beta(power), 3))
                        ))
               },
               
               "two_mean_unequal" = {
                   req(input$mean1, input$mean2, input$sd1, input$sd2)
                   validate(
                       need(input$sd1 > 0 & input$sd2 > 0, "Both SDs must be > 0."),
                       need(input$mean1 != input$mean2,     "Group means must differ.")
                   )
                   n <- ss_two_means_unequal(input$mean1, input$mean2, input$sd1, input$sd2, alpha, power)
                   list(total = n * 2, per_group = n, single_group = FALSE,
                        label   = "Two Sample Means (Unequal SD)",
                        formula = "n_per_group = ((Z_(α/2) + Z_β)² × (σ₁² + σ₂²)) / δ²",
                        note    = NULL,
                        params  = data.frame(
                            Parameter = c("Mean Group 1", "Mean Group 2",
                                          "SD Group 1 (σ₁)", "SD Group 2 (σ₂)",
                                          "Delta (δ)", "Alpha (α)", "Power (1-β)", "Z_(α/2)", "Z_β"),
                            Value     = c(input$mean1, input$mean2, input$sd1, input$sd2,
                                          abs(input$mean1 - input$mean2), alpha, power,
                                          round(z_alpha(alpha), 3), round(z_beta(power), 3))
                        ))
               },
               
               "paired" = {
                   req(input$mean1, input$mean2, input$sd_diff)
                   validate(
                       need(input$sd_diff > 0,          "SD of differences must be > 0."),
                       need(input$mean1 != input$mean2, "Means must differ.")
                   )
                   n <- ss_paired(input$mean1, input$mean2, input$sd_diff, alpha, power)
                   list(total = n, per_group = NA, single_group = TRUE,
                        label   = "Paired Samples",
                        formula = "n = ((Z_(α/2) + Z_β) × σ_d / δ)²",
                        note    = NULL,
                        params  = data.frame(
                            Parameter = c("Mean Before (μ₁)", "Mean After (μ₂)",
                                          "SD of Differences (σ_d)", "Delta (δ)",
                                          "Alpha (α)", "Power (1-β)", "Z_(α/2)", "Z_β"),
                            Value     = c(input$mean1, input$mean2, input$sd_diff,
                                          abs(input$mean1 - input$mean2), alpha, power,
                                          round(z_alpha(alpha), 3), round(z_beta(power), 3))
                        ))
               },
               
               "two_prop" = {
                   req(input$p1, input$p2)
                   validate(
                       need(input$p1 > 0 & input$p1 < 1, "p₁ must be between 0 and 1."),
                       need(input$p2 > 0 & input$p2 < 1, "p₂ must be between 0 and 1."),
                       need(input$p1 != input$p2,        "Proportions must differ.")
                   )
                   n <- ss_two_proportions(input$p1, input$p2, alpha, power)
                   list(total = n * 2, per_group = n, single_group = FALSE,
                        label   = "Two Sample Proportions",
                        formula = "n_per_group = ((Z_(α/2) + Z_β)² × (p₁(1-p₁) + p₂(1-p₂))) / δ²",
                        note    = NULL,
                        params  = data.frame(
                            Parameter = c("Proportion 1 (p₁)", "Proportion 2 (p₂)",
                                          "Delta (δ)", "Alpha (α)", "Power (1-β)", "Z_(α/2)", "Z_β"),
                            Value     = c(input$p1, input$p2, abs(input$p1 - input$p2),
                                          alpha, power,
                                          round(z_alpha(alpha), 3), round(z_beta(power), 3))
                        ))
               },
               
               "logistic" = {
                   req(input$p1, input$OR)
                   validate(
                       need(input$p1 > 0 & input$p1 < 1, "p₁ must be between 0 and 1."),
                       need(input$OR > 0 & input$OR != 1, "OR must be > 0 and ≠ 1.")
                   )
                   n <- ss_logistic(input$p1, input$OR, alpha, power)
                   list(total = n, per_group = NA, single_group = TRUE,
                        label   = "Logistic Regression",
                        formula = "n = (Z_(α/2) + Z_β)² / (p₁ × (1 - p₁) × ln(OR)²)",
                        note    = "⚠ Approximate — valid for single binary predictor only.",
                        params  = data.frame(
                            Parameter = c("Baseline Proportion (p₁)", "Odds Ratio (OR)",
                                          "ln(OR)", "Alpha (α)", "Power (1-β)", "Z_(α/2)", "Z_β"),
                            Value     = c(input$p1, input$OR, round(log(input$OR), 3),
                                          alpha, power,
                                          round(z_alpha(alpha), 3), round(z_beta(power), 3))
                        ))
               },
               
               "rr" = {
                   req(input$p1, input$RR)
                   p2 <- input$p1 * input$RR
                   validate(
                       need(input$p1 > 0 & input$p1 < 1, "p₁ must be between 0 and 1."),
                       need(input$RR > 0,                "RR must be > 0."),
                       need(p2 < 1,                      "p₂ = p₁ × RR exceeds 1. Reduce RR or p₁.")
                   )
                   n <- ss_relative_risk(input$p1, input$RR, alpha, power)
                   list(total = n * 2, per_group = n, single_group = FALSE,
                        label   = "Relative Risk",
                        formula = "p₂ = p₁ × RR\nn_per_group = ((Z_(α/2) + Z_β)² × (p₁(1-p₁) + p₂(1-p₂))) / δ²",
                        note    = NULL,
                        params  = data.frame(
                            Parameter = c("Baseline Proportion (p₁)", "Relative Risk (RR)",
                                          "Derived p₂", "Delta (δ)",
                                          "Alpha (α)", "Power (1-β)", "Z_(α/2)", "Z_β"),
                            Value     = c(input$p1, input$RR, round(p2, 3),
                                          round(abs(input$p1 - p2), 3), alpha, power,
                                          round(z_alpha(alpha), 3), round(z_beta(power), 3))
                        ))
               },
               
               "poisson" = {
                   req(input$rate1, input$rate2)
                   validate(
                       need(input$rate1 > 0 & input$rate2 > 0, "Both rates must be > 0."),
                       need(input$rate1 != input$rate2,         "Rates must differ.")
                   )
                   n <- ss_poisson(input$rate1, input$rate2, alpha, power)
                   list(total = n * 2, per_group = n, single_group = FALSE,
                        label   = "Poisson Regression",
                        formula = "n_per_group = ((Z_(α/2) + Z_β)² × (1/λ₁ + 1/λ₂)) / ln(λ₁/λ₂)²",
                        note    = NULL,
                        params  = data.frame(
                            Parameter = c("Rate Group 1 (λ₁)", "Rate Group 2 (λ₂)",
                                          "ln(λ₁/λ₂)", "Alpha (α)", "Power (1-β)", "Z_(α/2)", "Z_β"),
                            Value     = c(input$rate1, input$rate2,
                                          round(log(input$rate1 / input$rate2), 3),
                                          alpha, power,
                                          round(z_alpha(alpha), 3), round(z_beta(power), 3))
                        ))
               },
               
               "survival" = {
                   req(input$S1, input$S2)
                   validate(
                       need(input$S1 > 0 & input$S1 < 1, "S₁ must be between 0 and 1."),
                       need(input$S2 > 0 & input$S2 < 1, "S₂ must be between 0 and 1."),
                       need(input$S1 != input$S2,        "Survival proportions must differ.")
                   )
                   n <- ss_survival(input$S1, input$S2, alpha, power)
                   list(total = n * 2, per_group = n, single_group = FALSE,
                        label   = "Survival Analysis",
                        formula = "n_per_group = 2 × (Z_(α/2) + Z_β)² / ln(S₁/S₂)²",
                        note    = "⚠ Log-rank approx — assumes proportional hazards.",
                        params  = data.frame(
                            Parameter = c("Survival Group 1 (S₁)", "Survival Group 2 (S₂)",
                                          "ln(S₁/S₂)", "Alpha (α)", "Power (1-β)", "Z_(α/2)", "Z_β"),
                            Value     = c(input$S1, input$S2,
                                          round(log(input$S1 / input$S2), 3),
                                          alpha, power,
                                          round(z_alpha(alpha), 3), round(z_beta(power), 3))
                        ))
               }
        )
    }) |> bindEvent(input$calculate)
    
    # ── Outputs ──────────────────────────────────────────────────────
    output$total_n <- renderText({
        req(results()); format(results()$total, big.mark = ",")
    })
    
    output$adj_n <- renderText({
        req(results())
        nr  <- input$nonresponse
        raw <- results()$total
        if (nr > 0) paste0(format(adjust_nonresponse(raw, nr), big.mark = ","),
                           "  (+", round(nr * 100), "% NR)")
        else format(raw, big.mark = ",")
    })
    
    output$study_label <- renderText({
        req(results()); results()$label
    })
    
    output$formula_display <- renderUI({
        req(results())
        col  <- active_colors()
        dark <- isTRUE(input$dark_mode)
        r    <- results()
        
        bg_box  <- if (dark) "#2b2f3a" else "#f8fafb"
        fg_text <- if (dark) "#e8eaf0" else "#1c1c2e"
        ref_col <- if (dark) "#8a93a8" else "#6c757d"
        warn_col<- if (dark) "#ffc107" else "#7a5200"
        
        tagList(
            tags$pre(
                style = paste0(
                    "font-family:'Courier New',monospace; font-size:0.95em;",
                    "background:", bg_box, "; color:", fg_text, ";",
                    "padding:14px 16px;",
                    "border-left:4px solid ", col$primary, ";",
                    "border-radius:6px;",
                    "white-space:pre-wrap; word-break:break-word;",
                    "line-height:1.8; margin:0;"
                ),
                r$formula
            ),
            if (!is.null(r$note))
                tags$p(style = paste0("font-size:0.82em; color:", warn_col, ";
                               margin:8px 0 4px; line-height:1.5;"), r$note),
            tags$p(style = paste0("font-size:0.78em; color:", ref_col, ";
                             margin:6px 0 0; line-height:1.5;"),
                   "Ref: Lwanga & Lemeshow (WHO, 1991) · OpenEpi · Prof. Musonda, UNZA")
        )
    })
    
    output$interpretation <- renderUI({
        req(results())
        r    <- results()
        nr   <- input$nonresponse
        dark <- isTRUE(input$dark_mode)
        raw  <- r$total
        adj  <- if (nr > 0) adjust_nonresponse(raw, nr) else raw
        
        base <- if (r$single_group) {
            paste0("A total of ", tags$strong(format(raw, big.mark = ",")),
                   " participants are required for this ",
                   tags$em(r$label), " study (single group).")
        } else {
            paste0(tags$strong(format(r$per_group, big.mark = ",")),
                   " participants per group — ",
                   tags$strong(format(raw, big.mark = ",")),
                   " total across both groups.")
        }
        
        adj_note <- if (nr > 0) {
            tagList(
                tags$hr(style = "margin:10px 0;"),
                tags$p(style = "margin:0; font-size:0.9em;",
                       icon("triangle-exclamation"), " After ",
                       tags$strong(paste0(round(nr * 100), "%")),
                       " non-response adjustment: ",
                       tags$strong(format(adj, big.mark = ",")), " participants.")
            )
        }
        
        tags$div(
            style = "font-size:0.95em; line-height:1.8;",
            tags$p(style = "margin:0;", HTML(base)),
            adj_note
        )
    })
    
    output$params_table <- renderTable({
        req(results())
        nr  <- input$nonresponse
        raw <- results()$total
        df  <- results()$params
        if (nr > 0) {
            df <- rbind(df, data.frame(
                Parameter = c("Non-Response Rate", "Adjusted Total n"),
                Value     = c(nr, adjust_nonresponse(raw, nr))
            ))
        }
        df
    }, striped = TRUE, hover = TRUE, bordered = TRUE,
    width = "100%", align = "lr")
}

shinyApp(ui = ui, server = server)