library(shiny)
library(bslib)

APP_VERSION <- "v2.0 \u00b7 June 2026"

# ── Contaminant lookup table (name, RfD mg/kg/day) ───────────────────────────
contaminants <- data.frame(
  name = c("Custom", "Lead", "Arsenic", "Nitrates", "Fluoride",
           "Cadmium", "Mercury", "Chromium (VI)", "Benzene", "Chloroform"),
  rfd  = c(NA, 0.0035, 0.0003, 1.6, 0.06,
           0.001, 0.0003, 0.003, 0.004, 0.01),
  stringsAsFactors = FALSE
)

# ── ZDHS defaults ─────────────────────────────────────────────────────────────
zdhs_defaults <- list(
  C=0.05, IR=2, EF=350, ED=30, BW=65, AT=10950,
  ET=8, SA=5700, AF=0.07, ABS=0.1
)

# ── Validation limits ─────────────────────────────────────────────────────────
param_limits <- list(
  C  =list(min=0,   max=1e6,  warn="Must be > 0"),
  IR =list(min=0,   max=100,  warn="Typical: 0.5\u20135 L/day"),
  EF =list(min=1,   max=365,  warn="Must be 1\u2013365 days/year"),
  ED =list(min=0,   max=70,   warn="Must be 0\u201370 years"),
  BW =list(min=1,   max=300,  warn="Typical: 40\u2013120 kg"),
  AT =list(min=1,   max=25550,warn="Non-cancer: ED\u00d7365; Cancer: 25,550"),
  ET =list(min=0,   max=24,   warn="Must be 0\u201324 hr/day"),
  SA =list(min=0,   max=20000,warn="Typical adult: 3,000\u20138,000 cm\u00b2"),
  AF =list(min=0,   max=10,   warn="Typical: 0.01\u20131 mg/cm\u00b2"),
  ABS=list(min=0,   max=1,    warn="Must be 0\u20131")
)

# ── Mode parameter definitions ────────────────────────────────────────────────
mode_params <- list(
  Ingestion = list(
    list(id="C",  label="Contaminant concentration", unit="mg/L",       def=0.05),
    list(id="IR", label="Ingestion rate",             unit="L/day",      def=2),
    list(id="EF", label="Exposure frequency",         unit="days/year",  def=350),
    list(id="ED", label="Exposure duration",          unit="years",      def=30),
    list(id="BW", label="Body weight",                unit="kg",         def=65),
    list(id="AT", label="Averaging time",             unit="days",       def=10950)
  ),
  Inhalation = list(
    list(id="C",  label="Contaminant concentration", unit="mg/m\u00b3", def=0.01),
    list(id="IR", label="Inhalation rate",            unit="m\u00b3/hr", def=0.83),
    list(id="ET", label="Exposure time",              unit="hr/day",     def=8),
    list(id="EF", label="Exposure frequency",         unit="days/year",  def=250),
    list(id="ED", label="Exposure duration",          unit="years",      def=25),
    list(id="BW", label="Body weight",                unit="kg",         def=65),
    list(id="AT", label="Averaging time",             unit="days",       def=9125)
  ),
  `Skin Absorption` = list(
    list(id="C",   label="Contaminant concentration", unit="mg/cm\u00b2",def=0.001),
    list(id="SA",  label="Skin surface area",          unit="cm\u00b2",  def=5700),
    list(id="AF",  label="Skin adherence factor",      unit="mg/cm\u00b2",def=0.07),
    list(id="ABS", label="Absorption factor",          unit="0\u20131",  def=0.1),
    list(id="EF",  label="Exposure frequency",         unit="days/year", def=350),
    list(id="ED",  label="Exposure duration",          unit="years",     def=30),
    list(id="BW",  label="Body weight",                unit="kg",        def=65),
    list(id="AT",  label="Averaging time",             unit="days",      def=10950)
  )
)

# ── Helpers ───────────────────────────────────────────────────────────────────
tip <- function(id) {
  lim <- param_limits[[id]]; dv <- zdhs_defaults[[id]]
  tags$span(style="cursor:help;color:var(--bs-primary);margin-left:4px;font-size:0.75rem;",
            title=paste0("ZDHS default: ", dv, " \u00b7 ", lim$warn), "\u24d8")
}

param_inputs <- function(mode) {
  params <- mode_params[[mode]]
  inputs <- lapply(params, function(p)
    tags$div(class="param-cell",
             numericInput(paste0("param_", tolower(p$id)),
                          tagList(p$label, tags$small(class="text-muted ms-1", p$unit), tip(p$id)),
                          NA, min=0, step=NA, width="100%")))
  tags$div(class="param-grid", inputs)
}

fmt_num <- function(x, digits=7)
  format(round(x, digits), scientific=FALSE, big.mark=",")

# ── CDI calculation function ──────────────────────────────────────────────────
calc_cdi <- function(mode, v) {
  denom <- v[["BW"]] * v[["AT"]]
  if (is.null(denom) || denom == 0) return(NA)
  switch(mode,
         Ingestion        = (v$C * v$IR * v$EF * v$ED) / denom,
         Inhalation       = (v$C * v$IR * v$ET * v$EF * v$ED) / denom,
         `Skin Absorption`= (v$C * v$SA * v$AF * v$ABS * v$EF * v$ED) / denom,
         NA)
}

# ── Themes ────────────────────────────────────────────────────────────────────
themes <- list(
  "WHO Geneva"     = bs_theme(bg="#FFFFFF",fg="#1A1A2E",primary="#005A9E",secondary="#0085CA",success="#2E8B57",warning="#E07B00",danger="#CC2936",base_font=font_google("Noto Sans"),heading_font=font_google("Noto Sans"),code_font=font_google("Fira Mono")),
  "World Bank"     = bs_theme(bg="#F5F5F0",fg="#1C1C1C",primary="#003087",secondary="#B8860B",success="#2E7D32",warning="#F57C00",danger="#C62828",base_font=font_google("Open Sans"),heading_font=font_google("Open Sans"),code_font=font_google("Fira Mono")),
  "US Federal"     = bs_theme(bg="#F0F2F5",fg="#1B2A3B",primary="#005B6E",secondary="#4A8FA0",success="#2E6B3E",warning="#C8730A",danger="#A4262C",base_font=font_google("Source Sans 3"),heading_font=font_google("Source Sans 3"),code_font=font_google("Fira Mono")),
  "Corporate Dark" = bs_theme(bg="#1C2230",fg="#E8EBF0",primary="#4A9FD4",secondary="#7EC8B0",success="#56A063",warning="#E0A030",danger="#D95F5F",base_font=font_google("Inter"),heading_font=font_google("Inter"),code_font=font_google("Fira Mono")),
  "Academic"       = bs_theme(bg="#FAFAF7",fg="#2C2416",primary="#7B2D3E",secondary="#A0522D",success="#3A6B40",warning="#B8860B",danger="#8B0000",base_font=font_google("Lora"),heading_font=font_google("Lora"),code_font=font_google("Fira Mono")),
  "Microsoft"      = bs_theme(bg="#F3F2F1",fg="#201F1E",primary="#0078D4",secondary="#2B88D8",success="#107C10",warning="#FF8C00",danger="#D13438",base_font=font_google("Inter"),heading_font=font_google("Inter"),code_font=font_google("Fira Mono")),
  "IBM Carbon"     = bs_theme(bg="#F4F4F4",fg="#161616",primary="#0F62FE",secondary="#0043CE",success="#198038",warning="#F1C21B",danger="#DA1E28",base_font=font_google("IBM Plex Sans"),heading_font=font_google("IBM Plex Sans"),code_font=font_google("IBM Plex Mono")),
  "Google Material"= bs_theme(bg="#FFFFFF",fg="#202124",primary="#1A73E8",secondary="#34A853",success="#34A853",warning="#FBBC04",danger="#EA4335",base_font=font_google("Roboto"),heading_font=font_google("Roboto"),code_font=font_google("Roboto Mono")),
  "Midnight Pro"   = bs_theme(bg="#0D1117",fg="#C9D1D9",primary="#58A6FF",secondary="#3FB950",success="#3FB950",warning="#D29922",danger="#F85149",base_font=font_google("Inter"),heading_font=font_google("Inter"),code_font=font_google("Fira Mono")),
  "Teal Executive" = bs_theme(bg="#F7F9F9",fg="#0D2B2B",primary="#007A7A",secondary="#4DB6AC",success="#388E3C",warning="#F57C00",danger="#C62828",base_font=font_google("DM Sans"),heading_font=font_google("DM Sans"),code_font=font_google("Fira Mono"))
)

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- page_fillable(
  title = "Environmental Human Risk Assessment",
  theme = themes[["WHO Geneva"]],
  padding = "0.75rem",
  
  tags$head(tags$style(HTML("
    /* Title bar */
    .title-bar { display:flex; align-items:center; justify-content:space-between;
      padding:0.7rem 1.4rem; background-color:var(--bs-primary);
      margin-bottom:0.75rem; flex-shrink:0; }
    .title-bar .app-title    { font-size:1.05rem; font-weight:700; color:#fff; margin:0; }
    .title-bar .app-subtitle { font-size:0.72rem; color:rgba(255,255,255,0.78); margin:0.1rem 0 0; }
    .title-bar .theme-wrap   { display:flex; align-items:center; gap:0.6rem; }
    .title-bar label         { font-size:0.75rem; color:rgba(255,255,255,0.85); margin:0; white-space:nowrap; }
    .title-bar select        { font-size:0.75rem; padding:0.2rem 0.5rem; min-width:140px; }
    .title-bar .shiny-input-container { margin:0 !important; }
    /* Params */
    /* param grid and header inputs handled below */
    /* Results */
    .result-box { background:var(--bs-tertiary-bg); border-radius:0.4rem; padding:0.85rem 1rem; }
    .formula-code { font-family:monospace; font-size:0.8rem; opacity:0.8; }
    .risk-badge { font-size:0.8rem; padding:0.3rem 0.7rem; border-radius:0.4rem; }
    .interp-box { border-left:3px solid var(--bs-primary); background:var(--bs-tertiary-bg);
      border-radius:0 0.4rem 0.4rem 0; padding:0.7rem 1rem; font-size:0.82rem; line-height:1.65; }
    .progress-bar { transition:width 0.6s ease; }
    /* Sensitivity chart */
    .sens-bar-wrap { display:flex; align-items:center; gap:0.5rem; margin-bottom:0.35rem; }
    .sens-label { font-size:0.75rem; min-width:40px; text-align:right; }
    .sens-track { flex:1; height:14px; background:var(--bs-tertiary-bg); border-radius:3px; overflow:hidden; position:relative; }
    .sens-pos { height:100%; background:var(--bs-danger); opacity:0.7; border-radius:3px; }
    .sens-pct { font-size:0.72rem; min-width:38px; }
    /* Tabs */
    .tab-nav { display:flex; gap:0; margin-bottom:1rem; border-bottom:1px solid var(--bs-border-color); }
    .tab-btn { background:none; border:none; border-bottom:2px solid transparent;
      padding:0.4rem 0.9rem; font-size:0.8rem; cursor:pointer;
      color:var(--bs-secondary-color); margin-bottom:-1px;
      display:flex; align-items:center; gap:0.3rem; }
    .tab-btn:hover { color:var(--bs-primary); }
    .tab-btn.active { border-bottom-color:var(--bs-primary); color:var(--bs-primary); font-weight:600; }
    .tab-pane { display:none; } .tab-pane.active { display:block; }
    /* Reference */
    .ref-card { border-left:3px solid var(--bs-primary); background:var(--bs-secondary-bg);
      border-radius:0 0.4rem 0.4rem 0; padding:0.8rem 1rem; margin-bottom:0.75rem; }
    .ref-param { font-weight:600; font-size:0.85rem; margin-bottom:0.25rem; }
    .ref-desc  { font-size:0.8rem; margin-bottom:0.4rem; line-height:1.6; }
    .ref-source { font-size:0.78rem; margin-bottom:0.2rem; display:flex; gap:0.4rem; }
    .ref-badge { font-size:0.65rem; font-weight:700; padding:0.1rem 0.45rem;
      border-radius:0.25rem; white-space:nowrap; flex-shrink:0; margin-top:0.1rem; }
    .badge-zdhs { background:#0077b6; color:#fff; }
    .badge-codex{ background:#2d6a4f; color:#fff; }
    .cdi-box { background:var(--bs-primary); color:#fff;
      border-radius:0.5rem; padding:1rem 1.1rem; margin-bottom:1rem; }
    .cdi-box h6 { color:#fff; margin-bottom:0.4rem; font-size:0.9rem; }
    .cdi-box p  { font-size:0.8rem; opacity:0.9; margin:0; line-height:1.65; }
    .glossary-table th,.glossary-table td { font-size:0.78rem; vertical-align:top; padding:0.3rem 0.5rem; }
    .example-box { background:var(--bs-secondary-bg); border-radius:0.4rem; padding:0.9rem 1rem; margin-bottom:1rem; }
    .example-box h6 { font-size:0.85rem; font-weight:600; margin-bottom:0.5rem; }
    .example-row { display:flex; justify-content:space-between; font-size:0.8rem;
      padding:0.2rem 0; border-bottom:1px solid var(--bs-border-color); }
    .example-row:last-child { border-bottom:none; }
    /* History table */
    .history-table th,.history-table td { font-size:0.78rem; vertical-align:middle; padding:0.3rem 0.5rem; }
    /* About */
    .dev-avatar { width:60px; height:60px; border-radius:50%; background:var(--bs-primary);
      color:#fff; font-size:1.3rem; font-weight:700;
      display:flex; align-items:center; justify-content:center; flex-shrink:0; }
    .dev-row  { display:flex; align-items:center; gap:1rem; margin-bottom:1.1rem; }
    .dev-name { font-size:1rem; font-weight:700; margin:0; }
    .dev-degree { font-size:0.78rem; opacity:0.65; margin:0.1rem 0 0; }
    .dev-detail { display:flex; align-items:center; gap:0.5rem; font-size:0.82rem; margin-bottom:0.45rem; }
    .version-badge { display:inline-block; font-size:0.7rem; padding:0.2rem 0.6rem;
      background:var(--bs-secondary-bg); border:1px solid var(--bs-border-color);
      border-radius:0.3rem; margin-bottom:1rem; }
    .disclaimer { border-left:3px solid var(--bs-warning); background:var(--bs-secondary-bg);
      border-radius:0 0.4rem 0.4rem 0; padding:0.75rem 1rem;
      font-size:0.8rem; line-height:1.65; margin-top:0.75rem; }
    /* Misc */
    .btn-row { display:flex; gap:0.5rem; margin-top:0.5rem; }
    .btn-row .btn { flex:1; }
    /* Two-column flexbox layout */
    .app-columns { display:flex; gap:1rem; flex:1 1 0; min-height:0; overflow:hidden; }
    .left-panel  { flex:0 0 48%; min-width:0; overflow-y:auto;
      background:var(--bs-card-bg); border:1px solid var(--bs-border-color);
      border-radius:0.5rem; padding:0.75rem; box-sizing:border-box;
      align-self:flex-start; }
    .right-panel { flex:1 1 0; min-width:0; overflow:hidden;
      display:flex; flex-direction:column; }
    .right-panel .card { flex:1 1 0; min-height:0; }
    .lp-section  { margin-bottom:0.6rem; }
    .lp-label    { font-size:0.78rem; font-weight:600; margin-bottom:0.2rem; opacity:0.75; }
    .lp-grid2    { display:grid; grid-template-columns:1fr 1fr; gap:0 0.75rem; }
    .lp-section .shiny-input-container { margin-bottom:0 !important; }
    /* ── Unified typography & sizing system ── */
    /* All labels: same size, weight, spacing */
    .left-panel label,
    .left-panel .lp-label,
    .left-panel .control-label { font-size:0.82rem !important; font-weight:500 !important;
      margin-bottom:3px !important; line-height:1.3 !important; color:inherit !important; }
    /* All text inputs and selects: same height and font */
    .left-panel input.form-control,
    .left-panel select.form-select { font-size:0.82rem !important;
      padding:0.28rem 0.5rem !important;
      height:calc(1.5 * 0.82rem + 0.56rem + 2px) !important; }
    /* Remove extra Shiny wrapper margins */
    .left-panel .form-group,
    .left-panel .shiny-input-container { margin-bottom:0.55rem !important; width:100% !important; }
    /* Radio buttons consistent */
    .left-panel .radio-inline,
    .left-panel .shiny-options-group { font-size:0.82rem !important; }
    .left-panel .radio-inline + .radio-inline { margin-left:0.75rem !important; }
    /* Small unit tags next to labels */
    .left-panel label small,
    .left-panel .control-label small { font-size:0.72rem !important; font-weight:400 !important;
      opacity:0.65; margin-left:3px !important; }
    /* Tooltip icon */
    .left-panel label span[title] { font-size:0.72rem !important; }
    /* Param grid */
    .param-grid { display:grid; grid-template-columns:1fr 1fr; gap:0 0.75rem;
      margin-bottom:0.4rem; }
    .param-cell .form-group,
    .param-cell .shiny-input-container { width:100% !important; margin-bottom:0.5rem !important; }
    .fill-progress-track { height:3px; background:var(--bs-border-color);
      border-radius:2px; margin-bottom:0.25rem; overflow:hidden; }
    .fill-progress { height:100%; background:var(--bs-primary); border-radius:2px;
      transition:width 0.3s ease; }
    .fill-label { font-size:0.72rem; color:var(--bs-secondary-color); margin-bottom:0.4rem; }
    /* Print */
    @media print {
      .title-bar,.btn-row,.tab-nav,.card-header { -webkit-print-color-adjust:exact; }
      .tab-pane { display:block !important; page-break-before:always; }
      #pane-reference,#pane-about { display:none !important; }
    }
  "))),
  
  # ── Title bar ────────────────────────────────────────────────────────────────
  div(class="title-bar",
      div(
        p(class="app-title",    "Environmental Human Risk Assessment"),
        p(class="app-subtitle", "Chronic Daily Intake (CDI)")
      ),
      div(class="theme-wrap",
          input_dark_mode(id="dark_mode", mode="light"),
          tags$label(`for`="theme_select","Theme:"),
          tags$select(id="theme_select", class="form-select form-select-sm",
                      onchange="Shiny.setInputValue('theme_choice',this.value)",
                      lapply(names(themes), function(nm)
                        tags$option(value=nm, selected=if(nm=="WHO Geneva")"selected" else NULL, nm))
          )
      )
  ),
  
  # Two-column layout using plain flexbox — avoids bslib grid height issues
  div(class="app-columns",
      
      # ── LEFT ──────────────────────────────────────────────────────────────────
      div(class="left-panel",
          
          # Section: Exposure route
          div(class="lp-section",
              div(class="lp-label", bsicons::bs_icon("signpost-split"), " Exposure Route"),
              radioButtons("mode", NULL, choices=names(mode_params),
                           selected="Ingestion", inline=TRUE)
          ),
          
          # Section: Contaminant + RfD
          div(class="lp-section lp-grid2",
              div(
                div(class="lp-label", bsicons::bs_icon("search"), " Contaminant"),
                selectInput("contaminant", NULL,
                            choices=contaminants$name, selected="Custom", width="100%")
              ),
              div(
                div(class="lp-label", bsicons::bs_icon("shield-check"),
                    " RfD (mg/kg\u00b7day)"),
                numericInput("rfd", NULL, value=1.0, min=0, step=NA, width="100%")
              )
          ),
          
          # Fill progress
          div(style="margin:0.3rem 0 0.5rem;",
              div(class="fill-progress-track",
                  div(class="fill-progress", id="fill-bar", style="width:0%;")
              ),
              uiOutput("fill_status")
          ),
          
          # Parameter grid — renders naturally, no height tricks
          uiOutput("param_ui"),
          uiOutput("validation_ui"),
          
          # Buttons
          div(class="btn-row lp-section",
              actionButton("prefill",
                           tagList(bsicons::bs_icon("lightning-fill"), " ZDHS Defaults"),
                           class="btn-outline-secondary btn-sm"),
              actionButton("reset_btn",
                           tagList(bsicons::bs_icon("arrow-counterclockwise"), " Reset"),
                           class="btn-outline-secondary btn-sm"),
              actionButton("calculate",
                           tagList(bsicons::bs_icon("calculator-fill"), " Calculate"),
                           class="btn-primary btn-sm")
          )
      ),
      
      # ── RIGHT ─────────────────────────────────────────────────────────────────
      div(class="right-panel",
          card(fill=TRUE, full_screen=TRUE, style="height:100%;margin:0;",
               card_header(
                 div(class="tab-nav",
                     tags$button(id="tab-results",   class="tab-btn active",  onclick="switchTab('results')",
                                 bsicons::bs_icon("bar-chart-line")," Results"),
                     tags$button(id="tab-analysis",  class="tab-btn",         onclick="switchTab('analysis')",
                                 bsicons::bs_icon("graph-up")," Analysis"),
                     tags$button(id="tab-history",   class="tab-btn",         onclick="switchTab('history')",
                                 bsicons::bs_icon("clock-history")," History"),
                     tags$button(id="tab-reference", class="tab-btn",         onclick="switchTab('reference')",
                                 bsicons::bs_icon("journals")," Reference"),
                     tags$button(id="tab-about",     class="tab-btn",         onclick="switchTab('about')",
                                 bsicons::bs_icon("person-badge")," About")
                 )
               ),
               card_body(fillable=TRUE, style="overflow-y:auto;",
                         
                         # Results pane
                         div(id="pane-results", class="tab-pane active",
                             uiOutput("output_ui")
                         ),
                         
                         # Analysis pane (sensitivity + exposure pathway diagram)
                         div(id="pane-analysis", class="tab-pane",
                             uiOutput("analysis_ui")
                         ),
                         
                         # History pane
                         div(id="pane-history", class="tab-pane",
                             uiOutput("history_ui")
                         ),
                         
                         # Reference pane
                         div(id="pane-reference", class="tab-pane",
                             div(class="cdi-box",
                                 tags$h6(bsicons::bs_icon("info-circle")," What is CDI?"),
                                 tags$p("The Chronic Daily Intake (CDI) estimates the average daily dose of a
              chemical over a long-term exposure period (30 years non-cancer; 70 years
              carcinogens). Expressed in mg\u00b7kg\u207b\u00b9\u00b7day\u207b\u00b9.
              Hazard Quotient HQ = CDI \u00f7 RfD; HQ < 1 = acceptable non-carcinogenic
              risk (US EPA RAGS, 1989).")
                             ),
                             tags$h6(class="fw-bold mb-2", bsicons::bs_icon("card-list")," Glossary"),
                             tags$table(class="table table-sm table-bordered glossary-table mb-3",
                                        tags$thead(tags$tr(tags$th("Term"),tags$th("Full Name"),tags$th("Definition"))),
                                        tags$tbody(
                                          tags$tr(tags$td("CDI"), tags$td("Chronic Daily Intake"),   tags$td("Average daily contaminant dose over an exposure period")),
                                          tags$tr(tags$td("HQ"),  tags$td("Hazard Quotient"),        tags$td("CDI \u00f7 RfD; HQ < 1 = acceptable non-cancer risk")),
                                          tags$tr(tags$td("RfD"), tags$td("Reference Dose"),         tags$td("Acceptable daily intake; contaminant-specific (EPA IRIS)")),
                                          tags$tr(tags$td("RAGS"),tags$td("Risk Assessment Guidance"),tags$td("US EPA Superfund framework for human health risk assessment")),
                                          tags$tr(tags$td("ZDHS"),tags$td("Zambia DHS"),             tags$td("National survey providing Zambian population exposure defaults")),
                                          tags$tr(tags$td("MRL"), tags$td("Maximum Residue Limit"),  tags$td("Codex Alimentarius safe upper contaminant concentration"))
                                        )
                             ),
                             tags$h6(class="fw-bold mb-2", bsicons::bs_icon("table")," Parameter Reference"),
                             div(class="ref-card",
                                 p(class="ref-param","C \u2014 Contaminant Concentration"),
                                 p(class="ref-desc","Amount in exposure medium from environmental sampling. Units: mg/L (ingestion), mg/m\u00b3 (inhalation), mg/cm\u00b2 (dermal)."),
                                 div(class="ref-source",span(class="ref-badge badge-zdhs","ZDHS"),span("From field environmental sampling per ZDHS protocols.")),
                                 div(class="ref-source",span(class="ref-badge badge-codex","Codex"),span("MRLs per Codex Stan 193-1995."))),
                             div(class="ref-card",
                                 p(class="ref-param","IR \u2014 Ingestion / Inhalation Rate"),
                                 p(class="ref-desc","Volume of medium contacted per unit time. Reflects individual physiology and behaviour."),
                                 div(class="ref-source",span(class="ref-badge badge-zdhs","ZDHS"),span("ZDHS 2018: water \u223c2 L/day; air 0.83 m\u00b3/hr.")),
                                 div(class="ref-source",span(class="ref-badge badge-codex","Codex"),span("Codex CAC/RCP 1-1969: national per-capita consumption data."))),
                             div(class="ref-card",
                                 p(class="ref-param","EF \u2014 Exposure Frequency"),
                                 p(class="ref-desc","Days/year of contact. Values below 365 account for time away from source."),
                                 div(class="ref-source",span(class="ref-badge badge-zdhs","ZDHS"),span("350 days/year.")),
                                 div(class="ref-source",span(class="ref-badge badge-codex","Codex"),span("365 days/year for chronic dietary assessments."))),
                             div(class="ref-card",
                                 p(class="ref-param","ED \u2014 Exposure Duration"),
                                 p(class="ref-desc","Total years of exposure. Standard non-cancer = 30 years; cancer = 70 years."),
                                 div(class="ref-source",span(class="ref-badge badge-zdhs","ZDHS"),span("30 years adult residential.")),
                                 div(class="ref-source",span(class="ref-badge badge-codex","Codex"),span("70 years for lifetime carcinogenic risk."))),
                             div(class="ref-card",
                                 p(class="ref-param","BW \u2014 Body Weight"),
                                 p(class="ref-desc","Normalises dose per kg. Lower BW = higher CDI. Use population-specific values."),
                                 div(class="ref-source",span(class="ref-badge badge-zdhs","ZDHS"),span("Male 68 kg, female 62 kg (ZDHS 2018 adult means).")),
                                 div(class="ref-source",span(class="ref-badge badge-codex","Codex"),span("60 kg; US EPA: 70 kg."))),
                             div(class="ref-card",
                                 p(class="ref-param","AT \u2014 Averaging Time"),
                                 p(class="ref-desc","Period over which dose is averaged. Non-cancer: ED\u00d7365. Cancer: 25,550 days."),
                                 div(class="ref-source",span(class="ref-badge badge-zdhs","ZDHS"),span("Non-cancer: ED\u00d7365. Cancer: 25,550 days.")),
                                 div(class="ref-source",span(class="ref-badge badge-codex","Codex"),span("Full lifetime for chronic risk."))),
                             tags$h6(class="fw-bold mb-2 mt-3", bsicons::bs_icon("play-circle")," Worked Example \u2014 Zambia Water Ingestion"),
                             div(class="example-box",
                                 tags$h6("Scenario: Lead-contaminated borehole water, Lusaka community"),
                                 div(class="example-row",span("C \u2014 Lead"),       span("0.05 mg/L")),
                                 div(class="example-row",span("IR"),                  span("2 L/day")),
                                 div(class="example-row",span("EF"),                  span("350 days/year")),
                                 div(class="example-row",span("ED"),                  span("30 years")),
                                 div(class="example-row",span("BW"),                  span("65 kg (ZDHS)")),
                                 div(class="example-row",span("AT"),                  span("10,950 days")),
                                 tags$hr(class="my-2"),
                                 div(class="example-row fw-bold",
                                     span("CDI = (0.05\u00d72\u00d7350\u00d730)/(65\u00d710,950)"),
                                     span("= 0.000493 mg/kg\u00b7day")),
                                 div(class="example-row fw-bold",
                                     span("HQ (RfD lead = 0.0035)"), span("= 0.141 \u2192 Acceptable"))
                             ),
                             tags$p(class="text-muted",style="font-size:0.72rem;",
                                    "Sources: ZDHS 2018; Codex Stan 193-1995; US EPA RAGS Vol. I (1989); WHO EHC.")
                         ),
                         
                         # About pane
                         div(id="pane-about", class="tab-pane",
                             span(class="version-badge", APP_VERSION),
                             div(class="dev-row",
                                 div(class="dev-avatar","KC"),
                                 div(p(class="dev-name","Kabinga Chanda"),
                                     p(class="dev-degree","MSc Field Epidemiology \u00b7 School of Public Health"),
                                     p(class="dev-degree","University of Zambia \u00b7 Class of 2026"))
                             ),
                             div(class="result-box mb-3",
                                 div(class="dev-detail",bsicons::bs_icon("envelope-fill"),
                                     tags$a(href="mailto:kabingachanda16@gmail.com","kabingachanda16@gmail.com")),
                                 div(class="dev-detail",bsicons::bs_icon("mortarboard-fill"),span("MSc Field Epidemiology")),
                                 div(class="dev-detail",bsicons::bs_icon("hospital"),span("School of Public Health, UNZA")),
                                 div(class="dev-detail",bsicons::bs_icon("geo-alt-fill"),span("Lusaka, Zambia"))
                             ),
                             div(class="disclaimer",
                                 bsicons::bs_icon("exclamation-triangle-fill")," ",
                                 tags$strong("Academic use only.")," ",
                                 "Developed as an MSc Field Epidemiology academic project, UNZA (2026).
             For educational purposes only. Not for regulatory or clinical decisions."
                             )
                         )
               ),
               
               tags$script(HTML("
        function switchTab(name){
          ['results','analysis','history','reference','about'].forEach(function(t){
            document.getElementById('pane-'+t).classList.remove('active');
            document.getElementById('tab-'+t).classList.remove('active');
          });
          document.getElementById('pane-'+name).classList.add('active');
          document.getElementById('tab-'+name).classList.add('active');
        }
        // Enter key triggers calculate
        document.addEventListener('keydown', function(e){
          if(e.key==='Enter' && document.activeElement.tagName==='INPUT')
            document.getElementById('calculate').click();
        });
      "))
          )  # end card
      )  # end right-panel
  )  # end app-columns
)

# ── Server ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {
  
  # Theme
  observeEvent(input$theme_choice, {
    req(input$theme_choice %in% names(themes))
    session$setCurrentTheme(themes[[input$theme_choice]])
  })
  
  # Contaminant lookup → update RfD
  observeEvent(input$contaminant, {
    row <- contaminants[contaminants$name == input$contaminant, ]
    if (nrow(row) > 0 && !is.na(row$rfd))
      updateNumericInput(session, "rfd", value = row$rfd)
  })
  
  # Param inputs
  output$param_ui <- renderUI({
    mode <- if (is.null(input$mode)) "Ingestion" else input$mode
    param_inputs(mode)
  })
  
  # ZDHS prefill
  observeEvent(input$prefill, {
    mode <- if (is.null(input$mode)) "Ingestion" else input$mode
    lapply(mode_params[[mode]], function(p)
      updateNumericInput(session, paste0("param_", tolower(p$id)), value=p$def))
  })
  
  # Reset
  observeEvent(input$reset_btn, {
    mode <- if (is.null(input$mode)) "Ingestion" else input$mode
    lapply(mode_params[[mode]], function(p)
      updateNumericInput(session, paste0("param_", tolower(p$id)), value=NA))
  })
  
  # Collect values
  param_values <- reactive({
    mode   <- if (is.null(input$mode)) "Ingestion" else input$mode
    params <- mode_params[[mode]]
    vals   <- lapply(params, function(p) input[[paste0("param_", tolower(p$id))]])
    names(vals) <- sapply(params, `[[`, "id")
    vals
  })
  
  # Fill progress
  output$fill_status <- renderUI({
    v     <- param_values()
    total <- length(v)
    done  <- sum(!sapply(v, is.null) & !sapply(v, is.na))
    pct   <- round(done / total * 100)
    tagList(
      tags$script(paste0("document.getElementById('fill-bar').style.width='",pct,"%';")),
      tags$div(class="fill-label", done, " / ", total, " fields complete")
    )
  })
  
  # Inline validation
  output$validation_ui <- renderUI({
    v     <- param_values()
    filled <- Filter(Negate(is.na), v)
    msgs  <- lapply(names(filled), function(id) {
      val <- as.numeric(filled[[id]]); lim <- param_limits[[id]]
      if (!is.null(lim) && !is.na(val) && (val < lim$min || val > lim$max))
        div(class="text-warning", style="font-size:0.75rem;margin-bottom:2px;",
            bsicons::bs_icon("exclamation-circle")," ", id, ": ", lim$warn)
    })
    msgs <- Filter(Negate(is.null), msgs)
    if (length(msgs)==0) return(NULL)
    div(class="mt-1", msgs)
  })
  
  # CDI reactive
  cdi_result <- eventReactive(input$calculate, {
    mode <- if (is.null(input$mode)) "Ingestion" else input$mode
    v    <- param_values()
    if (any(sapply(v, is.na)))
      return(list(error="Please fill in all parameter fields."))
    v <- lapply(v, as.numeric)
    if ((v[["BW"]] * v[["AT"]]) == 0)
      return(list(error="BW and AT must be > 0."))
    warn_ids <- names(Filter(function(id){
      val <- v[[id]]; lim <- param_limits[[id]]
      !is.null(lim) && !is.na(val) && (val < lim$min || val > lim$max)
    }, names(v)))
    cdi <- calc_cdi(mode, v)
    list(cdi=cdi, hq=cdi/(input$rfd %||% 1), rfd=input$rfd %||% 1,
         mode=mode, v=v, warnings=warn_ids, error=NULL,
         contaminant=input$contaminant, ts=format(Sys.time(),"%H:%M:%S"))
  })
  
  # Calculation history (up to 10)
  history <- reactiveVal(list())
  observeEvent(input$calculate, {
    res <- cdi_result()
    if (!is.null(res$error)) return()
    h <- history()
    entry <- data.frame(
      Time        = res$ts,
      Mode        = res$mode,
      Contaminant = res$contaminant,
      CDI         = fmt_num(res$cdi, 6),
      HQ          = fmt_num(res$hq,  4),
      Risk        = if(res$cdi>0.001)"Elevated" else if(res$cdi>0.0001)"Moderate" else "Acceptable",
      stringsAsFactors=FALSE)
    history(c(list(entry), h)[seq_len(min(10, length(h)+1))])
  })
  
  # ── Results pane ─────────────────────────────────────────────────────────────
  output$output_ui <- renderUI({
    if (input$calculate == 0)
      return(div(class="d-flex flex-column align-items-center justify-content-center h-100 text-muted",
                 style="min-height:260px;gap:0.6rem;",
                 bsicons::bs_icon("calculator",size="2.5rem"),
                 p("Fill parameters and press",strong("Calculate"))))
    
    res <- cdi_result()
    if (!is.null(res$error))
      return(div(class="alert alert-warning",
                 bsicons::bs_icon("exclamation-triangle")," ",res$error))
    
    cdi <- res$cdi; hq <- res$hq
    formula_str <- switch(res$mode,
                          Ingestion        ="CDI = (C \u00d7 IR \u00d7 EF \u00d7 ED) / (BW \u00d7 AT)",
                          Inhalation       ="CDI = (C \u00d7 IR \u00d7 ET \u00d7 EF \u00d7 ED) / (BW \u00d7 AT)",
                          `Skin Absorption`="CDI = (C \u00d7 SA \u00d7 AF \u00d7 ABS \u00d7 EF \u00d7 ED) / (BW \u00d7 AT)")
    threshold <- 0.001
    risk <- if (cdi > threshold)
      list(label="Elevated Risk",    class="danger", icon="exclamation-triangle-fill",
           interp=paste0("CDI of ",fmt_num(cdi,6)," mg/kg\u00b7day EXCEEDS the 0.001 threshold. ",
                         "Significant non-carcinogenic risk. Remediation and further assessment recommended."))
    else if (cdi > threshold*0.1)
      list(label="Moderate Concern", class="warning",icon="exclamation-circle-fill",
           interp=paste0("CDI of ",fmt_num(cdi,6)," mg/kg\u00b7day is below threshold but elevated. ",
                         "Monitoring and periodic reassessment advised."))
    else
      list(label="Acceptable",       class="success",icon="check-circle-fill",
           interp=paste0("CDI of ",fmt_num(cdi,6)," mg/kg\u00b7day is well below the 0.001 threshold. ",
                         "Risk is acceptable under these exposure assumptions. Continue routine monitoring."))
    bar_pct <- min(100,(cdi/threshold)*100)
    comparisons <- lapply(names(res$v), function(id){
      val <- as.numeric(res$v[[id]]); dval <- zdhs_defaults[[id]]
      if (is.null(dval)||is.na(val)) return(NULL)
      ic <- if(val>dval)"arrow-up-circle-fill" else if(val<dval)"arrow-down-circle-fill" else "check-circle-fill"
      cl <- if(val>dval)"text-warning" else if(val<dval)"text-info" else "text-success"
      tags$span(class=paste("me-2",cl),style="font-size:0.75rem;",
                bsicons::bs_icon(ic)," ",id,": ",val,
                if(val!=dval) paste0(" (ZDHS: ",dval,")") else " (= ZDHS default)")
    })
    comparisons <- Filter(Negate(is.null), comparisons)
    
    tagList(
      if(length(res$warnings)>0)
        div(class="alert alert-warning py-1 mb-2",style="font-size:0.78rem;",
            bsicons::bs_icon("exclamation-triangle")," Out-of-range: ",
            paste(res$warnings,collapse=", "),". Review before interpreting.") else NULL,
      div(class="result-box mb-3", p(class="formula-code mb-0",formula_str)),
      layout_column_wrap(width=1/2,
                         div(class="result-box",
                             p(class="text-muted mb-1",style="font-size:0.72rem;","CHRONIC DAILY INTAKE"),
                             h4(class="mb-0",fmt_num(cdi,7)),
                             p(class="text-muted mb-0",style="font-size:0.7rem;","mg \u00b7 kg\u207b\u00b9 \u00b7 day\u207b\u00b9")),
                         div(class="result-box",
                             p(class="text-muted mb-1",style="font-size:0.72rem;",
                               paste0("HAZARD QUOTIENT (RfD = ",res$rfd,")")),
                             h4(class="mb-0",fmt_num(hq,4)),
                             p(class="text-muted mb-0",style="font-size:0.7rem;",
                               paste0("CDI \u00f7 ",res$contaminant," RfD")))
      ),
      div(class="mt-3",
          div(class="d-flex align-items-center justify-content-between mb-1",
              span(class="text-muted",style="font-size:0.8rem;","Risk level"),
              span(class=paste0("risk-badge bg-",risk$class,"-subtle text-",risk$class,"-emphasis"),
                   bsicons::bs_icon(risk$icon)," ",risk$label)),
          div(class="progress",style="height:8px;",
              div(class=paste0("progress-bar bg-",risk$class),
                  style=paste0("width:",round(bar_pct,1),"%;"),role="progressbar")),
          div(class="d-flex justify-content-between mt-1",
              span(class="text-muted",style="font-size:0.7rem;","0"),
              span(class="text-muted",style="font-size:0.7rem;","Threshold (0.001)"))
      ),
      div(class="interp-box mt-3", bsicons::bs_icon("chat-text")," ",risk$interp),
      tags$hr(),
      p(class="text-muted mb-1",style="font-size:0.72rem;font-weight:600;","INPUT vs ZDHS DEFAULTS"),
      div(style="font-size:0.75rem;line-height:2;",comparisons),
      tags$hr(),
      div(class="d-flex justify-content-between align-items-center",
          p(class="text-muted mb-0",style="font-size:0.7rem;",
            "HQ < 1 = acceptable non-carcinogenic risk. Use contaminant-specific RfD from EPA IRIS."),
          downloadButton("download_report","Export CSV",
                         class="btn-outline-secondary btn-sm",icon=icon("download")))
    )
  })
  
  # ── Analysis pane (sensitivity + pathway) ────────────────────────────────────
  output$analysis_ui <- renderUI({
    if (input$calculate == 0)
      return(div(class="text-muted text-center",style="padding:2rem;",
                 bsicons::bs_icon("graph-up",size="2rem"),
                 p("Run a calculation first to see analysis.")))
    res <- cdi_result()
    if (!is.null(res$error)) return(div(class="alert alert-warning",res$error))
    
    # Sensitivity: vary each param ±10%, compute % change in CDI
    base_cdi <- res$cdi
    sens <- lapply(names(res$v), function(id){
      v_up <- res$v; v_up[[id]] <- res$v[[id]] * 1.1
      v_dn <- res$v; v_dn[[id]] <- res$v[[id]] * 0.9
      cdi_up <- calc_cdi(res$mode, lapply(v_up, as.numeric))
      cdi_dn <- calc_cdi(res$mode, lapply(v_dn, as.numeric))
      pct_change <- if (!is.na(cdi_up) && base_cdi != 0)
        abs((cdi_up - cdi_dn) / (2 * base_cdi)) * 100 else 0
      list(id=id, pct=round(pct_change, 1))
    })
    sens <- sens[order(sapply(sens,"[[","pct"), decreasing=TRUE)]
    max_pct <- max(sapply(sens,"[[","pct"), 1)
    
    sens_bars <- lapply(sens, function(s){
      bar_w <- round(s$pct / max_pct * 100)
      div(class="sens-bar-wrap",
          span(class="sens-label", s$id),
          div(class="sens-track",
              div(class="sens-pos", style=paste0("width:",bar_w,"%;")))  ,
          span(class="sens-pct", paste0(s$pct,"%"))
      )
    })
    
    # Exposure pathway text diagram
    pathway <- switch(res$mode,
                      Ingestion        = "Source \u2192 Water / Food \u2192 Ingestion \u2192 GI Tract \u2192 Systemic circulation",
                      Inhalation       = "Source \u2192 Air \u2192 Inhalation \u2192 Lungs \u2192 Systemic circulation",
                      `Skin Absorption`= "Source \u2192 Soil / Surface \u2192 Skin contact \u2192 Dermal absorption \u2192 Systemic circulation"
    )
    
    tagList(
      tags$h6(class="fw-bold mb-3", bsicons::bs_icon("sliders")," Sensitivity Analysis"),
      p(style="font-size:0.8rem;color:var(--bs-secondary-color);margin-bottom:1rem;",
        "Each bar shows the % change in CDI when that parameter is varied by \u00b110%.
         Taller bars = greater influence on the result."),
      div(sens_bars),
      tags$hr(class="my-3"),
      tags$h6(class="fw-bold mb-2", bsicons::bs_icon("diagram-3")," Exposure Pathway"),
      div(class="result-box",
          p(style="font-family:monospace;font-size:0.82rem;margin:0;",pathway)),
      tags$hr(class="my-3"),
      tags$h6(class="fw-bold mb-2", bsicons::bs_icon("toggles")," Carcinogenic Risk"),
      p(style="font-size:0.8rem;color:var(--bs-secondary-color);",
        "For carcinogenic contaminants, risk = CDI \u00d7 Cancer Slope Factor (CSF).
         Enter CSF below (from EPA IRIS) to calculate lifetime cancer risk."),
      layout_column_wrap(width=1/2,
                         numericInput("csf","CSF (mg/kg\u00b7day)\u207b\u00b9",value=NA,min=0,step=NA,width="100%"),
                         div(class="result-box",style="margin-top:1.5rem;",
                             p(class="text-muted mb-1",style="font-size:0.72rem;","CANCER RISK"),
                             h5(class="mb-0", uiOutput("cancer_risk_val")),
                             p(class="text-muted mb-0",style="font-size:0.7rem;","Acceptable if < 1\u00d710\u207b\u2076"))
      )
    )
  })
  
  output$cancer_risk_val <- renderUI({
    res <- cdi_result()
    if (is.null(res) || !is.null(res$error) || is.na(input$csf) || input$csf <= 0)
      return(span("—"))
    cr <- res$cdi * input$csf
    span(format(cr, scientific=TRUE, digits=3))
  })
  
  # ── History pane ──────────────────────────────────────────────────────────────
  output$history_ui <- renderUI({
    h <- history()
    if (length(h)==0)
      return(div(class="text-muted text-center",style="padding:2rem;",
                 bsicons::bs_icon("clock-history",size="2rem"),
                 p("No calculations yet. Results will appear here.")))
    df <- do.call(rbind, h)
    tagList(
      p(style="font-size:0.8rem;color:var(--bs-secondary-color);",
        "Last ",nrow(df)," calculations (most recent first)."),
      div(style="overflow-x:auto;",
          tags$table(class="table table-sm table-hover history-table",
                     tags$thead(tags$tr(lapply(names(df), tags$th))),
                     tags$tbody(
                       apply(df, 1, function(row){
                         risk_class <- if(row["Risk"]=="Elevated") "table-danger"
                         else if(row["Risk"]=="Moderate") "table-warning"
                         else "table-success"
                         tags$tr(class=risk_class, lapply(row, tags$td))
                       })
                     )
          )
      ),
      div(class="d-flex justify-content-end mt-2",
          downloadButton("download_history","Export History CSV",
                         class="btn-outline-secondary btn-sm",icon=icon("download")))
    )
  })
  
  # ── Downloads ─────────────────────────────────────────────────────────────────
  output$download_report <- downloadHandler(
    filename=function() paste0("CDI_report_",format(Sys.Date(),"%Y%m%d"),".csv"),
    content=function(file){
      res <- cdi_result()
      if (is.null(res)||!is.null(res$error)){
        write.csv(data.frame(Error=res$error),file,row.names=FALSE); return()
      }
      df <- rbind(
        data.frame(Parameter=names(res$v), Value=unlist(res$v), stringsAsFactors=FALSE),
        data.frame(Parameter="Exposure Mode",        Value=res$mode),
        data.frame(Parameter="Contaminant",          Value=res$contaminant),
        data.frame(Parameter="RfD (mg/kg/day)",      Value=res$rfd),
        data.frame(Parameter="CDI (mg/kg/day)",      Value=round(res$cdi,8)),
        data.frame(Parameter="Hazard Quotient (HQ)", Value=round(res$hq,6)),
        data.frame(Parameter="Risk Classification",
                   Value=if(res$cdi>0.001)"Elevated Risk" else if(res$cdi>0.0001)"Moderate Concern" else "Acceptable"),
        data.frame(Parameter="Date",  Value=as.character(Sys.Date())),
        data.frame(Parameter="Tool",  Value=paste("EHRA App",APP_VERSION))
      )
      write.csv(df,file,row.names=FALSE)
    }
  )
  
  output$download_history <- downloadHandler(
    filename=function() paste0("CDI_history_",format(Sys.Date(),"%Y%m%d"),".csv"),
    content=function(file){
      h <- history()
      if(length(h)==0){ write.csv(data.frame(Message="No history"),file,row.names=FALSE); return() }
      write.csv(do.call(rbind,h),file,row.names=FALSE)
    }
  )
}

shinyApp(ui, server)
