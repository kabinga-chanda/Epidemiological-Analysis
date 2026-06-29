
# Load Libraries ----------------------------------------------------------

pacman::p_load(
    R2MLwiN,
    lme4,
    epiDisplay,
    sjPlot,
    tidyverse,
    naniar,
    labelled,
    gt,
    performance,
    gtsummary
)


# Load Data ---------------------------------------------------------------

data(bang)

# Data Preparations -------------------------------------------------------

data_clean <- bang |> 
    set_variable_labels(
        woman    = "Woman ID",
        district = "District ID",
        use      = "Contraceptive Use",
        use4     = "Contraceptive Method Choice",
        lc       = "Number of Living Children",
        age      = "Age of Woman (Mean-Centered)",
        urban    = "Type of Area",
        educ     = "Education Level",
        hindu    = "Religion",
        d_lit    = "District Literacy Rate",
        d_pray   = "District Religiosity",
        cons     = "Intercept Constant"
    ) |> 
    mutate(
        "contraceptive_use" = if_else(
            use == "Using", 1, 0
        ),
        
        lc = factor(
            case_when(
                lc == "One_child" ~ "One Child",
                lc == "Two_children" ~ "Two Children",
                lc ==  "Three_plus" ~  "More than 3 Children",
                lc == "None" ~ "No Children",
                TRUE ~ NA_character_
            ),
            levels = c("No Children",
                       "One Child",
                       "Two Children",
                       "More than 3 Children"
                       )
            
        ),
        
        use4 = factor(
            case_when(
                use4 == "Not_using_contraception" ~ "Not Using",
                use4 == "Sterilization"           ~ "Sterilization",
                use4 == "Modern_reversible_method" ~ "Modern Reversible",
                use4 == "Traditional_method"      ~ "Traditional",
                TRUE                              ~ NA_character_
            ),
            levels = c(
                "Not Using",
                "Sterilization",
                "Modern Reversible",
                "Traditional"
            )
        ),
        
        educ = factor(
            case_when(
                educ == "None"                ~ "No Education",
                educ == "Lower_primary"       ~ "Lower Primary",
                educ == "Upper_primary"       ~ "Upper Primary",
                educ == "Secondary_and_above" ~ "Secondary and Above",
                TRUE                          ~ NA_character_
            ),
            levels = c(
                "No Education",
                "Lower Primary",
                "Upper Primary",
                "Secondary and Above"
            )
        )
    )



# Table of Descriptives ---------------------------------------------------

data_clean |> 
    select(-c(woman,
              district,
              cons,
              use)) |> 
    mutate(
        contraceptive_use = factor(
            if_else(
                contraceptive_use == 1, "Using", "Not Using"
            ), 
            
            levels = c("Not Using", "Using")
        ),
        
        age = age + 30
    ) |> 
    tbl_summary(
        by = contraceptive_use,
        label = list(
            age      = "Age (Years)",
            urban    = "Type of Area",
            educ     = "Education Level",
            hindu    = "Religion",
            use4     = "Contraceptive Choice",
            lc       = "Number of Living Children"
            
        )
    ) |> 
    bold_labels() |> 
    add_p(
        use4 ~ fisher.test
    ) |> 
    bold_p() |> 
    as_gt() |> 
    tab_header(
        title = md("**Characteristics And Demographics**"),
        subtitle = md("Sub-sample from the 1989 Bangladesh <br> Fertility Survey (see Huq & Cleveland, 1990)")
    ) |> 
    tab_source_note(
        source_note = md("**Source:** 1989 Bangladesh Fertility Survey (Huq and Cleland, 1990)")
    ) |> 
    tab_options(
        # 1. Clear internal row lines
        table_body.hlines.color = "transparent",
        
        # 2. Set the top line (above the column labels) to black
        heading.border.bottom.color = "black",
        
        # 3. Set the line separating the column labels from the data to black
        column_labels.border.top.color = "black",
        column_labels.border.bottom.color = "black",
        
        # 4. Set the very bottom line of the data body to black
        table_body.border.bottom.color = "black"
    )




# Logistic Regression -----------------------------------------------------

intervention_model <- glmer(
    contraceptive_use ~ urban + age + lc + (1 | district), 
    data = data_clean, 
    family = binomial(link = "logit")
)


tab_model(
    intervention_model, 
    title = "Predictors of Contraceptive Use, Bangladesh",
    pred.labels = c(
        "(Intercept)" = "Baseline Intercept",
        "urban"       = "Type of Area (Urban)",
        "age"         = "Age of Woman",
        "lcOne Child" = "Living Children: One Child",
        "lcTwo Children" = "Living Children: Two Children",
        "lcMore than 3 Children" = "Living Children: 3+ Children"
    ),
    show.p = F,
    show.se = T
)


# Diagnose the Model ------------------------------------------------------





