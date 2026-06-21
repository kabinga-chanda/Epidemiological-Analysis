
# Load the Libraries ------------------------------------------------------

pacman::p_load(
    rio,
    tidyverse,
    performance,
    here,
    gt,
    gtsummary,
    flextable,
    see
)


# Load the Data -----------------------------------------------------------

polyps <- data.frame(
    male = c(0,0,1,0,1,0,0,1,1,1,1,0,1,1,1,1,0,1,1,0,0,1),
    trt = c(1,0,1,0,1,0,1,0,0,0,1,1,0,0,0,1,0,1,0,1,1,1),
    baspolyps = c(7,77,7,5,23,35,11,12,7,318,160,8,20,11,24,34,54,16,30,10,20,12),
    age = c(17,20,16,18,22,13,23,34,50,19,17,23,22,30,27,23,22,13,34,23,22,42),
    postpolyps = c(6,67,4,5,16,31,6,20,7,347,142,1,16,20,26,27,45,10,30,6,5,8)
)

polyps <- polyps |> 
    mutate(
        
        male = if_else(
            male == 1, "Famale", "Male"
        ),
        trt = if_else(
            trt == 1, "Active", "Placebo"
        )
    )


polyps |> 
    tbl_summary(
        by = trt,
        label = list(
            male       = "Gender",
            trt        = "Treatment",
            age        = "Age",
            baspolyps  = "Baseline Polyps",
            postpolyps = "Post Polyps"
        )
    ) |> 
    add_p(
        test = list(
            all_categorical() ~ "fisher.test",
            all_continuous()  ~ "wilcox.test"
        ),
        pvalue_fun = label_style_pvalue(digits = 3)
    ) |> 
    bold_p() |> 
    bold_labels()|> 
    as_gt() |> 
    tab_header (
        title = md("**Table 1. Baseline Characteristics by Treatment Group**")
    ) 
    




# Poisson Regression ------------------------------------------------------

muv_poison <- glm(
    postpolyps ~ male + trt + age+baspolyps,
    family = poisson,
    data = polyps
) 



muv_gt <- tbl_regression(
    muv_poison,
    exponentiate = T,
    label = list(
        male  = "Gender",
        trt = "Treatment",
        age = "Age (Years)",
        baspolyps = "Baseline Post-Polyps"
    ),
    pvalue_fun = label_style_pvalue(digits = 3)
) |> 
    bold_p() |> 
    bold_labels()

univ_poison <- tbl_uvregression(
    data = polyps,
    y = postpolyps,
    include = c(male, trt, age,baspolyps),
    method = glm,
    method.args = list(family = poisson),
    exponentiate = TRUE,
    label = list(
        male  = "Gender",
        trt = "Treatment",
        age = "Age (Years)",
        baspolyps = "Baseline Post-Polyps"
    ),
    pvalue_fun = label_style_pvalue(digits = 3)
) |> 
    bold_p() |> 
    bold_labels()

tbl_merge(
    tbls = list(univ_poison, muv_gt),
    tab_spanner = c("**Univariable**", "**Multi-Variable**"
) )|> as_gt() |> 
    tab_header(
        md("**Predictors of Post-Polyps**")
    )

check_overdispersion(muv_poison)





