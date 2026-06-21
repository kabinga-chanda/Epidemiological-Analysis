library(survival)
library(survminer)
library(tidyverse)
library(gtsummary)
library(gt)
library(flextable)
library(rio)



heroin <- import("D:/Field Epidemiology 2026/Introduction to Biostatistics II/heroin.dta")
couart <- import("D:/Field Epidemiology 2026/Introduction to Biostatistics II/couart910.dta")



# Heroin ------------------------------------------------------------------

fit <- surv_fit(Surv(time, status) ~ clinic, data = heroin)
    ggsurvplot(fit,
           risk.table = F,
           pval.method = F,
           pval = F,
           palette = c("steelblue",
                       "red"))


log_rank <- survdiff(Surv(time, status) ~ clinic, data = heroin)
cat("The Log-rank: ", format(log_rank$pvalue, scientific = F), "X-stat: ", log_rank$chisq, "\n")

cox <- coxph(Surv(time, status) ~  dose + prison,
             data = heroin )



tbl_regression(
    cox,
    exponentiate = T
) |> bold_p() |> 
    bold_labels() |> 
    as_flex_table()

ggplot(heroin) +
    geom_bar(
        aes(time)
    )


ggforest(cox)

pairwise_survdiff(Surv(time, status) ~ clinic, data = heroin)



# Poisson Regression ------------------------------------------------------

pois <- glm(art ~ kid5+mar+phd+ment, family = poisson, data = couart)
tbl_regression(
    pois,
    exponentiate = T
) |> 
    bold_labels() |> 
    bold_p() |> 
    as_flex_table()

# Check for Missingness
heroin |> 
    group_by(clinic) |> 
    summarise(
        missingness = sum(is.na(clinic)),
        missingness_pct = (sum(is.na(clinic))/ nrow(heroin))*100
    )


