pacman::p_load(
    rio,
    tidyverse,
    performance,
    here,
    gt,
    gtsummary,
    flextable,
    see,
    report,
    AER,
    DHARMa,
    sjPlot,
    epiR
)


path <- "D:/Field Epidemiology 2026/Introduction to Biostatistics II/Poisson Regression/Intro"
couart <- import(here(path, "couart910.dta"))

data_phd <- couart |> 
    mutate(
        fem = if_else(
            fem == 1, "Female", "Male",
        ),
        
        mar = if_else(
            mar == 1, "Yes", "No"
        )
    )


plot_frq(
    data_phd$art,
    type = "histogram",
    normal.curve = T,
    show.mean = T, 
    show.sd = T
)+
    theme(
        panel.background = element_blank()
    )

plot_frq(
    data_phd$kid5,
    type = "histogram",
    normal.curve = T,
    show.mean = T, 
    show.sd = T
)+
    theme(
        panel.background = element_blank()
    )

 plot_grpfrq(data_phd$art, 
            data_phd$mar,
            type = "violin")+
    theme(
        panel.background = element_blank()
    )+
     coord_flip()



 ggplot(
    data_phd,
    aes(x = art)
)+
    geom_histogram(
        aes(y =after_stat(density)),
        color  = "navy",
        fill = "steelblue",
        binwidth = 1,
        alpha = 0.3
    )+
    
    stat_function(
        fun = dnorm,
        args = list(
            mean = mean(data_phd$art, na.rm = T),
            sd = sd(data_phd$art, na.rm = T)
        ),
        color = "navy",
        size = 1
    )+
    theme(
        panel.background = element_blank()
    )




poistype = pois_phd <- glm(
    art ~ .,
    family = poisson,
    data = data_phd
)



tbl_regression(
    pois_phd,
    exponentiate = T,
    label = list(
        fem = "Gender",
        mar = "Married",
        kid5 = "Number of Children"
    ),
    pvalue_fun = label_style_pvalue(digits = 3)
) |>  bold_labels() |> 
    bold_p() |> 
    as_flex_table()

summary(data_phd$art)

hist(data_phd$art, 
     col = "steelblue",
     main = "Distribution of the number of Articles")


ggplot(data_phd)+
    geom_histogram (
        aes(art
            
            ),
        binwidth = 1,
        fill = "lavender",
        color = "black"
    )+
    theme(plot.background = element_blank(),
          panel.background = element_blank()
          )+
    labs(
        y = "Number of Articles",
        x = " ",
        title = "Number of Articles"
    )
dev.new() 
check_overdispersion(pois_phd)
model_performance(pois_phd)
check_collinearity(pois_phd)
report_effectsize(pois_phd)
check_convergence(pois_phd)
check_zeroinflation(pois_phd)
check_heteroscedasticity(pois_phd)
check_model(pois_phd)
check_zeroinflation(pois_phd)

anova(pois_phd)
summ <- summary(pois_phd)


over_dispersion <- round(summ$deviance/summ$df.residual, 3)

cat("Dispersion Parameter :", 
    over_dispersion, "\n")


AER::dispersiontest(pois_phd)



data_phd %>%
    group_by(mar ) %>%
    summarise(
        w_stat = shapiro.test(art)$statistic,
        p_val = format(shapiro.test(art)$p.value, scientific = F)
    ) |> 
    gt()


