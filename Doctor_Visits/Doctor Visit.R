

pacman::p_load(gt, gtsummary, MASS, AER, tidyverse, stringr, sjPlot, ggExtra)


data <- MASS::Insurance

DoctorVisits <- DoctorVisits |> 
  mutate(across(c(gender, private, freerepat, freepoor, nchronic, lchronic),
                ~ str_to_title(.x)))

DoctorVisits$illness <- round(DoctorVisits$illness, 1)




data("DoctorVisits")
?DoctorVisits






DoctorVisits <- DoctorVisits |> 
  mutate(across(c(gender, private, freerepat, freepoor, nchronic, lchronic),
                ~ str_to_title(.x)))



# Data Summary ------------------------------------------------------------

DoctorVisits_1 <- DoctorVisits |> 
  mutate(visit_stat = if_else(visits>0, "Visted", "Didn't Vist")) |> 
  select(-visits)



tbl_summary(DoctorVisits_1,
            by = visit_stat,
            include = everything(),
            type = list(illness ~ "continuous"),
            statistic = list(all_continuous() ~ "{median} ({p25}, {p75})",
                             all_categorical() ~  "{n} ({p}%)"),
            percent = "row",
            label = list(gender ~ "Gender",
                         age ~ "Age(In Years)",
                         income ~ "Annual income",
                         illness = "Number of illnesses in past 2 weeks",
                         reduced = "Number of days of reduced activity in past 2 weeks",
                         health ~ "General health",
                         private ~ "Private health insurance",
                         freepoor ~ "Free government health insurance",
                         nchronic ~ "Chronic condition not limiting activity",
                         lchronic ~ " Chronic condition limiting activity")) |> 
  bold_labels() |> 
  add_p() |> 
  add_overall() |> 
  bold_p() |> as_gt() |> 
  tab_options(
    table_body.hlines.style = "none") |> 
  tab_header(
    title = md("**Characterics & Demographics of Doctor visits in Australia**"),
    subtitle = "Poisson regression results from the 1977-1978 Australian Health Survey") |> 
  opt_align_table_header(align = "left")


# Models ------------------------------------------------------------------



model_pois <- glm(visits~., data = DoctorVisits, family = poisson)
dispersiontest(model_pois)
model_ng <- glm.nb(visits~., data = DoctorVisits)
model_quis <- glm(visits~., data = DoctorVisits, family = quasipoisson())




m1 <- tbl_regression(model_pois,
               exponentiate = T,
               label = list(gender ~ "Gender",
                            age ~ "Age(In Years)",
                            income ~ "Annual income",
                            illness = "Number of illnesses in past 2 weeks",
                            reduced = "Number of days of reduced activity in past 2 weeks",
                            health ~ "General health",
                            private ~ "Private health insurance",
                            freepoor ~ "Free government health insurance",
                            nchronic ~ "Chronic condition not limiting activity",
                            lchronic ~ " Chronic condition limiting activity"))|> 
  bold_labels() |> 
  bold_p()|>   
  add_glance_table(include = c(nobs, logLik, AIC, BIC))


m2 <- tbl_regression(model_ng,
               exponentiate = T,
               label = list(gender ~ "Gender",
                            age ~ "Age(In Years)",
                            income ~ "Annual income",
                            illness = "Number of illnesses in past 2 weeks",
                            reduced = "Number of days of reduced activity in past 2 weeks",
                            health ~ "General health",
                            private ~ "Private health insurance",
                            freepoor ~ "Free government health insurance",
                            nchronic ~ "Chronic condition not limiting activity",
                            lchronic ~ " Chronic condition limiting activity"))|> 
  bold_labels() |> 
  bold_p()|>   
  add_glance_table(include = c(nobs, logLik, AIC, BIC))

m3 <- tbl_regression(model_quis, 
               exponentiate = T,
               label = list(gender ~ "Gender",
                            age ~ "Age(In Years)",
                            income ~ "Annual income",
                            illness = "Number of illnesses in past 2 weeks",
                            reduced = "Number of days of reduced activity in past 2 weeks",
                            health ~ "General health",
                            private ~ "Private health insurance",
                            freepoor ~ "Free government health insurance",
                            nchronic ~ "Chronic condition not limiting activity",
                            lchronic ~ " Chronic condition limiting activity"))|> 
  bold_labels() |> 
  bold_p()|>   
  add_glance_table(include = c(nobs, logLik, AIC, BIC))




tbl_merge(
  tbls = list(m1,m2,m3),
  tab_spanner = c("**Poisson**", "**Negative Binomial**", "**Quasipoisson**")) |> 
  as_gt() |> 
  tab_header(
    title = md("**Determinants of Doctor Visits in Australia**"),
    subtitle = "Poisson regression results from the 1977-1978 Australian Health Survey") |> 
   tab_source_note(source_note = md("**Source:** Journal of Applied Econometrics Data Archive. *http://qed.econ.queensu.ca/jae/1997-v12.3/mullahy/*")) |> 
  tab_footnote(footnote = "The results indicate that overdispersion is present in the Poisson Model") |> 
  tab_options(table_body.hlines.style = "none") |> 
  opt_align_table_header(align = "left") 



ggplot(DoctorVisits) +
  geom_histogram(aes(x = visits), bins = 50)



pacman::p_load(pscl)

# Fit a zero-inflated Poisson (ZIP) model
zip_model <- zeroinfl(visits ~ age + gender + income + illness + reduced + health + private + freepoor + nchronic + lchronic | 1,
                      data = DoctorVisits,
                      dist = "poisson")




vuong(model_pois, zip_model)

DoctorVisits |> 
ggplot()+
  geom_density(aes(x = visits),
               stat = "density",
               fill = "red")





# Assumption Check --------------------------------------------------------

# 1. Equidispersion
dispersiontest(model_pois)

# 2. Independence (visual check)
acf(residuals(model_pois, type = "pearson"))

# 3. Linearity
car::residualPlots(model_pois)


# 4. Model specification
car::ceresPlots(model_pois)

# 5. Zero-inflation
testZeroInflation(simulateResiduals(model_pois))

# 6. Overall fit
plot(simulateResiduals(model_pois))
