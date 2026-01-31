### Load Libraries
library(tidyverse)
library(rio)
library(janitor)
library(gt)
library(gtsummary)
library(tmap)
library(sf)
library(ggpubr)
library(naniar)
library(tmap)
library(rnaturalearth)
library(rnaturalearthdata)
library(sjPlot)


### Load the dataset
data_age_stand <- import("D:/DSO_mac/R Portifolio/Suicide Cases/data_standardized_aga.xlsx")
data <- import("D:/DSO_mac/R Portifolio/Suicide Cases/data_country_age.xlsx")
alco_rate <- import("D:/DSO_mac/R Portifolio/Suicide Cases/alco_rate.xls")
urban_rate <- import("D:/DSO_mac/R Portifolio/Suicide Cases/urban_rate.xls")
unemploy_rate <- import("D:/DSO_mac/R Portifolio/Suicide Cases/unemploy_rate.xls")


## Data cleaning
data_clean <- data |>
    filter(Period == 2020) |>
    select("country" = Location,
           "continent" = ParentLocation,
           "year" = Period,
           "sex" = Dim1,
           "suicide_rate" = FactValueNumericHigh)



### Data_age_standardize
data_age_standardize <- data_age_stand |>
    filter(Period == "2020") |>
    select(
        country = Location,
        suicide_rate = `Both sexes`,
        year = Period,
        male = Male,
        female = Female
    ) |>
    mutate(
        suicide_rate = as.numeric(str_extract(suicide_rate, "^[0-9.]+")),
        male = as.numeric(str_extract(male, "^[0-9.]+")),
        female = as.numeric(str_extract(female, "^[0-9.]+"))
    ) |>
    pivot_longer(cols = c(female, male),
                 names_to = "sex",
                 values_to = "suicide_sex")

### Data for poisson regression
data_pois_outcome <- data_age_stand |>
    filter(Period == "2020") |>
    select(
        country = Location,
        suicide_rate = `Both sexes`,
        year = Period
    ) |>
    mutate(
        suicide_rate = as.numeric(str_extract(suicide_rate, "^[0-9.]+"))
    )



## Clean Alcohol Consumption Data
data_alcoh_rate <- alco_rate |>
    select(
        "country" = `Country Name`,
        "alco_consu" = `2020`
    )


## Clean Urbanization
data_urban_rate <- urban_rate |>
    select(
        "country" = `Country Name`,
        "urban_rate" = `2020`
    )


## Unemployement Rate
data_unemploy_rate <- unemploy_rate |>
    select(
        "country" = `Country Name`,
        "unemploy_rate" = `2020`
    )

### plot
data_age_standardize |>
    mutate(sex = if_else(sex == "female", "Female", "Male")) |>  # Fixed: "Famale" → "Female"
    tbl_summary(
        by = sex,
        include = c(suicide_sex),  # Removed 'sex' from include
        label = list(
            suicide_sex ~ "Age Standardized Suicide Rate"
        )
    ) |>
    add_p() |>
    bold_p() |>
    bold_labels() |>
    as_gt() |>  # Fixed: gt() → as_gt()
    tab_header(title = md("**Age Standardized Suicide Rate Per 100,000 by Sex**"))  # Added closing


###Merge all the explanatory variables
predictors <- data_alcoh_rate |>
    left_join(data_unemploy_rate, by = "country") |>
    left_join(data_urban_rate, by = "country")

## Impoutate Missing data
predictors_clean <- predictors |>
    drop_na()

### Standardise the country names
regional_groups <- c(
    "Africa Eastern and Southern", "Africa Western and Central", "Arab World",
    "Central Europe and the Baltics", "Caribbean small states",
    "East Asia & Pacific (excluding high income)", "Early-demographic dividend",
    "East Asia & Pacific", "Europe & Central Asia (excluding high income)",
    "Europe & Central Asia", "Euro area", "European Union",
    "Fragile and conflict affected situations", "High income",
    "Heavily indebted poor countries (HIPC)", "IBRD only", "IDA & IBRD total",
    "IDA total", "IDA blend", "IDA only",
    "Latin America & Caribbean (excluding high income)",
    "Latin America & Caribbean", "Least developed countries: UN classification",
    "Low income", "Lower middle income", "Low & middle income",
    "Late-demographic dividend", "Middle East, North Africa, Afghanistan & Pakistan",
    "Middle income", "Middle East, North Africa, Afghanistan & Pakistan (excluding high income)",
    "North America", "OECD members", "Other small states", "Pre-demographic dividend",
    "Pacific island small states", "Post-demographic dividend", "South Asia",
    "Sub-Saharan Africa (excluding high income)", "Sub-Saharan Africa", "Small states",
    "East Asia & Pacific (IDA & IBRD countries)",
    "Europe & Central Asia (IDA & IBRD countries)",
    "Latin America & the Caribbean (IDA & IBRD countries)",
    "Middle East, North Africa, Afghanistan & Pakistan (IDA & IBRD)",
    "South Asia (IDA & IBRD)", "Sub-Saharan Africa (IDA & IBRD countries)",
    "Upper middle income", "World"
)


predictors_clean <- predictors_clean |>
    filter(!country %in% regional_groups) |>
    mutate(
        country = case_when(
            # Match to WHO reference names
            country == "Bahamas, The" ~ "Bahamas",
            country == "Bolivia" ~ "Bolivia (Plurinational State of)",
            country == "Brunei Darussalam" ~ "Brunei Darussalam",
            country == "Congo, Dem. Rep." ~ "Democratic Republic of the Congo",
            country == "Congo, Rep." ~ "Congo",
            country == "Egypt, Arab Rep." ~ "Egypt",
            country == "Gambia, The" ~ "Gambia",
            country == "Iran, Islamic Rep." ~ "Iran (Islamic Republic of)",
            country == "Korea, Rep." ~ "Republic of Korea",
            country == "Korea, Dem. People's Rep." ~ "Democratic People's Republic of Korea",
            country == "Kyrgyz Republic" ~ "Kyrgyzstan",
            country == "Lao PDR" ~ "Lao People's Democratic Republic",
            country == "Moldova" ~ "Republic of Moldova",
            country == "Netherlands" ~ "Netherlands (Kingdom of the)",
            country == "Russian Federation" ~ "Russian Federation",
            country == "Slovak Republic" ~ "Slovakia",
            country == "Somalia, Fed. Rep." ~ "Somalia",
            country == "St. Lucia" ~ "Saint Lucia",
            country == "St. Vincent and the Grenadines" ~ "Saint Vincent and the Grenadines",
            country == "Syrian Arab Republic" ~ "Syrian Arab Republic",
            country == "Tanzania" ~ "United Republic of Tanzania",
            country == "Turkiye" ~ "Türkiye",
            country == "United Kingdom" ~ "United Kingdom of Great Britain and Northern Ireland",
            country == "United States" ~ "United States of America",
            country == "Venezuela, RB" ~ "Venezuela (Bolivarian Republic of)",
            country == "Viet Nam" ~ "Viet Nam",
            country == "Yemen, Rep." ~ "Yemen",
            TRUE ~ country  # Keep all other names as is
        )
    )

### Merge the outcome and predictors
merged_data_clean <- data_pois_outcome |>
    inner_join(predictors_clean, by = "country")


## violin plot
data_age_standardize |>
    mutate(sex = if_else(sex == "female", "Female", "Male")) |>
    ggplot(aes(x = sex, y = suicide_sex, fill = sex)) +
    geom_violin(alpha = 0.6, width = 0.8) +
    geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
    stat_compare_means(
        method = "wilcox.test",
        label = "p.format",
        size = 5,
        label.y.npc = 0.95
    ) +
    scale_fill_manual(values = c("Female" = "#E74C3C", "Male" = "#3498DB")) +
    labs(
        title = "Age-Standardized Suicide Rate by Sex",
        subtitle = "Comparison across countries (2020)",
        x = NULL,
        y = "Age-Standardized Suicide Rate (per 100,000)",
        caption = "Statistical test: Wilcoxon rank-sum test"
    ) +
    theme_minimal(base_size = 13) +
    theme(
        plot.title = element_text(face = "bold", size = 16, margin = margin(b = 5)),
        plot.subtitle = element_text(color = "grey40", size = 11, margin = margin(b = 15)),
        axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
        legend.position = "none",
        panel.grid = element_blank(),
        plot.margin = margin(20, 20, 20, 20),
        plot.background = element_blank()
    )

# Top 10 countries with highest suicide rates
data_age_standardize |>
    arrange(desc(suicide_rate)) |>
    head(10) |>
    ggplot(aes(x = reorder(country, suicide_rate), y = suicide_rate)) +
    geom_col(fill = "#E74C3C", alpha = 0.8) +
    coord_flip() +
    labs(
        title = "Top 10 Countries with Highest Suicide Rates",
        subtitle = "Age-standardized rates (2020)",
        x = NULL,
        y = "Suicide Rate (per 100,000 population)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", size = 15),
        axis.title.x = element_text(face = "bold"),
        panel.grid.major.y = element_blank()
    )


# Compare top and bottom 10
top_bottom <- bind_rows(
    data_age_standardize |> arrange(desc(suicide_rate)) |> head(10) |> mutate(group = "Highest"),
    data_age_standardize |> arrange(suicide_rate) |> head(10) |> mutate(group = "Lowest")
)

top_bottom |>
    ggplot(aes(x = reorder(country, suicide_rate), y = suicide_rate, fill = group)) +
    geom_col(alpha = 0.8) +
    coord_flip() +
    scale_fill_manual(values = c("Highest" = "#E74C3C", "Lowest" = "#27AE60")) +
    labs(
        title = "Countries with Highest and Lowest Suicide Rates",
        x = NULL,
        y = "Suicide Rate (per 100,000)",
        fill = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", size = 15),
        legend.position = "top",
        panel.grid.major.y = element_blank()
    )


### Spatial Epidemiology
# Load required libraries

# Step 1: Get world map data
world <- ne_countries(scale = "medium", returnclass = "sf")

# Step 2: Standardize country names in your data to match the map
merged_data_clean <- merged_data_clean |>
    mutate(
        country_map = case_when(
            country == "United States of America" ~ "United States",
            country == "United Kingdom of Great Britain and Northern Ireland" ~ "United Kingdom",
            country == "Iran (Islamic Republic of)" ~ "Iran",
            country == "Bolivia (Plurinational State of)" ~ "Bolivia",
            country == "Venezuela (Bolivarian Republic of)" ~ "Venezuela",
            country == "Republic of Korea" ~ "South Korea",
            country == "Democratic People's Republic of Korea" ~ "North Korea",
            country == "Lao People's Democratic Republic" ~ "Laos",
            country == "Syrian Arab Republic" ~ "Syria",
            country == "Republic of Moldova" ~ "Moldova",
            country == "Russian Federation" ~ "Russia",
            country == "Türkiye" ~ "Turkey",
            country == "Viet Nam" ~ "Vietnam",
            country == "United Republic of Tanzania" ~ "Tanzania",
            country == "Democratic Republic of the Congo" ~ "Democratic Republic of the Congo",
            country == "Congo" ~ "Republic of the Congo",
            country == "Cote d'Ivoire" ~ "Ivory Coast",
            country == "Czechia" ~ "Czech Republic",
            country == "Netherlands (Kingdom of the)" ~ "Netherlands",
            country == "occupied Palestinian territory, including east Jerusalem" ~ "Palestine",
            TRUE ~ country
        )
    )

# Step 3: Join your data with the world map
map_data <- world |>
    left_join(
        merged_data_clean |> select(country_map, suicide_rate),
        by = c("name" = "country_map")
    )

# Step 4: Create the map
tm_shape(map_data) +
    tm_polygons(
        "suicide_rate",
        palette = "-RdYlGn",  # Reversed: red = high, green = low
        style = "pretty",
        n = 6,
        border.col = "white",
        border.alpha = 0.3,
        lwd = 0.5,
        title = "Suicide Rate\n(per 100,000)",
        textNA = "No data",
        colorNA = "grey85",
        legend.reverse = TRUE
    ) +
    tm_layout(
        main.title = "Global Suicide Rates by Country (2020)",
        main.title.position = "center",
        main.title.size = 1.3,
        main.title.fontface = "bold",
        legend.position = c("left", "center"),
        legend.bg.color = "white",
        legend.bg.alpha = 0.9,
        legend.frame = TRUE,
        legend.title.size = 1.1,
        legend.text.size = 0.8,
        frame = FALSE,
        bg.color = "#f0f0f0"
    ) +
    tm_credits(
        "Source: WHO | Age-standardized rates per 100,000 population",
        position = c("center", "bottom"),
        size = 0.65,
        fontface = "italic"
    )

# Optional: Check how many countries matched
matched_countries <- sum(!is.na(map_data$suicide_rate))
total_countries <- nrow(merged_data_clean)
cat(sprintf("Matched %d out of %d countries in your data (%.1f%%)\n",
            matched_countries, total_countries,
            matched_countries/total_countries*100))

# Optional: Save the map as a high-quality image
# tmap_save(filename = "suicide_rates_world_map.png",
#           width = 12, height = 6, dpi = 300)



### Regression
model_data <- merged_data_clean %>%
    select(country, suicide_rate, alco_consu, unemploy_rate, urban_rate) %>%
    filter(suicide_rate > 0) %>%
    na.omit()

# Rename columns for presentation-quality labels
model_data <- model_data %>%
    rename(
        `Alcohol consumption (litres)` = alco_consu,
        `Unemployment rate (%)` = unemploy_rate,
        `Urban population (%)` = urban_rate
    )

cat("Final sample size:", nrow(model_data), "countries\n")

# ============================================================
# Step 1: Fit Gamma GLM Models
# ============================================================

# Multivariable model
gamma_model <- glm(
    suicide_rate ~ `Alcohol consumption (litres)` +
        `Unemployment rate (%)` +
        `Urban population (%)`,
    data = model_data,
    family = Gamma(link = "log")
)

# Univariable models
uni_alco <- glm(suicide_rate ~ `Alcohol consumption (litres)`,
                data = model_data, family = Gamma(link = "log"))
uni_unemploy <- glm(suicide_rate ~ `Unemployment rate (%)`,
                    data = model_data, family = Gamma(link = "log"))
uni_urban <- glm(suicide_rate ~ `Urban population (%)`,
                 data = model_data, family = Gamma(link = "log"))

# ============================================================
# Step 2: Create and Merge gtsummary Tables
# ============================================================

# Define a helper to standardize regression tables
make_tbl <- function(model) {
    tbl_regression(model, exponentiate = TRUE)
}

# Merge univariable models into one stack, then merge with multivariable
tbl_merged <- tbl_merge(
    tbls = list(
        tbl_stack(list(make_tbl(uni_alco), make_tbl(uni_unemploy), make_tbl(uni_urban))),
        make_tbl(gamma_model)
    ),
    tab_spanner = c("**Univariable Analysis**", "**Multivariable Analysis**")
)

# Final Styling: Fix headers and formatting
tbl_final <- tbl_merged %>%
    modify_header(
        update = list(
            estimate_1 ~ "**IRR**",
            estimate_2 ~ "**IRR**",
            label ~ "**Predictor Variable**"
        )
    ) %>%
    modify_caption("**Gamma Regression Models of Global Suicide Rates (2020)**") %>%
    bold_labels() %>%
    italicize_levels()

# Print Table
tbl_final



### Forest PLot
# ============================================================
# FULL WORKING CODE: Gamma GLM with Corrected Table Output
# ============================================================

# Load necessary libraries
library(gtsummary)
library(gt)
library(performance)
library(dplyr)

# --- Step 0: Prepare model data ---
# Ensure suicide_rate is strictly positive for Gamma regression
model_data <- merged_data_clean %>%
    select(country, suicide_rate, alco_consu, unemploy_rate, urban_rate) %>%
    filter(suicide_rate > 0) %>%
    na.omit()

# Rename columns for presentation-quality labels
model_data <- model_data %>%
    rename(
        `Alcohol consumption (litres)` = alco_consu,
        `Unemployment rate (%)` = unemploy_rate,
        `Urban population (%)` = urban_rate
    )

cat("Final sample size:", nrow(model_data), "countries\n")

# ============================================================
# Step 1: Fit Gamma GLM Models
# ============================================================

# Multivariable model
gamma_model <- glm(
    suicide_rate ~ `Alcohol consumption (litres)` +
        `Unemployment rate (%)` +
        `Urban population (%)`,
    data = model_data,
    family = Gamma(link = "log")
)

# Univariable models
uni_alco <- glm(suicide_rate ~ `Alcohol consumption (litres)`,
                data = model_data, family = Gamma(link = "log"))
uni_unemploy <- glm(suicide_rate ~ `Unemployment rate (%)`,
                    data = model_data, family = Gamma(link = "log"))
uni_urban <- glm(suicide_rate ~ `Urban population (%)`,
                 data = model_data, family = Gamma(link = "log"))

# ============================================================
# Step 2: Create and Merge gtsummary Tables
# ============================================================

# Define a helper to standardize regression tables
make_tbl <- function(model) {
    tbl_regression(model, exponentiate = TRUE)
}

# Merge univariable models into one stack, then merge with multivariable
tbl_merged <- tbl_merge(
    tbls = list(
        tbl_stack(list(make_tbl(uni_alco), make_tbl(uni_unemploy), make_tbl(uni_urban))),
        make_tbl(gamma_model)
    ),
    tab_spanner = c("**Univariable Analysis**", "**Multivariable Analysis**")
)

# Final Styling: Fix headers and formatting
tbl_final <- tbl_merged %>%
    modify_header(
        update = list(
            estimate_1 ~ "**IRR**",
            estimate_2 ~ "**IRR**",
            label ~ "**Predictor Variable**"
        )
    ) %>%
    modify_caption("**Gamma Regression Models of Global Suicide Rates (2020)**") %>%
    bold_labels() %>%
    italicize_levels()

# Print Table
tbl_final





library(sjPlot)





plot_model(gamma_model,
           type = "est",
           exponentiate = TRUE,
           show.values = TRUE,
           value.offset = .4,            # Increased slightly for better readability
           vline.color = "grey80",
           title = "Predictors of Global Suicide Rates (2020)",
           axis.labels = c(
               "Urban population (%)",
               "Unemployment rate (%)",
               "Alcohol consumption (litres)"
           ),
           colors = "Set1") +
    theme_minimal() +
    theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        # --- Remove all grid lines ---
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Add an axis line to keep the plot grounded
        axis.line.x = element_line(color = "black"),
        axis.ticks.x = element_line(color = "black")
    )



### Model Diagnostic:

dev_resid <- residuals(gamma_model, type = "deviance")

# Shapiro-Wilk test  (H0: residuals are normally distributed)
shapiro_result <- shapiro.test(dev_resid)
cat("\n--- Shapiro-Wilk Test (Deviance Residuals) ---\n")
cat("W =", round(shapiro_result$statistic, 4),
    " p-value =", round(shapiro_result$p.value, 4), "\n")
cat(ifelse(shapiro_result$p.value > 0.05,
           ">> Do not reject H0: residuals are approximately normal.\n",
           ">> Reject H0: residuals deviate from normality.\n   (Mild deviations are acceptable in GLMs with moderate n.)\n"))

# QQ plot of deviance residuals
qqnorm(dev_resid,
       main = "QQ Plot of Deviance Residuals (Gamma GLM)",
       xlab = "Theoretical Quantiles",
       ylab = "Deviance Residuals",
       pch  = 19, col = "#3498DB", cex = 0.8)
qqline(dev_resid, col = "red", lwd = 2)

# -----------------------------------------------------------------
# 6b. Homoscedasticity (constant variance)
# -----------------------------------------------------------------
# Scale-Location plot: sqrt(|deviance residuals|) vs fitted values.
# A flat loess line indicates constant variance.
fitted_vals   <- fitted(gamma_model)
scaled_resid  <- sqrt(abs(dev_resid))

plot(fitted_vals, scaled_resid,
     main  = "Scale-Location Plot (Homoscedasticity Check)",
     xlab  = "Fitted Values",
     ylab  = "sqrt(|Deviance Residuals|)",
     pch   = 19, col = "#3498DB", cex = 0.7)
lines(loess(scaled_resid ~ fitted_vals, span = 0.5),
      col = "red", lwd = 2)

# Breusch-Pagan test  (H0: homoscedasticity holds)
bp_result <- performance::test_homoscedasticity(gamma_model)
print(bp_result)

# -----------------------------------------------------------------
# 6c. Influential observations (Cook's distance)
# -----------------------------------------------------------------
cook      <- cooks.distance(gamma_model)
n         <- nrow(model_data)
threshold <- 4 / n          # Common rule-of-thumb cutoff

plot(seq_along(cook), cook,
     type = "h",
     main = "Cook's Distance by Observation",
     xlab = "Observation Index",
     ylab = "Cook's Distance",
     col  = ifelse(cook > threshold, "red", "#3498DB"))
abline(h = threshold, col = "red", lty = 2, lwd = 1.5)
text(x = which(cook > threshold),
     y = cook[cook > threshold],
     labels = which(cook > threshold),
     pos = 3, cex = 0.75, col = "red")

cat("\n--- Influential Observations (Cook's D >", round(threshold, 4), ") ---\n")
influential <- which(cook > threshold)
if (length(influential) > 0) {
    cat("Indices  :", influential, "\n")
    cat("Countries:", as.character(model_data$country[influential]), "\n")
} else {
    cat("None detected.\n")
}

# -----------------------------------------------------------------
# 6d. Multicollinearity (VIF)
# -----------------------------------------------------------------
# VIF < 5 = acceptable | 5-10 = moderate | > 10 = high concern
vif_values <- performance::check_collinearity(gamma_model)
print(vif_values)

cat("\n--- VIF Interpretation ---\n")
vif_df <- as.data.frame(vif_values)
for (i in seq_len(nrow(vif_df))) {
    v    <- vif_df$VIF[i]
    flag <- ifelse(v < 5, "OK",
                   ifelse(v < 10, "Moderate concern", "HIGH - consider removing"))
    cat(sprintf("  %-50s VIF = %5.2f  [%s]\n", vif_df$Term[i], v, flag))
}

# -----------------------------------------------------------------
# 6e. Overall diagnostics panel (performance package)
# -----------------------------------------------------------------
check_model(gamma_model)

# ============================================================
# Step 7: Merged Univariable + Multivariable gtsummary Table
# ============================================================

# Custom tidy function: manually exponentiates estimates and CIs.
# This avoids the label auto-detection bug that Gamma GLMs trigger
# when exponentiate = TRUE is used directly.
tidy_gamma_exp <- function(x, ...) {
    broom.helpers::tidy_parameters(x, ...) |>
        dplyr::mutate(
            estimate  = exp(estimate),
            conf.low  = exp(conf.low),
            conf.high = exp(conf.high)
        )
}

# --- Univariable table for each predictor ---
tbl_uni_alco <- tbl_regression(uni_alco,
                               exponentiate = FALSE,
                               show_p_stars = FALSE,
                               tidy_fun = tidy_gamma_exp) |>
    add_n()

tbl_uni_unemploy <- tbl_regression(uni_unemploy,
                                   exponentiate = FALSE,
                                   show_p_stars = FALSE,
                                   tidy_fun = tidy_gamma_exp) |>
    add_n()

tbl_uni_urban <- tbl_regression(uni_urban,
                                exponentiate = FALSE,
                                show_p_stars = FALSE,
                                tidy_fun = tidy_gamma_exp) |>
    add_n()

# --- Multivariable table ---
tbl_multi <- tbl_regression(gamma_model,
                            exponentiate = FALSE,
                            show_p_stars = FALSE,
                            tidy_fun = tidy_gamma_exp) |>
    add_n()

# --- Stack univariable tables, then merge with multivariable ---
tbl_uni_stacked <- tbl_stack(
    list(
        tbl_uni_alco,
        tbl_uni_unemploy,
        tbl_uni_urban
    ),
    quiet = TRUE
)

# Merge: Univariable | Multivariable side by side
tbl_merged <- tbl_merge(
    list(tbl_uni_stacked, tbl_multi),
    tab_spanner = c("Univariable", "Multivariable")
) |>
    modify_columns(estimate_1 ~ label = "Rate Ratio",
                   estimate_2 ~ label = "Rate Ratio") |>
    as_gt() |>
    tab_caption("Gamma Regression: Predictors of Suicide Rates (Rate Ratios with 95% CI)")

# Print the merged table
tbl_merged

# ============================================================
# Step 8: Forest Plot (Multivariable Model)
# ============================================================

tidy(gamma_model, exponentiate = TRUE, conf.int = TRUE) |>
    filter(term != "(Intercept)") |>
    ggplot(aes(x = estimate, y = reorder(term, estimate))) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red", linewidth = 1) +
    geom_point(size = 4, color = "#3498DB") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                   height = 0.2, linewidth = 1, color = "#3498DB") +
    labs(
        title    = "Gamma Regression: Predictors of Suicide Rates",
        subtitle = "Rate Ratios with 95% Confidence Intervals",
        x        = "Rate Ratio (RR)",
        y        = NULL,
        caption  = "Note: RR > 1 indicates positive association; RR < 1 indicates negative association"
    ) +
    scale_x_continuous(breaks = seq(0, 3, 0.5)) +
    theme_minimal(base_size = 13) +
    theme(
        plot.title         = element_text(face = "bold", size = 16),
        plot.subtitle      = element_text(color = "grey40"),
        axis.title.x       = element_text(face = "bold", margin = margin(t = 10)),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank()
    )
