# Multi-Country Excess Mortality Analysis
# LinkedIn-Ready Analysis: UK, Brazil, South Africa
#
# Analysis Options:
#   A) Reported COVID-19 deaths vs Excess deaths
#   B) Excess mortality per 100,000 population
#   C) Baseline (2015-2019) vs Pandemic period (2020-2024)

## 1. ENVIRONMENT SETUP

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
    tidyverse, readr, MASS, scales, patchwork,
    lubridate, viridis, ggrepel, countrycode
)

theme_set(theme_minimal(base_size = 12, base_family = "sans"))

# Professional color palette based on national flags
country_colors <- c(
    "United Kingdom" = "#012169",  # Royal blue from Union Jack
    "Brazil" = "#009B3A",           # Green from Brazilian flag
    "South Africa" = "#FFB81C",     # Gold from South African flag
    "Global Average" = "#d62728"
)

## 2. DATA IMPORT

cat("Loading World Mortality Dataset...\n")
url <- "https://raw.githubusercontent.com/akarlinsky/world_mortality/master/world_mortality.csv"
mort_raw <- read_csv(url, show_col_types = FALSE)

cat("Dataset loaded successfully.\n\n")

## 3. COUNTRY-SPECIFIC ANALYSIS FUNCTION

analyze_country <- function(data, country_name) {

    cat(paste0("Analyzing: ", country_name, "\n"))

    # Check what time_unit is available for this country
    available_units <- data %>%
        filter(country_name == !!country_name) %>%
        distinct(time_unit) %>%
        pull(time_unit)

    cat(paste0("  Available time units: ", paste(available_units, collapse = ", "), "\n"))

    # Prefer weekly, but use monthly if that's all that's available
    use_unit <- if("weekly" %in% available_units) "weekly" else "monthly"
    periods_per_year <- if(use_unit == "weekly") 52 else 12

    cat(paste0("  Using: ", use_unit, " data\n"))

    # Baseline period (2015-2019)
    baseline <- data %>%
        filter(
            country_name == !!country_name,
            time_unit == use_unit,
            year >= 2015,
            year <= 2019
        ) %>%
        mutate(
            t = year + time / periods_per_year,
            sin_annual = sin(2 * pi * time / periods_per_year),
            cos_annual = cos(2 * pi * time / periods_per_year),
            sin_semi = sin(4 * pi * time / periods_per_year),
            cos_semi = cos(4 * pi * time / periods_per_year)
        )

    if (nrow(baseline) == 0) {
        cat(paste0("  WARNING: No baseline data for ", country_name, "\n"))
        return(NULL)
    }

    # Fit model
    model <- tryCatch({
        glm.nb(deaths ~ t + sin_annual + cos_annual + sin_semi + cos_semi, data = baseline)
    }, error = function(e) {
        # Try simpler model if enhanced fails
        glm.nb(deaths ~ t + sin_annual + cos_annual, data = baseline)
    })

    # Full dataset
    full_data <- data %>%
        filter(
            country_name == !!country_name,
            time_unit == use_unit
        ) %>%
        mutate(
            t = year + time / periods_per_year,
            sin_annual = sin(2 * pi * time / periods_per_year),
            cos_annual = cos(2 * pi * time / periods_per_year),
            sin_semi = sin(4 * pi * time / periods_per_year),
            cos_semi = cos(4 * pi * time / periods_per_year),
            period = time,
            time_unit_used = use_unit
        ) %>%
        arrange(year, time)

    # Predictions
    pred <- predict(model, newdata = full_data, type = "link", se.fit = TRUE)

    full_data <- full_data %>%
        mutate(
            expected = exp(pred$fit),
            lower_95 = exp(pred$fit - 1.96 * pred$se.fit),
            upper_95 = exp(pred$fit + 1.96 * pred$se.fit),
            excess = deaths - expected,
            excess_lower = deaths - upper_95,
            excess_upper = deaths - lower_95,
            pct_excess = 100 * excess / expected,
            significant = deaths > upper_95 | deaths < lower_95,
            country = country_name
        )

    cat(paste0("  ✓ Complete: ", nrow(full_data), " ", use_unit, " periods analyzed\n"))

    return(list(
        data = full_data,
        model = model,
        baseline_deaths = mean(baseline$deaths, na.rm = TRUE),
        time_unit = use_unit
    ))
}

## 4. ANALYZE ALL COUNTRIES

# Check for Brazil variations in the dataset
cat("Checking for Brazil in dataset...\n")
brazil_check <- mort_raw %>%
    filter(grepl("brazil|Brasil", country_name, ignore.case = TRUE)) %>%
    distinct(country_name)

if (nrow(brazil_check) > 0) {
    cat("Found Brazil as:", brazil_check$country_name[1], "\n")
    brazil_name <- brazil_check$country_name[1]
} else {
    cat("WARNING: Brazil not found in dataset. Using 'Brazil' anyway.\n")
    brazil_name <- "Brazil"
}

countries <- c("United Kingdom", brazil_name, "South Africa")

results <- map(countries, ~analyze_country(mort_raw, .x))
names(results) <- c("United Kingdom", "Brazil", "South Africa")

# Extract data and standardize country names
all_data <- map_df(results, ~.x$data) %>%
    mutate(country = case_when(
        country == brazil_name ~ "Brazil",
        TRUE ~ country
    ))

# Get population data (approximate, for per capita calculations)
# Note: These are approximate 2020 populations in millions
population_data <- tibble(
    country = c("United Kingdom", "Brazil", "South Africa"),
    population = c(67.1, 212.6, 59.3)  # in millions
)

## 5. SUMMARY STATISTICS

cat("\n### COUNTRY COMPARISONS ###\n\n")

# Overall summary (2020+)
country_summary <- all_data %>%
    filter(year >= 2020) %>%
    group_by(country, time_unit_used) %>%
    summarise(
        `Total Deaths` = sum(deaths, na.rm = TRUE),
        `Expected Deaths` = round(sum(expected, na.rm = TRUE)),
        `Excess Deaths` = round(sum(excess, na.rm = TRUE)),
        `P-score (%)` = round(100 * sum(excess, na.rm = TRUE) /
                                  sum(expected, na.rm = TRUE), 1),
        `Periods Analyzed` = n(),
        `Periods Significant` = sum(significant, na.rm = TRUE),
        `% Periods Significant` = round(100 * mean(significant, na.rm = TRUE), 1),
        .groups = "drop"
    ) %>%
    left_join(population_data, by = "country") %>%
    mutate(
        `Excess per 100k` = round((`Excess Deaths` / population) / 10, 1)
    )

print(country_summary)
cat("\n")

# Annual comparison
annual_comparison <- all_data %>%
    filter(year >= 2020) %>%
    group_by(country, year) %>%
    summarise(
        excess = sum(excess, na.rm = TRUE),
        pct_excess = 100 * sum(excess, na.rm = TRUE) / sum(expected, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    left_join(population_data, by = "country") %>%
    mutate(
        excess_per_100k = (excess / population) / 10
    )

cat("### ANNUAL COMPARISON ###\n")
print(annual_comparison %>% dplyr::select(country, year, excess, pct_excess, excess_per_100k))
cat("\n")

## 6. STATISTICAL COMPARISON BETWEEN COUNTRIES

cat("### STATISTICAL SIGNIFICANCE TESTS ###\n\n")

# Prepare data for statistical testing
comparison_data <- all_data %>%
    filter(year >= 2020) %>%
    dplyr::select(country, year, time, excess, pct_excess, expected, deaths, significant)

# Test 1: Overall excess mortality differences (Kruskal-Wallis test)
cat("1. OVERALL EXCESS MORTALITY COMPARISON\n")
cat("   Test: Kruskal-Wallis (non-parametric)\n")
kw_test <- kruskal.test(excess ~ country, data = comparison_data)
cat(sprintf("   Chi-squared = %.2f, p-value = %.4f\n", kw_test$statistic, kw_test$p.value))
if (kw_test$p.value < 0.05) {
    cat("   ✓ Significant differences exist between countries (p < 0.05)\n\n")
} else {
    cat("   ✗ No significant differences between countries (p ≥ 0.05)\n\n")
}

# Test 2: Pairwise comparisons (Wilcoxon test with Bonferroni correction)
cat("2. PAIRWISE COMPARISONS (Wilcoxon Rank-Sum Tests)\n")
cat("   Bonferroni-adjusted p-values:\n\n")

countries_list <- unique(comparison_data$country)
pairwise_results <- data.frame(
    Comparison = character(),
    W_statistic = numeric(),
    p_value = numeric(),
    adj_p_value = numeric(),
    Significant = character(),
    stringsAsFactors = FALSE
)

n_comparisons <- choose(length(countries_list), 2)
comparison_idx <- 1

for (i in 1:(length(countries_list)-1)) {
    for (j in (i+1):length(countries_list)) {
        c1 <- countries_list[i]
        c2 <- countries_list[j]

        data1 <- comparison_data %>% filter(country == c1) %>% pull(excess)
        data2 <- comparison_data %>% filter(country == c2) %>% pull(excess)

        wtest <- wilcox.test(data1, data2, exact = FALSE)
        adj_p <- min(wtest$p.value * n_comparisons, 1)  # Bonferroni correction

        pairwise_results <- rbind(pairwise_results, data.frame(
            Comparison = paste(c1, "vs", c2),
            W_statistic = wtest$statistic,
            p_value = wtest$p.value,
            adj_p_value = adj_p,
            Significant = ifelse(adj_p < 0.05, "Yes ***",
                                 ifelse(adj_p < 0.10, "Yes *", "No"))
        ))
    }
}

print(pairwise_results, row.names = FALSE)
cat("\n")

# Test 3: Summary statistics table with confidence intervals
cat("3. SUMMARY STATISTICS WITH 95% CONFIDENCE INTERVALS\n\n")

summary_stats <- comparison_data %>%
    group_by(country) %>%
    summarise(
        `Mean Excess` = mean(excess, na.rm = TRUE),
        `Median Excess` = median(excess, na.rm = TRUE),
        `SD` = sd(excess, na.rm = TRUE),
        `SE` = sd(excess, na.rm = TRUE) / sqrt(n()),
        `CI Lower` = `Mean Excess` - 1.96 * `SE`,
        `CI Upper` = `Mean Excess` + 1.96 * `SE`,
        `N Periods` = n(),
        .groups = "drop"
    ) %>%
    mutate(across(where(is.numeric), ~round(., 1)))

print(summary_stats)
cat("\n")

# Test 4: Effect sizes (Cohen's d for practical significance)
cat("4. EFFECT SIZES (Cohen's d - Practical Significance)\n")
cat("   Interpretation: Small (0.2), Medium (0.5), Large (0.8)\n\n")

effect_sizes <- data.frame(
    Comparison = character(),
    Cohens_d = numeric(),
    Interpretation = character(),
    stringsAsFactors = FALSE
)

for (i in 1:(length(countries_list)-1)) {
    for (j in (i+1):length(countries_list)) {
        c1 <- countries_list[i]
        c2 <- countries_list[j]

        data1 <- comparison_data %>% filter(country == c1) %>% pull(excess)
        data2 <- comparison_data %>% filter(country == c2) %>% pull(excess)

        # Calculate Cohen's d
        pooled_sd <- sqrt(((length(data1)-1)*sd(data1)^2 + (length(data2)-1)*sd(data2)^2) /
                              (length(data1) + length(data2) - 2))
        cohens_d <- (mean(data1) - mean(data2)) / pooled_sd

        interpretation <- case_when(
            abs(cohens_d) < 0.2 ~ "Negligible",
            abs(cohens_d) < 0.5 ~ "Small",
            abs(cohens_d) < 0.8 ~ "Medium",
            TRUE ~ "Large"
        )

        effect_sizes <- rbind(effect_sizes, data.frame(
            Comparison = paste(c1, "vs", c2),
            Cohens_d = round(cohens_d, 2),
            Interpretation = interpretation
        ))
    }
}

print(effect_sizes, row.names = FALSE)
cat("\n")

# Test 5: Proportion of significant excess periods comparison
cat("5. PROPORTION OF SIGNIFICANT EXCESS PERIODS\n\n")

prop_comparison <- comparison_data %>%
    group_by(country) %>%
    summarise(
        `Total Periods` = n(),
        `Significant Periods` = sum(significant, na.rm = TRUE),
        `Proportion` = mean(significant, na.rm = TRUE),
        `Percentage` = round(100 * mean(significant, na.rm = TRUE), 1),
        .groups = "drop"
    )

print(prop_comparison)
cat("\n")

# Chi-square test for proportion differences
prop_table <- prop_comparison %>%
    dplyr::select(country, `Significant Periods`, `Total Periods`) %>%
    mutate(`Non-Significant` = `Total Periods` - `Significant Periods`) %>%
    dplyr::select(country, `Significant Periods`, `Non-Significant`)

# Create contingency table
contingency <- as.matrix(prop_table[,-1])
rownames(contingency) <- prop_table$country

chi_test <- chisq.test(contingency)
cat(sprintf("Chi-square test for proportion differences:\n"))
cat(sprintf("   Chi-squared = %.2f, p-value = %.4f\n", chi_test$statistic, chi_test$p.value))
if (chi_test$p.value < 0.05) {
    cat("   ✓ Significant difference in proportions between countries (p < 0.05)\n\n")
} else {
    cat("   ✗ No significant difference in proportions (p ≥ 0.05)\n\n")
}

cat("##############################\n\n")

## 7. OPTION A: COVID-19 DEATHS vs EXCESS DEATHS

cat("Preparing Option A visualizations...\n")

# Cumulative excess by country
cumulative_data <- all_data %>%
    filter(year >= 2020) %>%
    arrange(country, year, time) %>%
    group_by(country) %>%
    mutate(
        cumulative_excess = cumsum(excess)
    ) %>%
    ungroup()

# Plot A1: Cumulative excess deaths over time
p_a1 <- ggplot(cumulative_data, aes(x = t, y = cumulative_excess,
                                    color = country, fill = country)) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = country_colors, name = NULL) +
    scale_y_continuous(labels = comma_format(), breaks = pretty_breaks(n = 8)) +
    labs(
        title = "Cumulative Excess Deaths: UK vs Brazil vs South Africa",
        subtitle = "Total deaths above expected baseline (2015-2019) since 2020",
        x = "Year",
        y = "Cumulative excess deaths",
        caption = "Data: World Mortality Dataset | Note: Countries have different population sizes"
    ) +
    theme(
        legend.position = "top",
        plot.title = element_text(face = "bold", size = 14),
        panel.grid.minor = element_blank()
    )

# Plot A2: Weekly excess deaths (faceted)
p_a2 <- ggplot(all_data %>% filter(year >= 2020),
               aes(x = t, y = excess, fill = significant)) +
    geom_col(width = 0.019, alpha = 0.8) +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    facet_wrap(~country, scales = "free_y", ncol = 1) +
    scale_fill_manual(
        values = c("TRUE" = "#C0392B", "FALSE" = "#7F8C8D"),
        labels = c("Within expected", "Significant excess/deficit"),
        name = NULL
    ) +
    scale_y_continuous(labels = comma_format()) +
    labs(
        title = "Excess Mortality by Country",
        subtitle = "Red bars indicate statistically significant excess (95% CI) | Note: Brazil uses monthly data",
        x = "Year",
        y = "Excess deaths",
        caption = "Negative values = fewer deaths than expected"
    ) +
    theme(
        legend.position = "top",
        plot.title = element_text(face = "bold", size = 14),
        strip.text = element_text(face = "bold", size = 12),
        panel.grid.minor = element_blank()
    )

## 8. OPTION B: EXCESS PER 100,000 POPULATION

cat("Preparing Option B visualizations...\n")

# Calculate excess per 100k for all available countries
all_countries_summary <- mort_raw %>%
    filter(
        time_unit == "weekly",
        year >= 2020,
        !is.na(deaths)
    ) %>%
    group_by(country_name) %>%
    summarise(
        total_deaths = sum(deaths, na.rm = TRUE),
        weeks = n(),
        .groups = "drop"
    ) %>%
    filter(weeks >= 100)  # Only countries with substantial data

# For our focus countries with calculated excess
per_capita_data <- country_summary %>%
    dplyr::select(country, `Excess Deaths`, `Excess per 100k`, `P-score (%)`) %>%
    arrange(desc(`Excess per 100k`))

# Plot B1: Excess per 100k ranking
p_b1 <- ggplot(per_capita_data,
               aes(x = reorder(country, `Excess per 100k`),
                   y = `Excess per 100k`, fill = country)) +
    geom_col(show.legend = FALSE) +
    geom_text(aes(label = scales::comma(`Excess per 100k`, accuracy = 1)),
              hjust = -0.2, fontface = "bold", size = 4) +
    coord_flip() +
    scale_fill_manual(values = country_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
        title = "Excess Mortality per 100,000 Population (2020-2024)",
        subtitle = "Standardized comparison accounting for population size",
        x = NULL,
        y = "Excess deaths per 100,000 population",
        caption = "Data: World Mortality Dataset | Higher values = greater mortality impact"
    ) +
    theme(
        plot.title = element_text(face = "bold", size = 14),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 12, face = "bold")
    )

# Plot B2: P-score comparison
p_b2 <- ggplot(per_capita_data,
               aes(x = reorder(country, `P-score (%)`),
                   y = `P-score (%)`, fill = country)) +
    geom_col(show.legend = FALSE) +
    geom_text(aes(label = paste0(sprintf("%.1f", `P-score (%)`), "%")),
              hjust = -0.2, fontface = "bold", size = 4) +
    coord_flip() +
    scale_fill_manual(values = country_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
        title = "Percentage Excess Mortality (P-score)",
        subtitle = "% by which observed deaths exceeded expected (2020-2024)",
        x = NULL,
        y = "P-score (%)",
        caption = "P-score = (Excess Deaths / Expected Deaths) × 100"
    ) +
    theme(
        plot.title = element_text(face = "bold", size = 14),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 12, face = "bold")
    )

## 9. OPTION C: BASELINE vs PANDEMIC COMPARISON

cat("Preparing Option C visualizations...\n")

# Calculate baseline vs pandemic averages
baseline_vs_pandemic <- all_data %>%
    mutate(
        period = if_else(year < 2020, "Baseline (2015-2019)",
                         "Pandemic (2020-2024)")
    ) %>%
    group_by(country, period) %>%
    summarise(
        mean_period_deaths = mean(deaths, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    pivot_wider(names_from = period, values_from = mean_period_deaths) %>%
    mutate(
        increase = `Pandemic (2020-2024)` - `Baseline (2015-2019)`,
        pct_increase = 100 * increase / `Baseline (2015-2019)`
    )

cat("### BASELINE vs PANDEMIC ###\n")
print(baseline_vs_pandemic)
cat("\n")

# Plot C1: Baseline vs Pandemic side-by-side
period_comparison <- all_data %>%
    mutate(
        period = if_else(year < 2020, "Baseline\n(2015-2019)",
                         "Pandemic\n(2020-2024)")
    ) %>%
    group_by(country, period) %>%
    summarise(
        mean_deaths = mean(deaths, na.rm = TRUE),
        .groups = "drop"
    )

p_c1 <- ggplot(period_comparison,
               aes(x = country, y = mean_deaths, fill = period)) +
    geom_col(position = "dodge", width = 0.7) +
    geom_text(aes(label = scales::comma(round(mean_deaths))),
              position = position_dodge(width = 0.7),
              vjust = -0.5, fontface = "bold", size = 3.5) +
    scale_fill_manual(
        values = c("Baseline\n(2015-2019)" = "#95A5A6",
                   "Pandemic\n(2020-2024)" = "#E74C3C"),
        name = "Period"
    ) +
    scale_y_continuous(labels = comma_format(), expand = expansion(mult = c(0, 0.15))) +
    labs(
        title = "Mean Deaths per Period: Baseline vs Pandemic",
        subtitle = "Comparing pre-pandemic (2015-2019) to pandemic era (2020-2024) | Note: Brazil uses monthly data",
        x = NULL,
        y = "Mean deaths per period",
        caption = "Data: World Mortality Dataset | Shows average mortality per time period"
    ) +
    theme(
        legend.position = "top",
        plot.title = element_text(face = "bold", size = 14),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold")
    )

# Plot C2: Time series showing baseline and pandemic periods
p_c2 <- ggplot(all_data, aes(x = t, y = deaths, color = country)) +
    geom_line(alpha = 0.7, linewidth = 0.6) +
    geom_vline(xintercept = 2020, linetype = "dashed",
               color = "red", linewidth = 1) +
    annotate("rect", xmin = 2015, xmax = 2019.99,
             ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "blue") +
    annotate("rect", xmin = 2020, xmax = Inf,
             ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "red") +
    annotate("text", x = 2017.5, y = Inf, label = "BASELINE",
             vjust = 1.5, fontface = "bold", color = "blue") +
    annotate("text", x = 2022, y = Inf, label = "PANDEMIC",
             vjust = 1.5, fontface = "bold", color = "red") +
    facet_wrap(~country, scales = "free_y", ncol = 1) +
    scale_color_manual(values = country_colors, name = NULL) +
    scale_y_continuous(labels = comma_format()) +
    labs(
        title = "Deaths Over Time: Baseline vs Pandemic Periods",
        subtitle = "Blue shaded = baseline (2015-2019) | Red shaded = pandemic (2020-2024) | Note: Brazil uses monthly data",
        x = "Year",
        y = "Deaths per period",
        caption = "Data: World Mortality Dataset"
    ) +
    theme(
        legend.position = "none",
        plot.title = element_text(face = "bold", size = 14),
        strip.text = element_text(face = "bold", size = 12),
        panel.grid.minor = element_blank()
    )

## 10. COMBINED INSIGHT PLOT (LinkedIn Hero Image)

cat("Creating LinkedIn hero visualization...\n")

# Create a comprehensive 2x2 panel
hero_data <- annual_comparison %>%
    filter(year %in% c(2020, 2021, 2022, 2023))

# Panel 1: Annual excess per 100k
p_hero1 <- ggplot(hero_data,
                  aes(x = factor(year), y = excess_per_100k,
                      fill = country, group = country)) +
    geom_col(position = "dodge", width = 0.7) +
    scale_fill_manual(values = country_colors, name = NULL) +
    labs(
        title = "A) Excess Deaths per 100k by Year",
        x = "Year", y = "Per 100,000 pop."
    ) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", size = 11))

# Panel 2: P-score trends
p_hero2 <- ggplot(hero_data,
                  aes(x = year, y = pct_excess, color = country, group = country)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    scale_color_manual(values = country_colors, name = NULL) +
    labs(
        title = "B) P-score Trends Over Time",
        x = "Year", y = "P-score (%)"
    ) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", size = 11))

# Panel 3: Cumulative comparison
p_hero3 <- ggplot(cumulative_data %>% filter(year >= 2020),
                  aes(x = t, y = cumulative_excess, color = country)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = country_colors, name = NULL) +
    scale_y_continuous(labels = comma_format()) +
    labs(
        title = "C) Cumulative Excess Deaths",
        x = "Year", y = "Cumulative excess"
    ) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", size = 11))

# Panel 4: Total comparison
total_comparison <- country_summary %>%
    dplyr::select(country, `Excess Deaths`) %>%
    mutate(country_short = case_when(
        country == "United Kingdom" ~ "UK",
        country == "South Africa" ~ "South Africa",
        TRUE ~ country
    ))

p_hero4 <- ggplot(total_comparison,
                  aes(x = reorder(country_short, `Excess Deaths`),
                      y = `Excess Deaths`, fill = country)) +
    geom_col() +
    geom_text(aes(label = scales::comma(`Excess Deaths`)),
              hjust = -0.1, fontface = "bold", size = 3) +
    coord_flip() +
    scale_fill_manual(values = country_colors) +
    scale_y_continuous(labels = comma_format(), expand = expansion(mult = c(0, 0.15))) +
    labs(
        title = "D) Total Excess Deaths (2020-2024)",
        x = NULL, y = "Total excess deaths"
    ) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", size = 11),
          panel.grid.major.y = element_blank())

# Combine into hero plot
p_hero <- (p_hero1 | p_hero2) / (p_hero3 | p_hero4) +
    plot_annotation(
        title = "EXCESS MORTALITY COMPARISON: UK vs Brazil vs South Africa (2020-2024)",
        subtitle = "Multi-dimensional analysis of pandemic impact across three continents",
        caption = "Data: World Mortality Dataset (Karlinsky & Kobak) | Baseline: 2015-2019",
        theme = theme(
            plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
            plot.subtitle = element_text(size = 12, hjust = 0.5),
            plot.caption = element_text(size = 9, hjust = 1)
        )
    )

# Add legend separately
legend_data <- tibble(
    country = names(country_colors),
    y = 1:length(country_colors)
)

## 11. DISPLAY ALL PLOTS

cat("\n")
cat("##############################\n")
cat("DISPLAYING ALL VISUALIZATIONS\n")
cat("##############################\n\n")

cat("OPTION A: COVID-19 vs Excess Deaths\n")
print(p_a1)
print(p_a2)

cat("\nOPTION B: Per Capita Analysis\n")
print(p_b1)
print(p_b2)

cat("\nOPTION C: Baseline vs Pandemic\n")
print(p_c1)
print(p_c2)

cat("\nLINKEDIN HERO IMAGE\n")
print(p_hero)

## 12. KEY INSIGHTS FOR LINKEDIN

cat("\n")
cat("##############################\n")
cat("KEY INSIGHTS FOR LINKEDIN POST\n")
cat("##############################\n\n")

# Calculate key statistics
uk_stats <- country_summary %>% filter(country == "United Kingdom")
brazil_stats <- country_summary %>% filter(country == "Brazil")
sa_stats <- country_summary %>% filter(country == "South Africa")

cat("💡 KEY FINDINGS:\n\n")

cat("1. ABSOLUTE IMPACT (Total Excess Deaths):\n")
cat(sprintf("   • Brazil: %s excess deaths\n",
            scales::comma(brazil_stats$`Excess Deaths`)))
cat(sprintf("   • UK: %s excess deaths\n",
            scales::comma(uk_stats$`Excess Deaths`)))
cat(sprintf("   • South Africa: %s excess deaths\n\n",
            scales::comma(sa_stats$`Excess Deaths`)))

cat("2. RELATIVE IMPACT (Per 100,000 Population):\n")
cat(sprintf("   • South Africa: %.1f per 100k (Highest)\n",
            sa_stats$`Excess per 100k`))
cat(sprintf("   • UK: %.1f per 100k\n",
            uk_stats$`Excess per 100k`))
cat(sprintf("   • Brazil: %.1f per 100k\n\n",
            brazil_stats$`Excess per 100k`))

cat("3. P-SCORE (% Excess):\n")
cat(sprintf("   • Brazil: %.1f%% above expected\n",
            brazil_stats$`P-score (%)`))
cat(sprintf("   • UK: %.1f%% above expected\n",
            uk_stats$`P-score (%)`))
cat(sprintf("   • South Africa: %.1f%% above expected\n\n",
            sa_stats$`P-score (%)`))

cat("4. STATISTICAL SIGNIFICANCE:\n")
cat(sprintf("   • UK: %.1f%% of periods showed significant excess\n",
            uk_stats$`% Periods Significant`))
cat(sprintf("   • Brazil: %.1f%% of periods showed significant excess\n",
            brazil_stats$`% Periods Significant`))
cat(sprintf("   • South Africa: %.1f%% of periods showed significant excess\n\n",
            sa_stats$`% Periods Significant`))

cat("##############################\n")
cat("LINKEDIN POST SUGGESTIONS\n")
cat("##############################\n\n")

cat("📊 OPTION 1 - Surprising Finding:\n")
cat('"Why excess mortality tells a different story than reported COVID deaths"\n')
cat("👉 Hook: Population size matters. Here's what the data really shows...\n\n")

cat("📊 OPTION 2 - Data-Driven Insight:\n")
cat('"I analyzed excess mortality across 3 continents. The results surprised me."\n')
cat("👉 Use the hero image with key stats in the post text\n\n")

cat("📊 OPTION 3 - Educational:\n")
cat('"What is excess mortality and why does it matter in public health?"\n')
cat("👉 Explain baseline concept with real data from UK, Brazil, SA\n\n")

cat("##############################\n\n")

cat("ANALYSIS COMPLETE!\n\n")
cat("Available plot objects:\n")
cat("  p_a1, p_a2 - Option A plots\n")
cat("  p_b1, p_b2 - Option B plots\n")
cat("  p_c1, p_c2 - Option C plots\n")
cat("  p_hero - LinkedIn hero image (4-panel)\n")
cat("\nData objects:\n")
cat("  all_data - Combined country data\n")
cat("  country_summary - Overall statistics\n")
cat("  annual_comparison - Year-by-year comparison\n")
cat("##############################\n")