# =============================================================================
# GLOBAL NEONATAL MORTALITY — SPATIAL CLUSTER ANALYSIS
# =============================================================================
# Author  : Kabinga Chanda
# Date    : May 2026
# Data    : UNICEF Neonatal Mortality Estimates, 2024
# Methods : Global Moran's I + Local Moran's I (LISA)
# =============================================================================


# 1. LOAD LIBRARIES -----------------------------------------------------------

pacman::p_load(
    tidyverse,   # data manipulation and ggplot2
    rio,         # flexible data import
    spdep,       # spatial weights and Moran's I
    sf,          # spatial data handling
    gt,          # professional tables
    here         # relative file paths
)


# 2. LOAD DATA ----------------------------------------------------------------

# Neonatal mortality data
data <- import(here("Neonatal_Mortality_Rates_2025-1.xlsx"))

# World shapefile from spData (built-in, clean geometry)
data("world", package = "spData")

# Rename country column to match world shapefile
data <- data |>
    rename(country = Country.Name)

# Check column names and first few rows
names(data)
head(data, 5)


# 3. JOIN DATA TO SHAPEFILE ---------------------------------------------------

world_joint <- world |>
    left_join(
        data |> select(country, `2024.5`),
        by = join_by(name_long == country)
    )

# Check join quality — how many countries matched?
cat(
    "Total countries in shapefile :", nrow(world),          "\n",
    "Countries with mortality data:", sum(!is.na(world_joint$`2024.5`)), "\n",
    "Countries missing data        :", sum(is.na(world_joint$`2024.5`)),  "\n"
)

# Check for unmatched countries — fix name mismatches if needed
unmatched <- data |>
    filter(!country %in% world$name_long) |>
    pull(country)

cat("\nUnmatched country names in data:\n")
print(unmatched)


# 4. RAW CHOROPLETH MAP -------------------------------------------------------

map_raw <- world_joint |>
    ggplot() +
    geom_sf(
        aes(fill = `2024.5`),
        color     = "white",
        linewidth = 0.1
    ) +
    scale_fill_stepsn(
        colours = c(
            "#F7F4F9", "#C994C7", "#DF65B0",
            "#E7298A", "#CE1256", "#91003F"
        ),
        breaks   = c(0, 5, 10, 15, 20, 30, 40),
        limits   = c(0, 50),
        na.value = "#D9D9D9",
        name     = "Deaths per 1,000 live births",
        guide    = guide_colorsteps(
            barwidth       = unit(8, "cm"),
            barheight      = unit(0.4, "cm"),
            title.position = "top",
            title.hjust    = 0.5
        )
    ) +
    labs(
        title   = "Global Neonatal Mortality Rate, 2024",
        caption = "Source: UNICEF estimates, 2024. Grey = no data."
    ) +
    theme_void(base_family = "sans") +
    theme(
        plot.title      = element_text(size = 11, face = "bold",
                                       hjust = 0.5,
                                       margin = margin(b = 6)),
        plot.caption    = element_text(size = 7, color = "grey50",
                                       hjust = 0,
                                       margin = margin(t = 8)),
        legend.position = "bottom",
        legend.title    = element_text(size = 8, face = "bold"),
        legend.text     = element_text(size = 7),
        plot.margin     = margin(10, 10, 10, 10),
        plot.background = element_rect(fill = "white", color = NA)
    )

map_raw

ggsave(
    filename = here("outputs", "figures", "01_neonatal_mortality_raw.png"),
    plot     = map_raw,
    width    = 12,
    height   = 7,
    dpi      = 300,
    bg       = "white"
)


# 5. PREPARE DATA FOR SPATIAL ANALYSIS ----------------------------------------

# Remove countries with missing mortality data
# Moran's I cannot handle NA values
world_clean <- world_joint |>
    filter(!is.na(`2024.5`))

cat("Countries in spatial analysis:", nrow(world_clean), "\n")

# Verify geometry is valid — must be TRUE before building weights
cat("All geometries valid:", all(st_is_valid(world_clean)), "\n")

# If FALSE, fix with:
# world_clean <- st_make_valid(world_clean)


# 6. BUILD SPATIAL WEIGHTS MATRIX ---------------------------------------------

# Queen contiguity: countries sharing any border point are neighbours
# style = "W": row-standardised — each neighbour gets equal weight summing to 1
# zero.policy = TRUE: allows island nations with zero neighbours

nb    <- poly2nb(world_clean, queen = TRUE)
listw <- nb2listw(nb, style = "W", zero.policy = TRUE)

# Inspect neighbour structure
summary(nb)

# Identify island nations excluded from analysis
islands <- which(card(nb) == 0)
cat("\nCountries with zero neighbours (islands):\n")
print(world_clean$name_long[islands])


# 7. GLOBAL MORAN'S I ---------------------------------------------------------

moran_result <- moran.test(
    x             = world_clean$`2024.5`,
    listw         = listw,
    randomisation = TRUE,   # permutation test — appropriate for non-normal data
    zero.policy   = TRUE
)

print(moran_result)


# 8. GLOBAL MORAN'S I — RESULTS TABLE ----------------------------------------

moran_table <- tibble(
    Statistic = c(
        "Moran's I statistic",
        "Expected I (H\u2080)",
        "Variance",
        "Standard deviate (z)",
        "p-value",
        "Interpretation"
    ),
    Value = c(
        round(moran_result$estimate["Moran I statistic"], 4),
        round(moran_result$estimate["Expectation"],       4),
        round(moran_result$estimate["Variance"],          6),
        round(moran_result$statistic,                     4),
        "< 0.001",
        "Significant positive spatial autocorrelation"
    )
)

moran_gt <- moran_table |>
    gt() |>
    tab_header(
        title    = md("**Global Moran's I Test**"),
        subtitle = md("*Neonatal Mortality Rate \u2014 Worldwide, 2024*")
    ) |>
    tab_source_note(
        source_note = md(
            "Weights: Queen contiguity, row-standardised (W).
             Test: randomisation. Island nations excluded (zero neighbours)."
        )
    ) |>
    tab_style(
        style     = list(cell_fill(color = "#1F3864"),
                         cell_text(color = "white", weight = "bold")),
        locations = cells_column_labels()
    ) |>
    tab_style(
        style     = cell_fill(color = "#D6E4F0"),
        locations = cells_body(rows = seq(1, nrow(moran_table), 2))
    ) |>
    tab_style(
        style     = cell_text(weight = "bold", color = "#1F3864"),
        locations = cells_body(columns = Statistic)
    ) |>
    tab_style(
        style = list(cell_fill(color = "#D5F5E3"),
                     cell_text(color = "#1E5C2E", weight = "bold")),
        locations = cells_body(rows = Statistic == "Interpretation")
    ) |>
    tab_style(
        style = list(cell_fill(color = "#D5F5E3"),
                     cell_text(color = "#1E5C2E", weight = "bold")),
        locations = cells_body(columns = Value,
                               rows    = Statistic == "Interpretation")
    ) |>
    cols_width(Statistic ~ px(260), Value ~ px(340)) |>
    cols_align(align = "left",  columns = Statistic) |>
    cols_align(align = "right", columns = Value) |>
    opt_table_font(font = google_font("Source Sans Pro")) |>
    opt_stylize(style = 1, color = "blue") |>
    tab_options(
        table.width                = px(600),
        heading.title.font.size    = px(16),
        heading.subtitle.font.size = px(12),
        table.border.top.color     = "#1F3864",
        table.border.top.width     = px(3),
        source_notes.font.size     = px(10)
    )

moran_gt

# Save table
gtsave(moran_gt,
       filename = here("outputs", "tables", "morans_I_results.png"))


# 9. LOCAL MORAN'S I (LISA) ---------------------------------------------------

# Run LISA — returns one row per country
lisa <- localmoran(
    x           = world_clean$`2024.5`,
    listw       = listw,
    zero.policy = TRUE,
    na.action   = na.exclude
)

# Compute spatial lag:
# weighted average of each country's neighbours' mortality rates
world_clean$lag_neo_mort_rat <- lag.listw(
    listw,
    world_clean$`2024.5`,
    zero.policy = TRUE
)

# Mean — the dividing line between high and low for classification
mean_mort_rate   <- mean(world_clean$`2024.5`, na.rm = TRUE)
median_mort_rate <- median(world_clean$`2024.5`, na.rm = TRUE)

cat(
    "Mean mortality rate  :", round(mean_mort_rate,   2), "\n",
    "Median mortality rate:", round(median_mort_rate, 2), "\n"
)

# Attach LISA statistics and classify each country
world_clean <- world_clean |>
    mutate(
        lisa_i = lisa[, "Ii"],
        lisa_p = lisa[, "Pr(z != E(Ii))"],
        
        # Classification:
        # High-High — country is high AND neighbours are high → hotspot
        # Low-Low   — country is low  AND neighbours are low  → cold spot
        # High-Low  — country is high but neighbours are low  → spatial outlier
        # Low-High  — country is low  but neighbours are high → possible gap
        hotspot_type = case_when(
            lisa_p < 0.05 & `2024.5` > mean_mort_rate & lag_neo_mort_rat > mean_mort_rate ~ "High-High",
            lisa_p < 0.05 & `2024.5` < mean_mort_rate & lag_neo_mort_rat < mean_mort_rate ~ "Low-Low",
            lisa_p < 0.05 & `2024.5` > mean_mort_rate & lag_neo_mort_rat < mean_mort_rate ~ "High-Low",
            lisa_p < 0.05 & `2024.5` < mean_mort_rate & lag_neo_mort_rat > mean_mort_rate ~ "Low-High",
            TRUE                                                                           ~ "Not significant"
        ),
        
        # Fix factor order for consistent legend display
        hotspot_type = factor(
            hotspot_type,
            levels = c("High-High", "Low-Low",
                       "High-Low",  "Low-High",
                       "Not significant")
        )
    )

# Cluster count summary
cat("\nLISA cluster counts:\n")
world_clean |>
    st_drop_geometry() |>
    count(hotspot_type, sort = FALSE) |>
    print()


# 10. LISA MAP ----------------------------------------------------------------

lisa_pal <- c(
    "High-High"       = "#E74C3C",   # red    — hotspot
    "Low-Low"         = "#3498DB",   # blue   — cold spot
    "High-Low"        = "#F39C12",   # orange — isolated high country
    "Low-High"        = "#9B59B6",   # purple — isolated low country
    "Not significant" = "#D5D8DC"    # grey
)

map_lisa <- ggplot(world_clean) +
    geom_sf(
        aes(fill = hotspot_type),
        color     = "white",
        linewidth = 0.1
    ) +
    scale_fill_manual(
        values   = lisa_pal,
        name     = "LISA Cluster Type",
        na.value = "#D9D9D9",
        drop     = FALSE
    ) +
    labs(
        title    = "Global Neonatal Mortality \u2014 Spatial Clusters (LISA), 2024",
        subtitle = paste0(
            "Global Moran's I = ",
            round(moran_result$estimate["Moran I statistic"], 3),
            "  (p < 0.001)"
        ),
        caption  = paste0(
            "Source: UNICEF estimates, 2024. ",
            "Grey = no data or island nations excluded from analysis.\n",
            "High-High = significant cluster of high mortality. ",
            "Low-Low = significant cluster of low mortality."
        )
    ) +
    theme_void(base_family = "sans") +
    theme(
        plot.title      = element_text(size = 11, face = "bold",
                                       hjust = 0.5,
                                       margin = margin(b = 4)),
        plot.subtitle   = element_text(size = 9, hjust = 0.5,
                                       color = "grey40",
                                       margin = margin(b = 6)),
        plot.caption    = element_text(size = 7, color = "grey50",
                                       hjust = 0,
                                       margin = margin(t = 8)),
        legend.position = "bottom",
        legend.title    = element_text(size = 8, face = "bold"),
        legend.text     = element_text(size = 7),
        plot.margin     = margin(10, 10, 10, 10),
        plot.background = element_rect(fill = "white", color = NA)
    ) +
    guides(
        fill = guide_legend(
            nrow           = 1,
            title.position = "top",
            title.hjust    = 0.5
        )
    )

map_lisa

ggsave(
    filename = here("outputs", "figures", "02_neonatal_mortality_LISA.png"),
    plot     = map_lisa,
    width    = 12,
    height   = 7,
    dpi      = 300,
    bg       = "white"
)


# 11. SIGNIFICANT CLUSTERS TABLE ----------------------------------------------

clusters_table <- world_clean |>
    st_drop_geometry() |>
    filter(hotspot_type != "Not significant",
           !is.na(hotspot_type)) |>
    select(
        Country      = name_long,
        Continent    = continent,
        `Mort. Rate` = `2024.5`,
        `Neighbours' Rate` = lag_neo_mort_rat,
        `Cluster Type`     = hotspot_type,
        `Local I`          = lisa_i,
        `p-value`          = lisa_p
    ) |>
    mutate(
        across(where(is.numeric), ~ round(.x, 3))
    ) |>
    arrange(`Cluster Type`, desc(`Mort. Rate`))

clusters_gt <- clusters_table |>
    gt() |>
    tab_header(
        title    = md("**Significant LISA Clusters — Neonatal Mortality**"),
        subtitle = md("*Countries with statistically significant spatial association (p < 0.05)*")
    ) |>
    tab_style(
        style     = list(cell_fill(color = "#1F3864"),
                         cell_text(color = "white", weight = "bold")),
        locations = cells_column_labels()
    ) |>
    tab_style(
        style     = list(cell_fill(color = "#FADBD8"),
                         cell_text(color = "#7B1E1E", weight = "bold")),
        locations = cells_body(rows = `Cluster Type` == "High-High")
    ) |>
    tab_style(
        style     = list(cell_fill(color = "#D6EAF8"),
                         cell_text(color = "#1A5276", weight = "bold")),
        locations = cells_body(rows = `Cluster Type` == "Low-Low")
    ) |>
    tab_style(
        style     = cell_fill(color = "#FEF9E7"),
        locations = cells_body(rows = `Cluster Type` == "High-Low")
    ) |>
    tab_style(
        style     = cell_fill(color = "#F5EEF8"),
        locations = cells_body(rows = `Cluster Type` == "Low-High")
    ) |>
    cols_align(align = "left",   columns = c(Country, Continent, `Cluster Type`)) |>
    cols_align(align = "center", columns = c(`Mort. Rate`, `Neighbours' Rate`,
                                             `Local I`, `p-value`)) |>
    opt_table_font(font = google_font("Source Sans Pro")) |>
    tab_options(
        table.width                = pct(100),
        heading.title.font.size    = px(14),
        heading.subtitle.font.size = px(11),
        table.border.top.color     = "#1F3864",
        table.border.top.width     = px(3),
        row.striping.include_table_body = FALSE
    ) |>
    tab_source_note(
        md("Mort. Rate = neonatal deaths per 1,000 live births.
            Neighbours' Rate = spatially lagged mean of contiguous neighbours.")
    )

clusters_gt

gtsave(clusters_gt,
       filename = here("outputs", "tables", "LISA_significant_clusters.png"))


# =============================================================================
# END OF SCRIPT
# =============================================================================