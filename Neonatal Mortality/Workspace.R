
# Load Libraries ----------------------------------------------------------

pacman::p_load(
    tidyverse,
    rio,
    tmap,
    spdep,
    geodata,
    sf,
    scales,
    summarytools,
    gt,
    here
)


# load data ---------------------------------------------------------------
# path <- "D:/DSO_mac/R Portifolio/Neonatal Mortality/Neonatal Mortality"
# 
# data <- import(file.path(path, "Neonatal_Mortality_Rates_2025-1.xlsx"))
data <- import(here("Neonatal_Mortality_Rates_2025-1.xlsx"))

data("world")


# 2. Export your spatial object (replace 'my_geodata' and 'my_shapefile')
st_write(obj = world, dsn = "D:/DSO_mac/R Portifolio/Neonatal Mortality", layer = "world_map", driver = "ESRI Shapefile")


# Fix factor order for consistent legend
world_clean <- world_clean |>
    mutate(
        hotspot_type = factor(
            hotspot_type,
            levels = c("High-High", "Low-Low",
                       "High-Low",  "Low-High",
                       "Not significant")
        )
    )

# LISA colour palette
lisa_pal <- c(
    "High-High"       = "#E74C3C",
    "Low-Low"         = "#3498DB",
    "High-Low"        = "#F39C12",
    "Low-High"        = "#9B59B6",
    "Not significant" = "#D5D8DC"
)

# Map
ggplot(world_clean) +
    geom_sf(
        aes(fill = hotspot_type),
        color     = "black",
        linewidth = 0.1
    ) +
    scale_fill_manual(
        values   = lisa_pal,
        name     = "LISA Cluster Type",
        na.value = "#D9D9D9",
        drop     = FALSE
    ) +
    labs(
        title    = "Global Neonatal Mortality — Spatial Clusters (LISA), 2024",
        subtitle = paste0(
            "Global Moran's I = ",
            round(moran_result$estimate["Moran I statistic"], 3),
            "  (p < 0.001)"
        ),
        caption  = "Source: UNICEF estimates, 2024. Grey = no data or island nations excluded from analysis."
    ) +
    theme_void(base_family = "sans") +
    theme(
        plot.title       = element_text(size = 11, face = "bold",
                                        hjust = 0.5,
                                        margin = margin(b = 4)),
        plot.subtitle    = element_text(size = 9, hjust = 0.5,
                                        color = "grey40",
                                        margin = margin(b = 6)),
        plot.caption     = element_text(size = 7, color = "grey50",
                                        hjust = 0,
                                        margin = margin(t = 8)),
        legend.position  = "bottom",
        legend.title     = element_text(size = 8, face = "bold"),
        legend.text      = element_text(size = 7),
        plot.margin      = margin(10, 10, 10, 10),
        plot.background  = element_rect(fill = "white", color = NA)
    ) +
    guides(
        fill = guide_legend(
            nrow           = 1,
            title.position = "top",
            title.hjust    = 0.5
        )
    )

# Save
ggsave(
    filename = "outputs/figures/neonatal_LISA_global.png",
    width    = 12,
    height   = 7,
    dpi      = 300,
    bg       = "white"
)

### Joing the 

data <- data |> 
    rename(
        "country" = Country.Name
    )


world_joint <- world |> 
    left_join(
        data |> 
            select(country, `2024.5`),
        by = join_by(name_long == country)
    )

# Visualize the map

tm_shape(world_joint)+
    tm_fill("2024.5",
            pallete = "Reds",
            title = "Neonatal Mortality Globally /n (2024)")

world_joint |> 
    ggplot()+
    geom_sf(
        aes(
            fill = `2024.5`
        ),
        color = "black"
    )+
    scale_fill_gradientn(
        colours = c(
            "#F7F4F9", "#C994C7", "#DF65B0",
            "#E7298A", "#CE1256", "#91003F"
        ),
        breaks  = c(0, 5, 10, 15, 20, 25, 30),
        limits  = c(0, 40),
        na.value = "#D9D9D9",
        labels  = c("0", "5", "10", "15", "20", "25", "30+"),
        name    = "Neonatal mortality rate/n(per 1,000 live births)",
        guide   = guide_colorsteps(
            barwidth       = unit(8, "cm"),
            barheight      = unit(0.4, "cm"),
            title.position = "top",
            title.hjust    = 0.5,
            show.limits    = TRUE
        )
    )+
    labs(
        title  = "Global Neonatal Mortality/n 2024",
        caption  = "Source: UNICEF estimates, 2024. Grey indicates no data."
    )+
    theme_void(base_family = "sans")+
    theme(
        plot.title    = element_text(size = 11, face = "bold", hjust = 0.5,
                                     margin = margin(b = 6)),
        plot.caption  = element_text(size = 7, color = "grey50", hjust = 0,
                                     margin = margin(t = 8)),
        legend.position    = "bottom",
        legend.title       = element_text(size = 8, face = "bold"),
        legend.text        = element_text(size = 7),
        plot.margin        = margin(10, 10, 10, 10),
        plot.background    = element_rect(fill = "white", color = NA)
    )
    



# Global Morans I ---------------------------------------------------------

world_clean <- world_joint |> 
    filter(!is.na(`2024.5`))
all(st_is_valid(world_clean))

nb <- poly2nb(world_clean, queen = T)
listw <- nb2listw(nb, style = "W", zero.policy = T)
summary(nb)


moran_result <- moran.test(
    x = world_clean$`2024.5`,
    listw = listw,
    randomisation = T,
    zero.policy = T
    
)

print(moran_result)

cat("Moran's I = ",round(moran_result$statistic["Moran I statistic"],2),"/n",
    "Expected I = ", round(moran_result$estimate["Expectation"], 2), "/n",
    "p.value = ", round(moran_result$p.value, 2), "/n/n"
    )






# Build the results data frame
moran_table <- tibble(
    Statistic = c(
        "Moran's I",
        "Expected I (H₀)",
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

# Render the gt table
moran_table |>
    gt() |>
    tab_header(
        title    = md("**Global Moran's I Test**"),
        subtitle = md("*Neonatal Mortality Rate — Worldwide, 2024*")
    ) |>
    tab_source_note(
        source_note = md(
            "Weights: Queen contiguity, row-standardised (W). 
       Test: randomisation. Island nations excluded (zero neighbours)."
        )
    ) |>
    tab_style(
        style = list(
            cell_fill(color = "#1F3864"),
            cell_text(color = "white", weight = "bold")
        ),
        locations = cells_column_labels()
    ) |>
    tab_style(
        style = cell_fill(color = "#D6E4F0"),
        locations = cells_body(rows = seq(1, nrow(moran_table), 2))
    ) |>
    tab_style(
        style = cell_text(weight = "bold", color = "#1F3864"),
        locations = cells_body(
            columns = Statistic
        )
    ) |>
    tab_style(
        style = list(
            cell_fill(color = "#D5F5E3"),
            cell_text(color = "#1E5C2E", weight = "bold")
        ),
        locations = cells_body(rows = Statistic == "Interpretation")
    ) |>
    tab_style(
        style = list(
            cell_fill(color = "#D5F5E3"),
            cell_text(color = "#1E5C2E", weight = "bold")
        ),
        locations = cells_body(
            columns = Value,
            rows    = Statistic == "Interpretation"
        )
    ) |>
    cols_width(
        Statistic ~ px(260),
        Value     ~ px(340)
    ) |>
    cols_align(align = "left",  columns = Statistic) |>
    cols_align(align = "right", columns = Value) |>
    opt_table_font(font = google_font("Source Sans Pro")) |>
    opt_stylize(style = 1, color = "blue") |>
    tab_options(
        table.width              = px(600),
        heading.title.font.size  = px(16),
        heading.subtitle.font.size = px(12),
        table.border.top.color   = "#1F3864",
        table.border.top.width   = px(3),
        source_notes.font.size   = px(10)
    )


# Local Morans I


# ── Step 1: Run LISA ───────────────────────────────────────────────────────
lisa <- localmoran(
    x           = world_clean$`2024.5`,
    listw       = listw,
    zero.policy = TRUE,
    na.action   = na.exclude
)

# ── Step 2: Compute spatial lag ────────────────────────────────────────────
# lag.listw gives each country the weighted average
# of its NEIGHBOURS' neonatal mortality rates
world_clean$lag_neo_mort_rat <- lag.listw(
    listw,
    world_clean$`2024.5`,
    zero.policy = TRUE          # ← must be TRUE not empty
)

# ── Step 3: Compute mean — the dividing line between high and low ──────────
mean_mort_rate   <- mean(world_clean$`2024.5`, na.rm = TRUE)
median_mort_rate <- median(world_clean$`2024.5`, na.rm = TRUE)

cat("Mean mortality rate   =", round(mean_mort_rate, 2), "/n",
    "Median mortality rate =", round(median_mort_rate, 2), "/n/n")

# ── Step 4: Attach LISA results and classify ───────────────────────────────
# Note: only ONE mutate() call
# Note: use `2024.5` (the column) not mean_mort_rate (a scalar)

world_clean <- world_clean |>
    mutate(
        lisa_i = lisa[, "Ii"],
        lisa_p = lisa[, "Pr(z != E(Ii))"],
        
        hotspot_type = case_when(
            # High-High: country is HIGH and neighbours are HIGH
            lisa_p < 0.05 & `2024.5` > mean_mort_rate & lag_neo_mort_rat > mean_mort_rate ~ "High-High",
            
            # Low-Low: country is LOW and neighbours are LOW
            lisa_p < 0.05 & `2024.5` < mean_mort_rate & lag_neo_mort_rat < mean_mort_rate ~ "Low-Low",
            
            # High-Low: country is HIGH but neighbours are LOW (isolated hotspot)
            lisa_p < 0.05 & `2024.5` > mean_mort_rate & lag_neo_mort_rat < mean_mort_rate ~ "High-Low",
            
            # Low-High: country is LOW but neighbours are HIGH (possible surveillance gap)
            lisa_p < 0.05 & `2024.5` < mean_mort_rate & lag_neo_mort_rat > mean_mort_rate ~ "Low-High",
            
            # Everything else
            TRUE ~ "Not significant"
        )
    )

# ── Step 5: Check the classification counts ────────────────────────────────
world_clean |>
    st_drop_geometry() |>
    count(hotspot_type, sort = TRUE)








lisa_pal <- c(
    "High-High"       = "#E74C3C",   # red   — high lifeExp, high neighbours
    "Low-Low"         = "#3498DB",   # blue  — low lifeExp,  low neighbours
    "High-Low"        = "#F39C12",   # orange — isolated high country
    "Low-High"        = "#9B59B6",   # purple — isolated low country
    "Not significant" = "#D5D8DC"    # grey
)

tm_shape(world_clean) +
    tm_fill(
        col     = "hotspot_type",
        palette = lisa_pal,
        title   = "LISA Cluster Type"
    ) +
    tm_borders(col = "white", lwd = 0.5) +
    tm_layout(
        main.title      = "Spatial Clustering of Life Expectancy-Africa",
        main.title.size = 1.1,
        legend.outside  = F,
        frame           = FALSE
    ) +
    tm_credits(
        paste0("Global Moran's I = ",
               round(moran_result$estimate["Moran I statistic"], 3),
               "  (p = ", round(moran_result$p.value, 4), ")"),
        position = c("left", "bottom"),
        size = 0.7
    )
