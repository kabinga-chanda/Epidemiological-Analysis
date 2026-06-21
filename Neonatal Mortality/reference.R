pacman::p_load(
    "spData",
    tidyverse,
    sf,
    tmap,
    spdep,
    geodata,
    RDHISdata
)


# Load it
data("world")

# Inspect
names(world)
glimpse(world)

world |>  
    group_by(continent) |> 
    summarise(
        num_country = sum(!is.na(name_long))) |> 
    arrange(desc(num_country
    )) |> 
    ggplot()+
    geom_col(
        aes(y = num_country,
            x = fct_reorder(continent,
                            num_country)),
        fill = "lavender",
        color = "black")+
    coord_flip()+
    labs(
        x = "Continent",
        y = "Number of Coutries",
        title = "Countients By Number of Countries"
    )+
    theme(
        panel.background = element_blank(),
        axis.line.x = element_line(color = "black"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold",
                                  size = 14)
        
    )

# Visualization -----------------------------------------------------------

tm_shape(world)+
    tm_fill(
        "lifeExp",
        pallete = "Y10rRd",
        title = "Life Expectancy"
    )+
    tm_borders(col = "white", lwd = 0.3) +
    tm_layout(main.title = "Global Life Expectancy — Raw Values",
              frame      = FALSE)

# Filter out Africa -------------------------------------------------------

africa <- world |> 
    filter(continent == "Africa")


africa |> 
    st_drop_geometry() |> 
    filter(is.na(lifeExp)) |> 
    select(name_long)


africa_clean <- africa |> filter(!is.na(lifeExp))
all(st_is_valid(africa_clean))

nb <- poly2nb (africa_clean, queen = T)
listw <- nb2listw(nb, style = "W", zero.policy = T)
summary(nb)


islands <- which(card(nb) == 0)
africa_clean$name_long[islands]


# Global Moran's I --------------------------------------------------------

moran_result <-  moran.test(
    x = africa_clean$lifeExp,
    listw = listw,
    randomisation = T,
    zero.policy = T
)

print(moran_result)

cat(
    "Moran's I  =", round(moran_result$estimate["Moran I statistic"], 3), "\n",
    "Expected I =", round(moran_result$estimate["Expectation"], 4),        "\n",
    "p-value    =", format(moran_result$p.value, scientific = F),                       "\n\n"
)


moran.plot(
    x           = africa_clean$lifeExp,
    listw       = listw,
    zero.policy = TRUE,
    labels      = africa_clean$name_long,
    xlab        = "Life Expectancy (standardised)",
    ylab        = "Spatial Lag — Neighbours' Average Life Expectancy",
    main        = "Moran Scatterplot — African Life Expectancy"
)

# Moran's I (LISA) --------------------------------------------------------

lisa <- localmoran(
    x = africa_clean$lifeExp,
    listw = listw,
    zero.policy = T,
    na.action = na.exclude
)

africa_clean$lag_lifeExp <- lag.listw(listw, africa_clean$lifeExp, zero.policy = TRUE)
mean_lifeExp <- mean(africa_clean$lifeExp, na.rm = TRUE)
median_lifeExp <- median(africa_clean$lifeExp, na.rm = T)

cat("Mean African life expectancy:", round(mean_lifeExp, 1), "years")
cat("Median African Life Expectancy:", round(median_lifeExp, 1), "years")




africa_clean <- africa_clean %>%
    mutate(
        lisa_i = lisa[, "Ii"],
        lisa_p = lisa[, "Pr(z != E(Ii))"],
        
        hotspot_type = case_when(
            lisa_p < 0.05 & lifeExp > mean_lifeExp & lag_lifeExp > mean_lifeExp ~ "High-High",
            lisa_p < 0.05 & lifeExp < mean_lifeExp & lag_lifeExp < mean_lifeExp ~ "Low-Low",
            lisa_p < 0.05 & lifeExp > mean_lifeExp & lag_lifeExp < mean_lifeExp ~ "High-Low",
            lisa_p < 0.05 & lifeExp < mean_lifeExp & lag_lifeExp > mean_lifeExp ~ "Low-High",
            TRUE                                                                 ~ "Not significant"
        )
    )



africa_clean %>%
    st_drop_geometry() %>%
    count(hotspot_type, sort = TRUE)







lisa_pal <- c(
    "High-High"       = "#E74C3C",   # red   — high lifeExp, high neighbours
    "Low-Low"         = "#3498DB",   # blue  — low lifeExp,  low neighbours
    "High-Low"        = "#F39C12",   # orange — isolated high country
    "Low-High"        = "#9B59B6",   # purple — isolated low country
    "Not significant" = "#D5D8DC"    # grey
)

tm_shape(africa_clean) +
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

zambia <- gadm(country = "ZMB", level = 2, path = tempdir()) |>
    st_as_sf()


nrow(zambia)
plot(st_geometry(zambia))

zambia |> 
    group_by(NAME_1) |> 
    summarise(
        num_dist = sum(!is.na(NAME_2))
    ) |> 
    arrange(desc(num_dist))



