

# Install pak — a better package installer for Windows
install.packages("pak")


# ── Load every session ─────────────────────────────────────


pacman::p_load(smerc, 
               sf,
               dplyr, 
               ggplot2, 
               readr, 
               terra,
               geodata,
               tmap)

library(tmap)

install.packages("sf")
install.packages("tmap", type = "binary") 

library(sf)



# Load data ---------------------------------------------------------------



print(pm25)


# Convert kg/m³ to µg/m³ (multiply by 1 billion)
pm25_ugm3 <- pm25 * 1e9

# Calculate annual mean across all 365 days
pm25_annual <- mean(pm25_ugm3)

# Check values look realistic (Zambia typically 15–60 µg/m³)
summary(pm25_annual)
plot(pm25_annual, main = "Annual Mean PM2.5 Zambia 2025 (µg/m³)")



zambia_dist <- gadm(country = "ZMB", level = 2, path = tempdir())
zambia_dist_sf <- sf::st_as_sf(zambia_dist)
plot(zambia_dist_sf)
cat("Number of districts:", nrow(zambia_dist_sf), "\n")
names(zambia_dist_sf)
head(zambia_dist_sf[, c("NAME_1", "NAME_2")])


# Extract mean PM2.5 per district
zonal <- terra::extract(
    pm25_annual,
    vect(zambia_dist_sf),
    fun   = mean,
    na.rm = TRUE
)

# Attach to districts
zambia_dist_sf$pm25_mean <- round(zonal[, 2], 2)

# View results — paste this to me!
# Simple fix — use as.data.frame instead
result_table <- zambia_dist_sf |>
    sf::st_drop_geometry() |>
    dplyr::select(NAME_1, NAME_2, pm25_mean) |>
    dplyr::arrange(desc(pm25_mean))

# Print top 20
print(result_table[1:20, ], row.names = FALSE)

# Make sure extraction worked
summary(zambia_dist_sf$pm25_mean)
range(zambia_dist_sf$pm25_mean, na.rm = TRUE)




df <- as.data.frame(zambia_dist_sf[, c("NAME_1", "NAME_2", "pm25_mean")])
df <- df[order(-df$pm25_mean), ]
head(df, 20)



# Satscan -----------------------------------------------------------------

library(smerc)
library(sf)

# Project to UTM Zone 36S (metres — required by smerc)
zambia_proj <- sf::st_transform(zambia_dist_sf, crs = 32736)

# Get district centroids
centroids <- sf::st_centroid(zambia_proj)
coords    <- sf::st_coordinates(centroids)

# Define cases — districts exceeding WHO IT3 (15 µg/m³)
zambia_dist_sf$cases <- as.integer(zambia_dist_sf$pm25_mean > 15)
cat("Districts exceeding 15 µg/m³:", sum(zambia_dist_sf$cases), "\n")

# You need population — add it
# For now use equal population as placeholder



# Use continuous PM2.5 values
zambia_dist_sf$cases_cont <- as.integer(round(zambia_dist_sf$pm25_mean * 100))
zambia_dist_sf$pop        <- rep(100000L, nrow(zambia_dist_sf))

# Check it looks right
cat("Min cases:", min(zambia_dist_sf$cases_cont), "\n")
cat("Max cases:", max(zambia_dist_sf$cases_cont), "\n")
cat("Sample values:", head(zambia_dist_sf$cases_cont), "\n")









# Project and get coordinates
zambia_proj <- sf::st_transform(zambia_dist_sf, crs = 32736)
centroids   <- sf::st_centroid(zambia_proj)
coords      <- sf::st_coordinates(centroids)

# Run scan
set.seed(123)
result <- scan.test(
    coords = coords,
    cases  = zambia_dist_sf$cases_cont,
    pop    = zambia_dist_sf$pop,
    nsim   = 999,
    alpha  = 0.05,
    ubpop  = 0.5
)

summary(result)


pm25_mean = 16.64 µg/m³  →  cases_cont = 1664
pm25_mean = 10.46 µg/m³  →  cases_cont = 1046
pm25_mean =  6.80 µg/m³  →  cases_cont =  680

# Step 1 - Create continuous cases
zambia_dist_sf$cases_cont <- as.integer(round(zambia_dist_sf$pm25_mean * 100))
zambia_dist_sf$pop        <- rep(100000L, nrow(zambia_dist_sf))

# Step 2 - Check values
cat("Min:", min(zambia_dist_sf$cases_cont), "\n")
cat("Max:", max(zambia_dist_sf$cases_cont), "\n")


# Step 3 - Project and get coordinates
zambia_proj <- sf::st_transform(zambia_dist_sf, crs = 32736)
centroids   <- sf::st_centroid(zambia_proj)
coords      <- sf::st_coordinates(centroids)


# Step 4 - Run scan
set.seed(123)
result <- scan.test(
    coords = coords,
    cases  = zambia_dist_sf$cases_cont,
    pop    = zambia_dist_sf$pop,
    nsim   = 999,
    alpha  = 0.05,
    ubpop  = 0.5
)

summary(result)





# See which districts fall in each cluster
for (i in 1:length(result$clusters)) {
    cl  <- result$clusters[[i]]
    ids <- cl$locids
    cat("\n=== Cluster", i, "===\n")
    cat("P-value:", cl$pvalue, "\n")
    cat("Rate Ratio:", round(cl$smr, 2), "\n")
    cat("Districts:\n")
    print(zambia_dist_sf[ids, c("NAME_1", "NAME_2", "pm25_mean")] |>
              sf::st_drop_geometry())
}






library(tmap)

# Tag districts with cluster membership
zambia_dist_sf$cluster <- "No Cluster"

sig <- which(sapply(result$clusters, function(x) x$pvalue <= 0.05))

for (i in sig) {
    ids <- result$clusters[[i]]$locids
    zambia_dist_sf$cluster[ids] <- paste0("Cluster ", i)
}

# Interactive map
tmap_mode("view")

tm_shape(zambia_dist_sf) +
    tm_polygons(
        col        = "cluster",
        palette    = c("No Cluster" = "grey80",
                       "Cluster 1"  = "red",
                       "Cluster 2"  = "orange", 
                       "Cluster 3"  = "yellow3",
                       "Cluster 4"  = "blue"),
        title      = "PM2.5 Clusters",
        popup.vars = c("Province" = "NAME_1",
                       "District" = "NAME_2",
                       "PM2.5"    = "pm25_mean",
                       "Cluster"  = "cluster")
    ) +
    tm_layout(title = "PM2.5 Spatial Clusters — Zambia 2025")





library(ggplot2)
library(sf)
library(dplyr)

# ── Tag clusters ───────────────────────────────────────────
zambia_dist_sf$cluster <- "No Cluster"
for (i in 1:4) {
    ids <- result$clusters[[i]]$locids
    zambia_dist_sf$cluster[ids] <- paste0("Cluster ", i)
}

# Add RR labels for legend
zambia_dist_sf$cluster_label <- case_when(
    zambia_dist_sf$cluster == "Cluster 1" ~ "Cluster 1 (RR=1.40, p=0.001)",
    zambia_dist_sf$cluster == "Cluster 2" ~ "Cluster 2 (RR=1.20, p=0.001)",
    zambia_dist_sf$cluster == "Cluster 3" ~ "Cluster 3 (RR=1.10, p=0.001)",
    zambia_dist_sf$cluster == "Cluster 4" ~ "Cluster 4 (RR=1.20, p=0.001)",
    TRUE ~ "No Significant Cluster"
)

# ── Neighbour countries outline ────────────────────────────
library(geodata)
neighbours <- c("COD", "TZA", "MWI", "MOZ", "ZWE", 
                "BWA", "NAM", "AGO")
neighbour_sf <- do.call(rbind, lapply(neighbours, function(c) {
    st_as_sf(gadm(country = c, level = 0, path = tempdir()))
}))

# ── Zambia outline ─────────────────────────────────────────
zambia_outline <- st_union(zambia_dist_sf)

# ── Colour palette ─────────────────────────────────────────
cluster_colours <- c(
    "Cluster 1 (RR=1.40, p=0.001)" = "#D7191C",
    "Cluster 2 (RR=1.20, p=0.001)" = "#F46D43",
    "Cluster 3 (RR=1.10, p=0.001)" = "#FDAE61",
    "Cluster 4 (RR=1.20, p=0.001)" = "#2166AC",
    "No Significant Cluster"        = "#D9D9D9"
)

# ── Build map ──────────────────────────────────────────────
map <- ggplot() +
    
    # Neighbour countries
    geom_sf(data = neighbour_sf,
            fill = "#F5F5F0", colour = "#BBBBBB", linewidth = 0.3) +
    
    # District polygons coloured by cluster
    geom_sf(data = zambia_dist_sf,
            aes(fill = cluster_label),
            colour = "white", linewidth = 0.15) +
    
    # Zambia border outline
    geom_sf(data = zambia_outline,
            fill = NA, colour = "#333333", linewidth = 0.8) +
    
    # Colours
    scale_fill_manual(
        values = cluster_colours,
        name   = "PM2.5 Spatial Clusters",
        guide  = guide_legend(
            title.position = "top",
            keywidth  = unit(0.5, "cm"),
            keyheight = unit(0.5, "cm")
        )
    ) +
    
    # Map extent — Zambia + buffer
    coord_sf(xlim = c(20, 36), ylim = c(-19, -7), expand = FALSE) +
    
    # Labels
    labs(
        title    = "Spatial Clustering of PM2.5 in Zambia, 2025",
        subtitle = "Kulldorff Spatial Scan Statistic | CAMS Satellite Data",
        caption  = "Data: Copernicus Atmosphere Monitoring Service (CAMS), 2025\nMethod: Kulldorff Scan Statistic (smerc), 999 Monte Carlo simulations\nProjection: WGS 84"
    ) +
    
    # Clean theme
    theme_void(base_family = "sans") +
    theme(
        # Title
        plot.title    = element_text(size = 14, face = "bold",
                                     hjust = 0.5, margin = margin(b = 4)),
        plot.subtitle = element_text(size = 10, hjust = 0.5,
                                     color = "grey40", margin = margin(b = 8)),
        plot.caption  = element_text(size = 7, color = "grey50",
                                     hjust = 0, margin = margin(t = 8)),
        
        # Legend
        # Move legend to bottom
        legend.position   = "bottom",
        legend.box        = "horizontal",
        legend.title      = element_text(size = 9, face = "bold",
                                         hjust = 0.5),
        legend.text       = element_text(size = 8.5),
        legend.margin     = margin(t = 5, b = 5),
        legend.key        = element_rect(colour = "grey60",
                                         linewidth = 0.3),
        
        # Panel
        plot.background  = element_rect(fill = "white", colour = NA),
        plot.margin      = margin(10, 10, 10, 10)
    )

print(map)





library(ggplot2)
library(sf)
library(dplyr)
library(geodata)

# ── Neighbour countries ────────────────────────────────────
neighbours <- c("COD", "TZA", "MWI", "MOZ", "ZWE",
                "BWA", "NAM", "AGO")
neighbour_sf <- do.call(rbind, lapply(neighbours, function(c) {
    st_as_sf(gadm(country = c, level = 0, path = tempdir()))
}))

# ── Zambia outline ─────────────────────────────────────────
zambia_outline <- st_union(zambia_dist_sf)

# ── Tag clusters ───────────────────────────────────────────
zambia_dist_sf$cluster <- "Not significant"
for (i in 1:4) {
    ids <- result$clusters[[i]]$locids
    zambia_dist_sf$cluster[ids] <- paste0("Cluster ", i)
}

# ── Factor levels for legend order ────────────────────────
zambia_dist_sf$cluster <- factor(
    zambia_dist_sf$cluster,
    levels = c("Cluster 1", "Cluster 2",
               "Cluster 3", "Cluster 4",
               "Not significant")
)

# ── Colours matching reference style ──────────────────────
cluster_colours <- c(
    "Cluster 1"       = "#D7191C",   # strong red
    "Cluster 2"       = "#F17720",   # orange
    "Cluster 3"       = "#FDCC8A",   # light orange
    "Cluster 4"       = "#2166AC",   # blue
    "Not significant" = "#D3D3D3"    # grey
)

# ── Build map ──────────────────────────────────────────────
map <- ggplot() +
    
    # Neighbour countries background
    geom_sf(data = neighbour_sf,
            fill = "#EEEEEE", colour = "#AAAAAA",
            linewidth = 0.3) +
    
    # Districts coloured by cluster
    geom_sf(data = zambia_dist_sf,
            aes(fill = cluster),
            colour = "white", linewidth = 0.2) +
    
    # Bold Zambia border
    geom_sf(data = zambia_outline,
            fill = NA, colour = "#222222",
            linewidth = 0.9) +
    
    # Colours + legend
    scale_fill_manual(
        values = cluster_colours,
        name   = "PM2.5 Cluster",
        drop   = FALSE,
        guide  = guide_legend(
            direction      = "horizontal",
            title.position = "top",
            title.hjust    = 0.5,
            nrow           = 1,
            keywidth       = unit(1.0, "cm"),
            keyheight      = unit(0.5, "cm"),
            label.position = "bottom"
        )
    ) +
    
    # Map extent
    coord_sf(xlim = c(20, 36),
             ylim = c(-19.5, -6.5),
             expand = FALSE) +
    
    # Titles
    labs(
        title    = "Spatial Clustering of PM\u2082.\u2085 in Zambia, 2025",
        subtitle = "Kulldorff Spatial Scan Statistic  |  All clusters: p = 0.001",
        caption  = "Source: Copernicus Atmosphere Monitoring Service (CAMS), 2025.\nGrey = no significant cluster detected.\nMethod: Kulldorff Scan Statistic, 999 Monte Carlo simulations (smerc, R)."
    ) +
    
    # Clean theme matching reference
    theme_void(base_size = 11, base_family = "sans") +
    theme(
        # Titles
        plot.title    = element_text(size   = 15, face = "bold",
                                     hjust  = 0.5,
                                     margin = margin(b = 4)),
        plot.subtitle = element_text(size   = 10, hjust = 0.5,
                                     colour = "grey40",
                                     margin = margin(b = 10)),
        plot.caption  = element_text(size   = 7.5, hjust = 0,
                                     colour = "grey50",
                                     margin = margin(t = 10)),
        
        # Legend at bottom — matching reference map style
        legend.position   = "bottom",
        legend.box        = "horizontal",
        legend.title      = element_text(size = 9, face = "bold",
                                         hjust = 0.5),
        legend.text       = element_text(size = 8.5),
        legend.margin     = margin(t = 10, b = 5),
        legend.key        = element_rect(colour = "grey60",
                                         linewidth = 0.3),
        legend.spacing.x  = unit(0.8, "cm"),
        
        # Background
        plot.background  = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "#EAF3FB", colour = NA),
        plot.margin      = margin(10, 15, 10, 15)
    )

print(map)

