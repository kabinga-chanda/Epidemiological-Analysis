# ============================================================
#  SPATIAL CLUSTERING OF PM2.5 IN ZAMBIA — COMPLETE SCRIPT
#  MSc Epidemiology & Biostatistics | University of Zambia
#  Data: CAMS Satellite PM2.5, 2025
#  Method: Kulldorff Spatial Scan Statistic (smerc)
#  Author: [Your Name]
#  Date: 2025
# ============================================================


# ============================================================
# SECTION 1: INSTALL AND LOAD PACKAGES
# ============================================================
# We only need to install packages once. After that, we just
# load them with library() at the start of every session.
# ─────────────────────────────────────────────────────────────

# Install all required packages (run once only)
# install.packages(c(
#   "terra",    # read and process raster/satellite data
#   "sf",       # handle vector spatial data (shapefiles, polygons)
#   "geodata",  # download country boundary data automatically
#   "smerc",    # Kulldorff spatial scan statistic
#   "tmap",     # interactive maps
#   "ggplot2",  # publication-quality static maps
#   "dplyr",    # data manipulation
#   "readr"     # read/write CSV files
# ))

# Load all packages
library(terra)    # raster processing
library(sf)       # vector spatial data
library(geodata)  # country boundaries
library(smerc)    # spatial scan statistic
library(tmap)     # interactive mapping
library(ggplot2)  # static mapping
library(dplyr)    # data wrangling
library(readr)    # CSV import/export


# ============================================================
# SECTION 2: LOAD YOUR PM2.5 SATELLITE DATA
# ============================================================
# The CAMS NetCDF file contains daily PM2.5 values for the
# whole of Zambia from January 1 to December 31, 2025.
# Each of the 365 layers = one day of PM2.5 measurements.
# Units are in kg/m³ — we will convert to µg/m³ later.
# ─────────────────────────────────────────────────────────────

# Load the NetCDF file downloaded from Copernicus CAMS
# terra::rast() reads raster files including NetCDF format
pm25_raw <- rast("D:/DSO_mac/R Portifolio/Partuculate Matter Zambia/PM2.5/data/3872a2a4de5c3da03a28113624103e74/data_sfc.nc")

# Inspect the loaded data
# This tells us: dimensions, resolution, extent, units, time steps
print(pm25_raw)

# What you should see:
# - 365 layers (one per day)
# - varname: pm2p5
# - unit: kg m**-3
# - extent covering Zambia


# ============================================================
# SECTION 3: CONVERT UNITS AND COMPUTE ANNUAL MEAN
# ============================================================
# CAMS stores PM2.5 in kg/m³ but all health guidelines use
# µg/m³ (micrograms per cubic metre). To convert:
#   1 kg/m³ = 1,000,000,000 µg/m³ (multiply by 1 billion)
#
# We then average all 365 daily layers to get one annual
# mean value per pixel across Zambia.
# ─────────────────────────────────────────────────────────────

# Step 3a: Convert units from kg/m³ to µg/m³
pm25_ugm3 <- pm25_raw * 1e9   # 1e9 = 1,000,000,000

# Step 3b: Compute annual mean across all 365 days
# mean() on a SpatRaster averages across all layers (days)
pm25_annual <- mean(pm25_ugm3)

# Check values are realistic (Zambia typically 6–20 µg/m³)
summary(pm25_annual)

# Quick visual check of the raster
plot(pm25_annual,
     main = "Annual Mean PM2.5 Zambia 2025 (µg/m³)",
     col  = hcl.colors(50, "YlOrRd", rev = TRUE))


# ============================================================
# SECTION 4: LOAD ZAMBIA DISTRICT BOUNDARIES
# ============================================================
# We need a shapefile of Zambia's 115 districts to:
# (a) extract average PM2.5 per district (zonal statistics)
# (b) use as the spatial units for the scan statistic
#
# gadm() downloads boundaries directly from the GADM database.
# level = 2 gives us districts (level 1 = provinces)
# ─────────────────────────────────────────────────────────────

# Download Zambia district boundaries automatically
# ZMB = ISO3 country code for Zambia
# level = 2 = district level
zambia_dist <- gadm(country = "ZMB", level = 2, path = tempdir())

# Convert from terra SpatVector to sf object
# sf format is easier to work with for mapping and analysis
zambia_dist_sf <- st_as_sf(zambia_dist)

# Verify what we have
cat("Number of districts:", nrow(zambia_dist_sf), "\n")
cat("Column names:", names(zambia_dist_sf), "\n")

# Preview district and province names
# NAME_1 = Province, NAME_2 = District
head(zambia_dist_sf[, c("NAME_1", "NAME_2")])

# Quick map to confirm districts look correct
plot(st_geometry(zambia_dist_sf),
     main = "Zambia Districts (n=115)")


# ============================================================
# SECTION 5: ZONAL STATISTICS — EXTRACT PM2.5 PER DISTRICT
# ============================================================
# Zonal statistics = calculate summary values (e.g. mean)
# of a raster within each polygon zone (district).
#
# This is equivalent to the "Zonal Statistics" tool in QGIS.
# terra::extract() is the R equivalent.
#
# We extract the mean PM2.5 value for every pixel that falls
# within each district polygon, then average those values.
# ─────────────────────────────────────────────────────────────

# Extract mean PM2.5 per district
# vect() converts sf object to terra format for extraction
# fun = mean → compute average of all pixels per district
# na.rm = TRUE → ignore missing values (e.g. cloud cover)
# exact = TRUE → account for partial pixels at borders
zonal_stats <- terra::extract(
  pm25_annual,
  vect(zambia_dist_sf),
  fun   = mean,
  na.rm = TRUE,
  exact = TRUE
)

# Attach PM2.5 values back to the district shapefile
zambia_dist_sf$pm25_mean <- round(zonal_stats[, 2], 2)

# Check results
summary(zambia_dist_sf$pm25_mean)

# View top 20 highest PM2.5 districts
results_table <- zambia_dist_sf |>
  st_drop_geometry() |>
  select(NAME_1, NAME_2, pm25_mean) |>
  arrange(desc(pm25_mean))

print(results_table, n = 20)

# How many districts exceed WHO thresholds?
cat("\nDistricts exceeding WHO guideline (5 µg/m³):",
    sum(zambia_dist_sf$pm25_mean > 5), "of 115\n")

cat("Districts exceeding WHO Interim Target 3 (15 µg/m³):",
    sum(zambia_dist_sf$pm25_mean > 15), "of 115\n")


# ============================================================
# SECTION 6: PREPARE DATA FOR SPATIAL SCAN
# ============================================================
# The smerc scan.test() function requires:
# (1) coords  — projected X,Y coordinates (in metres, not degrees)
# (2) cases   — integer count of "events" per district
# (3) pop     — population at risk per district
#
# IMPORTANT: smerc needs coordinates in metres (projected CRS),
# NOT in degrees (geographic CRS like WGS84).
# We use UTM Zone 36S (EPSG:32736) which is correct for Zambia.
#
# For cases: we multiply PM2.5 × 100 to convert continuous
# values to integers while preserving relative differences.
# This lets us use PM2.5 directly as the outcome variable.
# ─────────────────────────────────────────────────────────────

# Step 6a: Project districts to UTM Zone 36S (metres)
# This converts from degrees (WGS84) to metres (UTM)
zambia_proj <- st_transform(zambia_dist_sf, crs = 32736)

# Step 6b: Compute district centroids
# The scan places a circular window centred on each district
# centroid, so we need the centre point of each district
centroids <- st_centroid(zambia_proj)
coords    <- st_coordinates(centroids)  # extract X, Y as matrix

# Verify coordinates are now in metres (large numbers)
head(coords)

# Step 6c: Create integer cases from continuous PM2.5
# Multiply by 100 to preserve 2 decimal places as integer
# e.g. 16.64 µg/m³ → 1664 cases
zambia_dist_sf$cases_cont <- as.integer(
  round(zambia_dist_sf$pm25_mean * 100)
)

# Step 6d: Set equal population for each district
# Using equal population means the scan identifies clusters
# based purely on PM2.5 concentration differences
# Replace with real population data if available
zambia_dist_sf$pop <- rep(100000L, nrow(zambia_dist_sf))

# Quick sanity check
cat("Cases range:", range(zambia_dist_sf$cases_cont), "\n")
cat("Total districts:", nrow(zambia_dist_sf), "\n")


# ============================================================
# SECTION 7: RUN THE KULLDORFF SPATIAL SCAN STATISTIC
# ============================================================
# The Kulldorff scan statistic works by:
# 1. Placing a circular window at every district location
# 2. Varying the window size from 1 district up to 50% of
#    the total population (controlled by ubpop = 0.5)
# 3. Computing a likelihood ratio for each window:
#    how much more PM2.5 is inside vs outside the window?
# 4. The window with the highest likelihood ratio =
#    the most likely cluster
# 5. Statistical significance is tested using Monte Carlo
#    simulation (nsim = 999 random datasets generated and
#    compared to the observed data)
# 6. If p < 0.05 → cluster is statistically significant
#
# set.seed() ensures reproducible results
# ─────────────────────────────────────────────────────────────

set.seed(123)  # for reproducibility

result <- scan.test(
  coords = coords,              # district centroids (metres)
  cases  = zambia_dist_sf$cases_cont, # PM2.5 × 100 as integer
  pop    = zambia_dist_sf$pop,  # population per district
  nsim   = 999,                 # Monte Carlo simulations
  alpha  = 0.05,                # significance threshold
  ubpop  = 0.5                  # max 50% of population in window
)

# View summary of all detected clusters
summary(result)

# Column meanings:
# nregions  = number of districts in the cluster
# max_dist  = radius of the cluster (metres)
# cases     = observed PM2.5 (× 100) inside cluster
# ex        = expected PM2.5 if no clustering
# rr        = rate ratio (observed/expected) — >1 means hotspot
# stat      = likelihood ratio test statistic
# p         = Monte Carlo p-value


# ============================================================
# SECTION 8: INTERPRET CLUSTER RESULTS
# ============================================================
# Extract detailed information about each significant cluster:
# which districts are in it, what province, PM2.5 values, etc.
# ─────────────────────────────────────────────────────────────

# Print detailed cluster membership
cat("\n========================================\n")
cat("  PM2.5 SPATIAL CLUSTER RESULTS\n")
cat("  Zambia, 2025 | Kulldorff Scan Statistic\n")
cat("========================================\n")

for (i in seq_along(result$clusters)) {
  cl <- result$clusters[[i]]

  if (cl$pvalue > 0.05) next  # skip non-significant clusters

  cat("\n--- Cluster", i, "---\n")
  cat("Significant:", cl$pvalue <= 0.05, "\n")
  cat("P-value:", cl$pvalue, "\n")
  cat("Rate Ratio (RR):", round(cl$smr, 3),
      "→ PM2.5 is", round((cl$smr - 1) * 100, 0),
      "% higher than expected\n")
  cat("Number of districts:", length(cl$locids), "\n")
  cat("Cluster radius:", round(cl$max_dist / 1000, 0), "km\n")

  # Get district details for this cluster
  dist_in_cluster <- zambia_dist_sf[cl$locids, ] |>
    st_drop_geometry() |>
    select(NAME_1, NAME_2, pm25_mean) |>
    arrange(desc(pm25_mean))

  cat("Provinces included:",
      paste(unique(dist_in_cluster$NAME_1), collapse = ", "), "\n")
  cat("Districts and PM2.5 values:\n")
  print(dist_in_cluster, row.names = FALSE)
}


# ============================================================
# SECTION 9: TAG DISTRICTS WITH CLUSTER MEMBERSHIP
# ============================================================
# Add a cluster label to each district in the shapefile
# so we can colour them on a map
# ─────────────────────────────────────────────────────────────

# Start with all districts as "No Significant Cluster"
zambia_dist_sf$cluster <- "No Significant Cluster"

# Tag significant clusters
sig_clusters <- which(
  sapply(result$clusters, function(x) x$pvalue <= 0.05)
)

for (i in sig_clusters) {
  ids <- result$clusters[[i]]$locids
  zambia_dist_sf$cluster[ids] <- paste0("Cluster ", i)
}

# Add RR and p-value to legend labels
zambia_dist_sf$cluster_label <- case_when(
  zambia_dist_sf$cluster == "Cluster 1" ~
    paste0("Cluster 1 (RR=",
           round(result$clusters[[1]]$smr, 2), ", p=0.001)"),
  zambia_dist_sf$cluster == "Cluster 2" ~
    paste0("Cluster 2 (RR=",
           round(result$clusters[[2]]$smr, 2), ", p=0.001)"),
  zambia_dist_sf$cluster == "Cluster 3" ~
    paste0("Cluster 3 (RR=",
           round(result$clusters[[3]]$smr, 2), ", p=0.001)"),
  zambia_dist_sf$cluster == "Cluster 4" ~
    paste0("Cluster 4 (RR=",
           round(result$clusters[[4]]$smr, 2), ", p=0.001)"),
  TRUE ~ "No Significant Cluster"
)

# Check counts
table(zambia_dist_sf$cluster)


# ============================================================
# SECTION 10: PUBLICATION-QUALITY MAP (ggplot2)
# ============================================================
# Build a clean, professional map suitable for MSc dissertation
# and journal submission at 300 DPI resolution
# ─────────────────────────────────────────────────────────────

# Step 10a: Download neighbouring countries for context
neighbours <- c("COD", "TZA", "MWI", "MOZ",
                "ZWE", "BWA", "NAM", "AGO")

neighbour_sf <- do.call(rbind, lapply(neighbours, function(c) {
  st_as_sf(gadm(country = c, level = 0, path = tempdir()))
}))

# Step 10b: Create Zambia outer boundary
zambia_outline <- st_union(zambia_dist_sf)

# Step 10c: Define colour palette
cluster_colours <- c(
  "Cluster 1 (RR=1.4, p=0.001)"  = "#D7191C",  # dark red
  "Cluster 2 (RR=1.2, p=0.001)"  = "#F46D43",  # orange
  "Cluster 3 (RR=1.1, p=0.001)"  = "#FDAE61",  # light orange
  "Cluster 4 (RR=1.2, p=0.001)"  = "#2166AC",  # blue
  "No Significant Cluster"        = "#D3D3D3"   # grey
)

# Update cluster_label to match colour palette keys exactly
zambia_dist_sf$cluster_label <- case_when(
  zambia_dist_sf$cluster == "Cluster 1" ~ "Cluster 1 (RR=1.4, p=0.001)",
  zambia_dist_sf$cluster == "Cluster 2" ~ "Cluster 2 (RR=1.2, p=0.001)",
  zambia_dist_sf$cluster == "Cluster 3" ~ "Cluster 3 (RR=1.1, p=0.001)",
  zambia_dist_sf$cluster == "Cluster 4" ~ "Cluster 4 (RR=1.2, p=0.001)",
  TRUE ~ "No Significant Cluster"
)

# Step 10d: Build the map
map <- ggplot() +

  # Layer 1: Neighbouring countries (background context)
  geom_sf(data     = neighbour_sf,
          fill     = "#F0F0F0",
          colour   = "black",
          linewidth = 0.3) +

  # Layer 2: Zambia districts coloured by cluster
  geom_sf(data     = zambia_dist_sf,
          aes(fill = cluster_label),
          colour   = "black",
          linewidth = 0.2) +

  # Layer 3: Bold Zambia border on top
  geom_sf(data     = zambia_outline,
          fill     = NA,
          colour   = "#222222",
          linewidth = 0.9) +

  # Colour scale with legend
  scale_fill_manual(
    values = cluster_colours,
    name   = "PM2.5 Spatial Clusters",
    drop   = FALSE,
    guide  = guide_legend(
      title.position = "top",
      title.hjust    = 0.5,
      keywidth       = unit(0.6, "cm"),
      keyheight      = unit(0.5, "cm")
    )
  ) +

  # Set map extent to show all of Zambia with small buffer
  coord_sf(xlim = c(20, 36), ylim = c(-19.5, -6.5),
           expand = FALSE) +

  # Titles and caption
  labs(
    title    = "Spatial Clustering of PM2.5 in Zambia, 2025",
    subtitle = paste0("Kulldorff Spatial Scan Statistic  |  ",
                      "CAMS Satellite Data  |  ",
                      "999 Monte Carlo Simulations"),
    caption  = paste0(
      "Data: Copernicus Atmosphere Monitoring Service (CAMS), 2025.\n",
      "Method: Kulldorff Spatial Scan Statistic (smerc package, R).\n",
      "Grey districts = no statistically significant cluster (p > 0.05).\n",
      "Projection: WGS 84 | RR = Rate Ratio (Observed/Expected PM2.5)"
    )
  )  +

  # Clean publication theme
  theme_void(base_size = 11, base_family = "sans") +
  theme(
    # Main title
    plot.title    = element_text(size   = 15,
                                 face   = "bold",
                                 hjust  = 0.5,
                                 margin = margin(b = 4)),
    # Subtitle
    plot.subtitle = element_text(size   = 9,
                                 hjust  = 0.5,
                                 colour = "grey40",
                                 margin = margin(b = 10)),
    # Caption (bottom left)
    plot.caption  = element_text(size   = 7,
                                 hjust  = 0,
                                 colour = "grey50",
                                 margin = margin(t = 10)),
    # Legend at bottom
    legend.position  = "bottom",
    legend.box       = "horizontal",
    legend.title     = element_text(size = 9, face = "bold",
                                    hjust = 0.5),
    legend.text      = element_text(size = 8),
    legend.margin    = margin(t = 8, b = 5),
    legend.key       = element_rect(colour   = "grey70",
                                    linewidth = 0.3),
    # Light blue ocean/background
    panel.background = element_rect(fill = "#EAF3FB",
                                    colour = NA),
    # White page background
    plot.background  = element_rect(fill = "white",   colour = NA),
    plot.margin      = margin(10, 15, 10, 15)
  )
# Display map
print(map)


# ============================================================
# SECTION 11: SAVE ALL OUTPUTS
# ============================================================
# Save the map as a high-resolution PNG (300 DPI) suitable
# for dissertation and journal submission.
# Save results tables as CSV for reporting.
# ─────────────────────────────────────────────────────────────

# Define output folder (change to your path)
out_dir <- "D:/DSO_mac/R Portifolio/Partuculate Matter Zambia/PM2.5/"

# Save map at 300 DPI
ggsave(
  filename = paste0(out_dir, "zambia_pm25_clusters_2025.png"),
  plot     = map,
  width    = 10,
  height   = 9,
  dpi      = 300,
  bg       = "white"
)
cat("Map saved!\n")

# Save cluster summary table
cluster_summary <- do.call(rbind, lapply(
  seq_along(result$clusters), function(i) {
    cl <- result$clusters[[i]]
    data.frame(
      cluster      = i,
      n_districts  = length(cl$locids),
      provinces    = paste(
        unique(zambia_dist_sf$NAME_1[cl$locids]),
        collapse   = "; "),
      obs_pm25     = round(cl$cases / 100, 2),
      exp_pm25     = round(cl$expected / 100, 2),
      rate_ratio   = round(cl$smr, 3),
      pvalue       = cl$pvalue,
      significant  = cl$pvalue <= 0.05
    )
  }))

write_csv(cluster_summary,
          paste0(out_dir, "cluster_summary.csv"))
cat("Cluster summary saved!\n")

# Save district-level data with PM2.5 and cluster tags
district_export <- zambia_dist_sf |>
  st_drop_geometry() |>
  select(NAME_1, NAME_2, pm25_mean, cluster, cluster_label)

write_csv(district_export,
          paste0(out_dir, "districts_pm25_clusters.csv"))
cat("District data saved!\n")

# ============================================================
# END OF SCRIPT
# ============================================================
# Output files:
#   zambia_pm25_clusters_2025.png  — publication map (300 DPI)
#   cluster_summary.csv            — cluster statistics table
#   districts_pm25_clusters.csv    — district PM2.5 + clusters
# ============================================================
