#final figure EOT map for 2027 stop-MDA
# Figure 7 Togo
###final code with centroids
library(geodata)
library(sf)
library(ggplot2)
library(terra)
library(dplyr)
library(patchwork)
library(ggforce)
library(ggnewscale)
library(scatterpie)

###############################################################################
# 2. Download Togo Administrative Boundaries
###############################################################################
# Level 2 = Prefecture
togo_prefectures <- gadm(country = "TGO", level = 2, path = tempdir())
# Level 1 = Region (for thicker region borders)
togo_regions     <- gadm(country = "TGO", level = 1, path = tempdir())

# Convert SpatVector to sf
togo_sf         <- st_as_sf(togo_prefectures)
togo_regions_sf <- st_as_sf(togo_regions)

###############################################################################
# 3. Merge Selected Prefectures & Assign Likelihood
###############################################################################
togo_sf_B <- togo_sf %>%
  mutate(merged_prefecture = case_when(
    NAME_2 %in% c("Cinkassé", "Tône") ~ "Tône",      # Merge Cinkassé into Tône
    NAME_2 %in% c("Naki-Ouest", "Kpendjal") ~ "Kpendjal",
    NAME_2 %in% c("Oti", "Oti-Sud") ~ "Oti",
    NAME_2 %in% c("Mô", "Sotouboua") ~ "Sotouboua",
    NAME_2 %in% c("Agoe-Nyive", "Lomé", "Golfe") ~ "Golfe",
    TRUE ~ NAME_2  # Keep others as is
  )) %>%
  mutate(likelihood = factor(
    case_when(
      merged_prefecture %in% c("Assoli") ~ "Very Likely",
      merged_prefecture %in% c("Bas-Mono", "Danyi", "Akébou", "Bimah", "Zio", "Tandjouaré") ~ "Likely",
      merged_prefecture %in% c("Doufelgou", "Moyen-Mono", "Agou", "Blitta", "Avé", "Est-Mono", "Yoto", "Tchaudjo", "Anié") ~ "Possibly",
      merged_prefecture %in% c("Kpendjal", "Wawa", "Ogou", "Kozah", "Tchamba", "Haho", "Tône", "Kpélé") ~ "Unlikely",
      merged_prefecture %in% c("Bassar", "Dankpen", "Oti", "Kéran", "Sotouboua", "Amou", "Kloto") ~ "Very Unlikely",
      merged_prefecture %in% c("Golfe", "Lacs", "Vo") ~ "Non-endemic",
      TRUE ~ NA_character_
    ),
    levels = c("Very Likely", "Likely", "Possibly", "Unlikely", "Very Unlikely", "Non-endemic")
  ))

# Optional: add a placeholder to ensure "Very Unlikely" appears in the legend
very_unlikely_placeholder <- togo_sf_B[1, ]
very_unlikely_placeholder$merged_prefecture <- "Placeholder"
very_unlikely_placeholder$likelihood <- "Very Unlikely"
very_unlikely_placeholder$geometry <- st_sfc(st_point())
st_crs(very_unlikely_placeholder) <- st_crs(togo_sf_B)

togo_merged_sf_B <- rbind(togo_sf_B, very_unlikely_placeholder) %>%
  group_by(merged_prefecture, likelihood) %>%
  summarize(geometry = st_union(geometry), .groups = "drop")

###############################################################################
# 4. Calculate Centroids and Retain Prefecture Names
###############################################################################
togo_centroids_B <- togo_merged_sf_B %>%
  mutate(geometry = st_centroid(geometry)) %>%
  st_drop_geometry() %>%
  bind_cols(as.data.frame(st_coordinates(togo_merged_sf_B %>% st_centroid())) %>%
              rename(x_center = X, y_center = Y)) %>%
  rename(prefecture = merged_prefecture)

###############################################################################
# 5. Prepare EOT Data (prefecture_data)
###############################################################################
prefecture_data <- tribble(
  ~prefecture,     ~count_lt5, ~count_5_19, ~count_20_59, ~count_60_89, ~count_ge90,
  "Kpendjal",      1,          0,          2,           5,           1,
  "Oti",           3,          0,          2,           7,           17,
  "Tandjouaré",    0,          0,          0,           1,           4,
  "Tône",          1,          0,          0,           1,           3,
  "Assoli",        0,          0,          0,           0,           3,
  "Bassar",        2,          6,          0,           2,           10,
  "Bimah",         0,          0,          0,           1,           3,
  "Dankpen",       2,          3,          0,           1,           13,
  "Doufelgou",     0,          1,          0,           2,           2,
  "Kéran",         10,         3,          0,           0,           0,
  "Kozah",         1,          2,          0,           0,           24,
  "Blitta",        0,          0,          0,           10,          18,
  "Sotouboua",     9,          3,          0,           13,          10,
  "Tchamba",       1,          0,          0,           6,           9,
  "Tchaudjo",      0,          1,          0,           2,           5,
  "Agou",          0,          0,          0,           3,           10,
  "Akébou",        0,          0,          0,           1,           1,
  "Amou",          3,          0,          0,           3,           2,
  "Anié",          0,          0,          0,           3,           9,
  "Danyi",         0,          0,          0,           2,           3,
  "Est-Mono",      0,          0,          0,           5,           10,
  "Haho",          0,          3,          0,           2,           13,
  "Kloto",         3,          0,          0,           0,           2,
  "Kpélé",         1,          0,          0,           0,           1,
  "Moyen-Mono",    0,          0,          0,           4,           2,
  "Ogou",          0,          2,          0,           16,          19,
  "Wawa",          2,          0,          0,           4,           4,
  "Avé",           0,          0,          1,           0,           3,
  "Bas-Mono",      0,          0,          0,           2,           0,
  "Golfe",         0,          0,          0,           0,           0,
  "Lacs",          0,          0,          0,           0,           0,
  "Vo",            0,          0,          0,           0,           0,
  "Yoto",          0,          0,          0,           7,           15,
  "Zio",           0,          0,          0,           1,           11
)

# Convert counts to proportions (optional, so pies sum to 1)
prefecture_data <- prefecture_data %>%
  rowwise() %>%
  mutate(total_villages = sum(c_across(starts_with("count_")))) %>%
  mutate(
    p_lt5   = count_lt5   / total_villages,
    p_5_19  = count_5_19  / total_villages,
    p_20_59 = count_20_59 / total_villages,
    p_60_89 = count_60_89 / total_villages,
    p_ge90  = count_ge90  / total_villages
  ) %>%
  ungroup()

###############################################################################
# 6.1. Rename Proportion Columns for the Pie Legend
###############################################################################
prefecture_data <- prefecture_data %>%
  rename(
    `<5` = p_lt5,
    `5-19` = p_5_19,
    `20-59` = p_20_59,
    `60-89` = p_60_89,
    `≥90` = p_ge90
  )

###############################################################################
# 7. Create Centroids for Pies & Join EOT Data
###############################################################################
togo_centroids <- togo_merged_sf_B %>%
  st_centroid() %>%
  cbind(st_coordinates(.)) %>%
  as.data.frame() %>%
  rename(x_center = X, y_center = Y) %>%
  rename(prefecture = merged_prefecture)

plot_data <- togo_centroids %>%
  select(prefecture, x_center, y_center) %>%
  left_join(prefecture_data, by = "prefecture") %>%
  mutate(
    pie_x = x_center,
    pie_y = y_center
  )

# Ensure total_villages is in prefecture_data; if not, recalc:
plot_data <- plot_data %>%
  mutate(village_count_cat = case_when(
    .data$total_villages < 5 ~ "<5",
    .data$total_villages >= 5 & .data$total_villages < 20 ~ "5-19",
    .data$total_villages >= 20 & .data$total_villages < 60 ~ "≥20"
  )) %>%
  mutate(pie_radius = case_when(
    village_count_cat == "<5"   ~ 0.06,
    village_count_cat == "5-19"  ~ 0.09,
    village_count_cat == "≥20" ~ 0.12
  ))

###############################################################################
# 8. Create Dummy Data for the Size Legend (placed outside map extent)
###############################################################################
x_dummy <- min(plot_data$x_center) - 1
y_dummy <- min(plot_data$y_center) - 1
dummy_size <- data.frame(pie_radius = c(0.08, 0.12, 0.16),
                         x_dummy = x_dummy,
                         y_dummy = y_dummy)

###############################################################################
# 9. Create the Plot with Three Legends
###############################################################################
library(ggrepel)
down_offset   <- 0.026 
pie_radius_pt <- 8 
my_map <- ggplot() +
  # A) Base map: merged prefectures colored by likelihood
  geom_sf(
    data = togo_merged_sf_B, 
    aes(fill = likelihood), 
    color = "white", 
    linewidth = 1
  ) +
  scale_fill_manual(
    name = "Prefecture Elimination of Transmission (EOT)",
    values = c(
      "Very Likely"   = "#1a9850",
      "Likely"        = "#a6d96a",
      "Possibly"      = "#fee08b",
      "Unlikely"      = "#fdae61",
      "Very Unlikely" = "#f46d43",
      "Non-endemic"   = "grey"
    )
  ) +
  # B) Regions border (Level 1)
  geom_sf(
    data = togo_regions_sf, 
    fill = NA, 
    color = "black", 
    linewidth = 1.2
  ) +
  # Insert a new fill scale for the pie colors
  new_scale_fill() +
  # C) Segments connecting centroids to pie positions
  geom_segment(
    data = plot_data,
    aes(x = x_center, y = y_center, xend = pie_x, yend = pie_y),
    color = "black", linetype = "dotted"
  ) +
  # D) Pies: show EOT proportions per prefecture with size based on village sample size
  geom_scatterpie(
    data = plot_data,
    aes(x = pie_x, y = pie_y, r = pie_radius, group = prefecture),
    cols = c("<5", "5-19", "20-59", "60-89", "≥90"),
    color = "black", alpha = 0.8
  ) +
  scale_fill_manual(
    name = "Village Proportion of EOT Probabilities (%)",
    values = c(
      "<5"   = "red",
      "5-19"  = "orange",
      "20-59" = "yellow",
      "60-89" = "lightgreen",
      "≥90"  = "darkgreen"
    )
  ) +
  # E) Dummy layer for the size legend
  geom_point(
    data = dummy_size,
    aes(x = x_dummy, y = y_dummy, size = pie_radius),
    alpha = 0, 
    show.legend = TRUE
  ) +
  scale_size_continuous(
    name = "Number of Villages Surveyed",
    breaks = c(0.08, 0.12, 0.16),
    labels = c("<5", "5-19", "≥20")
  ) +
  guides(size = guide_legend(override.aes = list(size = c(5,7.5,10.5),alpha = 1, fill = "grey", color = "black", shape = 21))) +
  # F) Label each merged prefecture using centroids
  geom_text_repel(
    data = togo_centroids_B,
    aes(x = x_center, y = y_center - down_offset, label = prefecture),
    size = 3, colour = "black", fontface = "bold",
    bg.color = "white", bg.r = 0.15,
    direction = "y",          # only move vertically
    nudge_y   = -down_offset, # start below the point
    vjust     = 0.7,            # anchor text baseline above
    point.padding = unit(pie_radius_pt, "pt"),
    box.padding   = 0.15,
    min.segment.length = 0,
    segment.size = 0.25,
    seed = 5,
    force = 2, force_pull = 0  # bias away from the point
  )+
  theme_void() +
  labs(title = "", fill = NULL) +
  coord_sf() +
  theme(
    plot.title = element_text(hjust = -0.1, vjust = 0.5, size = 14, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8)
  )

###############################################################################
# 10. Display the Plot
###############################################################################
print(my_map)

###############################################################################
# 11. Save the Plot (Optional)
###############################################################################
ggsave("Figure 7.pdf",device = cairo_pdf, plot = my_map, width = 16, height = 8, dpi = 300)

####ALSO CREATE MAKE WITH COLOUR BLIND SAFE PATTERNS
# 9. Create the Plot with Three Legends
###############################################################################
library(ggrepel)
down_offset   <- 0.026 
pie_radius_pt <- 8 
my_map <- ggplot() +
  # A) Base map: merged prefectures colored by likelihood
  geom_sf(
    data = togo_merged_sf_B, 
    aes(fill = likelihood), 
    color ="black", 
    linewidth = 1
  ) +
  scale_fill_manual(
    name = "Prefecture Elimination of Transmission (EOT)",
    values = c(
      "Very Likely"   = "#08519c",
      "Likely"        = "#3182bd",
      "Possibly"      = "#6baed6",
      "Unlikely"      = "#c6dbef",
      "Very Unlikely" = "#f7fbff",
      "Non-endemic"   = "grey"
    )
  ) +
  # B) Regions border (Level 1)
  geom_sf(
    data = togo_regions_sf, 
    fill = NA, 
    color = "black", 
    linewidth = 1.5
  ) +
  # Insert a new fill scale for the pie colors
  new_scale_fill() +
  # C) Segments connecting centroids to pie positions
  geom_segment(
    data = plot_data,
    aes(x = x_center, y = y_center, xend = pie_x, yend = pie_y),
    color = "black", linetype = "dotted"
  ) +
  # D) Pies: show EOT proportions per prefecture with size based on village sample size
  geom_scatterpie(
    data = plot_data,
    aes(x = pie_x, y = pie_y, r = pie_radius, group = prefecture),
    cols = c("<5", "5-19", "20-59", "60-89", "≥90"),
    color = "black", alpha = 0.8
  ) +
  scale_fill_manual(
    name = "Village Proportion of EOT Probabilities (%)",
    values = c(
      "<5"   = "#f7fbff",
      "5-19"  = "#c6dbef",
      "20-59" = "#6baed6",
      "60-89" = "#3182bd",
      "≥90"  = "#08519c"
    )
  ) +
  # E) Dummy layer for the size legend
  geom_point(
    data = dummy_size,
    aes(x = x_dummy, y = y_dummy, size = pie_radius),
    alpha = 0, 
    show.legend = TRUE
  ) +
  scale_size_continuous(
    name = "Number of Villages Surveyed",
    breaks = c(0.08, 0.12, 0.16),
    labels = c("<5", "5-19", "≥20")
  ) +
  guides(size = guide_legend(override.aes = list(size = c(5,7.5,10.5),alpha = 1, fill = "grey", color = "black", shape = 21))) +
  # F) Label each merged prefecture using centroids
  geom_text_repel(
    data = togo_centroids_B,
    aes(x = x_center, y = y_center - down_offset, label = prefecture),
    size = 3, colour = "black", fontface = "bold",
    bg.color = "white", bg.r = 0.15,
    direction = "y",          # only move vertically
    nudge_y   = -down_offset, # start below the point
    vjust     = 0.7,            # anchor text baseline above
    point.padding = unit(pie_radius_pt, "pt"),
    box.padding   = 0.15,
    min.segment.length = 0,
    segment.size = 0.25,
    seed = 5,
    force = 2, force_pull = 0  # bias away from the point
  )+
  theme_void() +
  labs(title = "", fill = NULL) +
  coord_sf() +
  theme(
    plot.title = element_text(hjust = -0.1, vjust = 0.5, size = 14, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8)
  )

###############################################################################
# 10. Display the Plot
###############################################################################
print(my_map)

###############################################################################
# 11. Save the Plot (Optional)
###############################################################################
ggsave("Figure 7.pdf", device = cairo_pdf, plot = my_map, width = 16, height = 8, dpi = 300)

