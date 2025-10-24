#Figure 1 of Supplementary Material
#figure of baseline microfilarial endemicity
# Load necessary packages
library(geodata)
library(sf)
library(ggplot2)
library(dplyr)
library(patchwork)

# 1. Download Togo administrative boundaries (level 2 for prefectures, level 1 for regions)
togo_prefectures <- gadm(country = "TGO", level = 2, path = tempdir())
togo_regions <- gadm(country = "TGO", level = 1, path = tempdir())

# 2. Convert both SpatVector objects to sf (simple features)
togo_sf <- st_as_sf(togo_prefectures)
togo_regions_sf <- st_as_sf(togo_regions)

# 3. Create Plot A -------------------------------------------------------------

# Define merging for prefectures in Plot A
togo_sf_A <- togo_sf %>%
  mutate(merged_prefecture = case_when(
    NAME_2 %in% c("Cinkassé", "Tône") ~ "Tône",  # Merge Cinkassé into Tône
    NAME_2 %in% c("Naki-Ouest", "Kpendjal") ~ "Kpendjal",  # Merge Naki-Ouest into Kpendjal
    NAME_2 %in% c("Oti", "Oti-Sud") ~ "Oti",  # Merge Oti-Sud into Oti
    NAME_2 %in% c("Mô", "Sotouboua") ~ "Sotouboua",  # Merge Mô into Sotouboua
    NAME_2 %in% c("Agoe-Nyive", "Lomé", "Golfe") ~ "Golfe",  # Merge Agoe-Nyive, Lomé, Golfe into Golfe
    TRUE ~ NAME_2  # Keep others as they are
  ))

# Assign baseline prevalence classifications for Plot A
togo_sf_A <- togo_sf_A %>%
  mutate(prevalence = factor(case_when(
    merged_prefecture %in% c("Bas-Mono", "Golfe", "Lacs", "Vo") ~ "Non-endemic",
    merged_prefecture %in% c( "Dankpen", "Kloto", "Avé", "Zio","Tandjouaré") ~ "Hypoendemic",
    merged_prefecture %in% c("Assoli", "Oti", "Tchamba", "Danyi", "Kpélé") ~ "Mesoendemic",
    merged_prefecture %in% c("Kpendjal","Tône", "Bassar", "Bimah", "Doufelgou", "Kozah", "Blitta", "Yoto", "Sotouboua",
                             "Tchaudjo", "Agou", "Akébou", "Amou", "Anié", "Est-Mono", "Haho", "Moyen-Mono", 
                             "Ogou", "Wawa") ~ "Hyperendemic",
    merged_prefecture %in% c("Kéran") ~ "Holoendemic",
    TRUE ~ NA_character_
  ), levels = c("Non-endemic", "Hypoendemic", "Mesoendemic", "Hyperendemic", "Holoendemic")))  # Specify levels

# Placeholder for Plot A (in case it's needed for mapping aesthetics)
very_unlikely_placeholder_A <- togo_sf_A[1, ]
very_unlikely_placeholder_A$merged_prefecture <- "Placeholder"
very_unlikely_placeholder_A$prevalence <- "Holoendemic"
very_unlikely_placeholder_A$geometry <- st_sfc(st_point())  # Empty geometry

# Ensure CRS (Coordinate Reference System) is consistent
st_crs(very_unlikely_placeholder_A) <- st_crs(togo_sf_A)

# Merge placeholder with real data
togo_merged_sf_A <- rbind(togo_sf_A, very_unlikely_placeholder_A) %>%
  group_by(merged_prefecture, prevalence) %>%
  summarize(geometry = st_union(geometry), .groups = "drop")  # Merge polygons

# Calculate centroids for labels in Plot A
togo_centroids_A <- st_centroid(togo_merged_sf_A)

# Define colors for Plot A (baseline prevalence)
prevalence_colors_A <- c(
  "Non-endemic" = "#a6d96a",  # Light Green
  "Hypoendemic" = "#d9ef8b",  # Yellow
  "Mesoendemic" = "#fee08b",  # Orange
  "Hyperendemic" = "#fdae61",  # Red
  "Holoendemic" = "#d73027"  # Dark Red
)

# Create Plot A
plot_A <- ggplot() +
  geom_sf(data = togo_merged_sf_A, aes(fill = prevalence), color = "white", size = 0.2) +  # Prefecture layer
  geom_sf(data = togo_regions_sf, fill = NA, color = "black", linewidth = 1.2) +  # Region borders
  geom_text(data = togo_centroids_A, aes(geometry = geometry, label = merged_prefecture), 
            stat = "sf_coordinates", size = 3, color = "black", fontface = "bold") +  # Prefecture names
  scale_fill_manual(values = prevalence_colors_A, na.value = "white", drop = FALSE) +  # Colors for prevalence
  theme_void() +  # Clean theme
  labs(title = "B", fill = "Baseline endemicity") +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8)
  ) +
  coord_sf()

# 4. Create Plot B -------------------------------------------------------------

# Define merging for prefectures in Plot B
togo_sf_B <- togo_sf %>%
  mutate(merged_prefecture = case_when(
    NAME_2 %in% c("Cinkassé", "Tône") ~ "Tône",  # Merge Cinkassé into Tône
    NAME_2 %in% c("Naki-Ouest", "Kpendjal") ~ "Kpendjal",  # Merge Naki-Ouest into Kpendjal
    NAME_2 %in% c("Oti", "Oti-Sud") ~ "Oti",  # Merge Oti-Sud into Oti
    NAME_2 %in% c("Mô", "Sotouboua") ~ "Sotouboua",  # Merge Mô into Sotouboua
    NAME_2 %in% c("Agoe-Nyive", "Lomé", "Golfe") ~ "Golfe",  # Merge Agoe-Nyive, Lomé, Golfe into Golfe
    TRUE ~ NAME_2  # Keep others as they are
  ))

# Assign baseline prevalence classifications for Plot B
togo_sf_B <- togo_sf_B %>%
  mutate(prevalence = factor(case_when(
    merged_prefecture %in% c("Golfe", "Lacs", "Vo") ~ "Non-endemic",
    merged_prefecture %in% c("Tandjouaré") ~ "Hypoendemic",
    merged_prefecture %in% c("Assoli","Bas-Mono", "Zio") ~ "Mesoendemic",
    merged_prefecture %in% c("Kpendjal","Tône", "Bimah", "Doufelgou", "Blitta", "Yoto", "Tchaudjo", "Agou", "Akébou", 
                             "Amou", "Anié", "Est-Mono", "Haho", "Moyen-Mono", "Ogou", "Wawa", "Tchamba", 
                             "Danyi","Oti", "Kloto", "Kpélé", "Avé") ~ "Hyperendemic",
    merged_prefecture %in% c("Kéran", "Bassar", "Dankpen", "Kozah", "Sotouboua") ~ "Holoendemic",
    TRUE ~ NA_character_
  ), levels = c("Non-endemic", "Hypoendemic", "Mesoendemic", "Hyperendemic", "Holoendemic")))
# Placeholder for Plot B
very_unlikely_placeholder_B <- togo_sf_B[1, ]
very_unlikely_placeholder_B$merged_prefecture <- "Placeholder"
very_unlikely_placeholder_B$prevalence <- "Holoendemic"
very_unlikely_placeholder_B$geometry <- st_sfc(st_point())  # Empty geometry

# Ensure CRS (Coordinate Reference System) is consistent for Plot B
st_crs(very_unlikely_placeholder_B) <- st_crs(togo_sf_B)

# Combine placeholder with real data for Plot B
togo_merged_sf_B <- rbind(togo_sf_B, very_unlikely_placeholder_B) %>%
  group_by(merged_prefecture, prevalence) %>%
  summarize(geometry = st_union(geometry), .groups = "drop")  # Merge polygons

# Calculate centroids for labels in Plot B
togo_centroids_B <- st_centroid(togo_merged_sf_B)

# Define colors for Plot B (similar to Plot A)
prevalence_colors_B <- c(
  "Non-endemic" = "#a6d96a",  # Light Green
  "Hypoendemic" = "#d9ef8b",  # Yellow
  "Mesoendemic" = "#fee08b",  # Orange
  "Hyperendemic" = "#fdae61",  # Red
  "Holoendemic" = "#d73027"  # Dark Red
)

# Create Plot B
plot_B <- ggplot() +
  geom_sf(data = togo_merged_sf_B, aes(fill = prevalence), color = "white", size = 0.2) +  # Prefecture layer
  geom_sf(data = togo_regions_sf, fill = NA, color = "black", linewidth = 1.2) +  # Region borders
  geom_text(data = togo_centroids_B, aes(geometry = geometry, label = merged_prefecture), 
            stat = "sf_coordinates", size = 3, color = "black", fontface = "bold") +  # Prefecture names
  scale_fill_manual(values = prevalence_colors_B, na.value = "white", drop = FALSE) +  # Colors for prevalence
  theme_void() +  # Clean theme for publication
  labs(title = "C", fill = "Baseline endemicity") +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8)
  ) +
  coord_sf()

# 5. Combine Plot A and Plot B -------------------------------------------------

# Use patchwork to display the plots side by side
combined_plot <- plot_A + plot_B + plot_layout(ncol = 2)

# Display the combined plot
print(combined_plot)

# Optionally, save the combined plot to a file
#ggsave("togo_combined_baseline_prevalence_A_B.png", plot = combined_plot, width = 16, height = 8, dpi = 300)

# Assume togo_sf and togo_regions_sf are already created as in previous steps

# 1. Create Plot C -------------------------------------------------------------

# Define merging for prefectures in Plot C and assign prevalence levels
togo_sf_C <- togo_sf %>%
  mutate(merged_prefecture = case_when(
    NAME_2 %in% c("Cinkassé", "Tône") ~ "Tône",  # Merge Cinkassé into Tône
    NAME_2 %in% c("Naki-Ouest", "Kpendjal") ~ "Kpendjal",  # Merge Naki-Ouest into Kpendjal
    NAME_2 %in% c("Oti", "Oti-Sud") ~ "Oti",  # Merge Oti-Sud into Oti
    NAME_2 %in% c("Mô", "Sotouboua") ~ "Sotouboua",  # Merge Mô into Sotouboua
    NAME_2 %in% c("Agoe-Nyive", "Lomé", "Golfe") ~ "Golfe",  # Merge Agoe-Nyive, Lomé, Golfe into Golfe
    TRUE ~ NAME_2  # Keep others as they are
  ))

# Assign prevalence classifications for Plot C (Non-endemic to Hypoendemic, Mesoendemic, Hyper- to Holoendemic)
togo_sf_C <- togo_sf_C %>%
  mutate(prevalence = factor(case_when(
    merged_prefecture %in% c("Golfe", "Lacs", "Vo") ~ "Non-endemic to Hypoendemic",  # Combined non-endemic and hypoendemic
    merged_prefecture %in% c("Assoli", "Bas-Mono", "Zio", "Danyi", "Kloto", "Kpélé", "Agou", "Yoto", "Avé") ~ "Mesoendemic",
    merged_prefecture %in% c("Kpendjal", "Tône", "Tandjouaré", "Bimah", "Doufelgou", "Blitta", "Tchaudjo", "Akébou", 
                             "Amou", "Anié", "Est-Mono", "Haho", "Moyen-Mono", "Ogou", "Wawa", "Tchamba", "Kéran", 
                             "Oti", "Bassar", "Dankpen", "Kozah", "Sotouboua") ~ "Hyper- to Holoendemic",
    TRUE ~ NA_character_
  ), levels = c("Non-endemic to Hypoendemic", "Mesoendemic", "Hyper- to Holoendemic")))  # Specify levels

# Merge polygons for Plot C
togo_merged_sf_C <- togo_sf_C %>%
  group_by(merged_prefecture, prevalence) %>%
  summarize(geometry = st_union(geometry), .groups = "drop")

# Calculate centroids for Plot C labels
togo_centroids_C <- st_centroid(togo_merged_sf_C)

# Define colors for Plot C
prevalence_colors_C <- c(
  "Non-endemic to Hypoendemic" = "#a6d96a",  # Light Green
  "Mesoendemic" = "#fee08b",  # Orange
  "Hyper- to Holoendemic" = "#fdae61"  # Dark Red
)

# Create Plot C
plot_C <- ggplot() +
  geom_sf(data = togo_merged_sf_C, aes(fill = prevalence), color = "white", size = 0.2) +  # Prefecture layer
  geom_sf(data = togo_regions_sf, fill = NA, color = "black", linewidth = 1.2) +  # Region borders
  geom_text(data = togo_centroids_C, aes(geometry = geometry, label = merged_prefecture), 
            stat = "sf_coordinates", size = 3, color = "black", fontface = "bold") +  # Prefecture names
  scale_fill_manual(values = prevalence_colors_C, na.value = "white", drop = FALSE) +  # Colors for prevalence
  theme_void() +  # Clean theme
  labs(title = "A", fill = "Baseline endemicity") +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8)
  ) +
  coord_sf()

# 2. Create Plot D -------------------------------------------------------------

# Define merging for prefectures in Plot D and assign prevalence levels
togo_sf_D <- togo_sf %>%
  mutate(merged_prefecture = case_when(
    NAME_2 %in% c("Cinkassé", "Tône") ~ "Tône",  # Merge Cinkassé into Tône
    NAME_2 %in% c("Naki-Ouest", "Kpendjal") ~ "Kpendjal",  # Merge Naki-Ouest into Kpendjal
    NAME_2 %in% c("Oti", "Oti-Sud") ~ "Oti",  # Merge Oti-Sud into Oti
    NAME_2 %in% c("Mô", "Sotouboua") ~ "Sotouboua",  # Merge Mô into Sotouboua
    NAME_2 %in% c("Agoe-Nyive", "Lomé", "Golfe") ~ "Golfe",  # Merge Agoe-Nyive, Lomé, Golfe into Golfe
    TRUE ~ NAME_2  # Keep others as they are
  ))

# Assign prevalence classifications for Plot D (Original levels)
togo_sf_D <- togo_sf_D %>%
  mutate(prevalence = factor(case_when(
    merged_prefecture %in% c("Golfe", "Lacs","Vo") ~ "Non-endemic to Hypoendemic",
    merged_prefecture %in% c("Assoli","Bas-Mono", "Zio") ~ "Mesoendemic",
    merged_prefecture %in% c("Tône", "Tandjouaré", "Kpendjal", "Bimah", "Doufelgou", "Blitta", "Yoto", "Tchaudjo", 
                             "Agou", "Akébou", "Amou", "Anié", "Est-Mono", "Haho", "Moyen-Mono", "Ogou", "Wawa", 
                             "Tchamba","Oti", "Danyi", "Kloto", "Kpélé", "Avé") ~ "Hyperendemic",
    merged_prefecture %in% c("Kéran", "Bassar", "Dankpen", "Kozah", "Sotouboua") ~ "Holoendemic",
    TRUE ~ NA_character_
  ), levels = c("Non-endemic to Hypoendemic", "Mesoendemic", "Hyperendemic", "Holoendemic")))

# Merge polygons for Plot D
togo_merged_sf_D <- togo_sf_D %>%
  group_by(merged_prefecture, prevalence) %>%
  summarize(geometry = st_union(geometry), .groups = "drop")

# Calculate centroids for Plot D labels
togo_centroids_D <- st_centroid(togo_merged_sf_D)

# Define colors for Plot D (same as previous plots)
prevalence_colors_D <- c(
  "Non-endemic to Hypoendemic" = "#a6d96a",  # Light Green
  "Mesoendemic" = "#fee08b",  # Orange
  "Hyperendemic" = "#fdae61",  # Red
  "Holoendemic" = "#d73027"  # Dark Red
)

# Create Plot D
plot_D <- ggplot() +
  geom_sf(data = togo_merged_sf_D, aes(fill = prevalence), color = "white", size = 0.2) +  # Prefecture layer
  geom_sf(data = togo_regions_sf, fill = NA, color = "black", linewidth = 1.2) +  # Region borders
  geom_text(data = togo_centroids_D, aes(geometry = geometry, label = merged_prefecture), 
            stat = "sf_coordinates", size = 3, color = "black", fontface = "bold") +  # Prefecture names
  scale_fill_manual(values = prevalence_colors_D, na.value = "white", drop = FALSE) +  # Colors for prevalence
  theme_void() +  # Clean theme
  labs(title = "D", fill = "Baseline endemicity") +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8)
  ) +
  coord_sf()

# 3. Combine Plots A, B, C, and D ---------------------------------------------

# Use patchwork to display the plots in a 2x2 grid
combined_plotI <- plot_C + plot_A + plot_layout(ncol = 2)
combined_plotF <- plot_B + plot_D + plot_layout(ncol = 2)

# Display the combined plot
print(combined_plotI)
print(combined_plotF)
combined_plots <- (plot_C + plot_A) / (plot_B + plot_D)


#Figure 2 of Supplementary Material
#Now figure of OCP and SIZ
library(geodata)
library(sf)
library(ggplot2)
library(dplyr)
library(patchwork)

# 1. Download Togo administrative boundaries (prefectures and regions)
togo_prefectures <- gadm(country = "TGO", level = 2, path = tempdir())
togo_regions <- gadm(country = "TGO", level = 1, path = tempdir())

# Convert both SpatVector objects to sf (simple features)
togo_sf <- st_as_sf(togo_prefectures)
togo_regions_sf <- st_as_sf(togo_regions)

# 2. Map 1: Onchocerciasis Control Programme Areas

# Assign the control programme phase to each prefecture
togo_sf_ocp <- togo_sf %>%
  mutate(ocp_phase = case_when(
    NAME_2 %in% c("Cinkassé", "Tône", "Tandjouaré") ~ "Phase II",
    NAME_2 %in% c("Kpendjal", "Oti", "Naki-Ouest","Oti-Sud", "Kéran", "Dankpen", "Doufelgou", "Bimah", 
                  "Kozah", "Assoli", "Bassar") ~ "Phase III East",
    NAME_2 %in% c("Golfe", "Lomé", "Lacs", "Agoe-Nyive") ~ "Not Included",  # Set to grey without adding a new category
    TRUE ~ "Southern Extension"  # Remaining areas, 
  ))

# 3. Create the map for Onchocerciasis Control Programme Areas

# Define colors for the control programme areas
ocp_colors <- c(
  "Phase II" = "#fdae61",           # Muted Orange
  "Phase III East" = "#fee08b",     # Muted Yellow
  "Southern Extension" = "#d73027", # Dark Red
  "Not Included" = "grey70"                 # Grey for Mô, Sotouboua, and Tchaudjo
)

# Calculate centroids for labeling
togo_centroids_ocp <- st_centroid(togo_sf_ocp)

# Plot the Onchocerciasis Control Programme areas
plot_ocp <- ggplot() +
  geom_sf(data = filter(togo_sf_ocp, ocp_phase != "Grey"), aes(fill = ocp_phase), color = "white", size = 0.2) +  # Prefecture layer excluding grey regions from the legend
  geom_sf(data = filter(togo_sf_ocp, ocp_phase == "Grey"), fill = "#d73027", color = "white", size = 0.2) +  # Grey regions without legend
  geom_sf(data = togo_regions_sf, fill = NA, color = "black", linewidth = 1.2) +  # Region borders
  geom_text(data = togo_centroids_ocp, aes(geometry = geometry, label = NAME_2), 
            stat = "sf_coordinates", size = 3, color = "black", fontface = "bold") +  # Prefecture names
  scale_fill_manual(values = ocp_colors, na.value = "white", drop = TRUE,  breaks = c("Phase II", "Phase III East","Southern Extension","Not Included")) +  # Colors for OCP phases
  theme_void() +  # Clean theme
  labs(title = "A", fill = "Onchocerciasis Control Programme in West Africa (OCP) Phase") +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8)
  ) +
  coord_sf()

# 4. Map 2: Special Intervention Zones (SIZ)

# Assign Yes/No for Special Intervention Zones (SIZ), including Oti-Sud
togo_sf_siz <- togo_sf %>%
  mutate(siz_status = case_when(
    NAME_2 %in% c("Kpendjal", "Oti", "Oti-Sud", "Kéran", 
                  "Dankpen", "Doufelgou", "Bimah", "Kozah", "Assoli", "Bassar", 
                  "Tchaudjo", "Sotouboua", "Mô") ~ "Yes",
    TRUE ~ "No"
  ))

# Define colors for SIZ map (muted yellow and white)
siz_colors <- c(
  "Yes" = "#e0e0b3",  # Soft, muted yellow for Yes
  "No" = "#f7f7f7"  # Light grey for No
)

# Calculate centroids for SIZ labels
togo_centroids_siz <- st_centroid(togo_sf_siz)

# Plot Special Intervention Zones (SIZ)
plot_siz <- ggplot() +
  geom_sf(data = togo_sf_siz, aes(fill = siz_status), color = "white", size = 0.2) +  # Prefecture layer
  geom_sf(data = togo_regions_sf, fill = NA, color = "black", linewidth = 1.2) +  # Region borders
  geom_text(data = togo_centroids_siz, aes(geometry = geometry, label = NAME_2), 
            stat = "sf_coordinates", size = 3, color = "black", fontface = "bold") +  # Prefecture names
  scale_fill_manual(values = siz_colors, na.value = "white", drop = FALSE, 
                    breaks = c("Yes", "No")) +  # Order "Yes" first in the legend
  theme_void() +  # Clean theme
  labs(title = "B", fill = "Special Intervention Zones (SIZ)") +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8)
  ) +
  coord_sf()

# 5. Combine the two maps ------------------------------------------------------

# Use patchwork to combine the two maps side by side
combined_maps <- plot_ocp + plot_siz + plot_layout(ncol = 2)
combined_maps



#Figure 3 of Supplementary Material
#11.4 plot of crude prevalence of SIZ and non-SIZ zones over time
#For the below plot, it makes more sense to use the later created jitter E8 database, to see all the points.
library(ggplot2)
library(tidyr)
# Replace NA values in the SIZ column with "No"
E8$SIZ <- replace_na(E8$SIZ, "No")
# Create the plot with strict axis limits and remove padding
# Create the plot with vertical lines and labels for key events
E8$SIZ <- factor(E8$SIZ, levels = c("Yes", "No"))
c1 <- ggplot(E8, aes(x = Survey_year, y = PREV_cr, colour = SIZ)) +
  geom_point(size = 3.5) +  # Increase point size
  scale_y_continuous(
    name = "Crude Microfilarial Prevalence (%)", 
    breaks = seq(0, 100, by = 10), 
    limits = c(0, 100), 
    expand = c(0, 0)  # Remove padding on the y-axis
  ) +
  scale_x_continuous(
    name = "Year", 
    breaks = seq(1970, 2020, by = 5), 
    limits = c(1975, 2020), 
    expand = c(0, 0)  # Remove padding on the x-axis
  ) +
  scale_color_discrete(name = "Village within SIZ") +
  # Add vertical lines for key events
  geom_vline(xintercept = 2002, linetype = "dashed", color = "blue", size = 1) +
  geom_vline(xintercept = 2007, linetype = "dashed", color = "red", size = 1) +
  
  # Annotate the key events
  annotate("text", x = 2002, y = 98, label = "OCP Closure & SIZ Start", color = "blue", angle = 90, vjust = -0.5, hjust = 1, size = 4.5) +
  annotate("text", x = 2007, y = 98, label = "End of Vector Control in SIZ", color = "red", angle = 90, vjust = -0.5, hjust = 1, size = 4.5) +
  
  # Add theme and styles
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),  # Center the title
    axis.text.y = element_text(size = 17),
    axis.text.x = element_text(size = 17),
    axis.title = element_text(size = 17),
    legend.text = element_text(size = 15),  # Increase SIZ legend text size
    legend.title = element_text(size = 17)  # Increase SIZ legend title size
  )

# Display the plot
c1

#Figure 4 of Supplementary Material 
library(ggplot2)
library(dplyr)
library(patchwork)  # For arranging plots side-by-side

# Data subsets for each period
# Replace NA values in the SIZ column with "No"
#E3$SIZ <- replace_na(E3$SIZ, "No")
subset_1998_2002 <- E8 %>% filter(Survey_year >= 1998, Survey_year <= 2002)
subset_2007_2011 <- E8 %>% filter(Survey_year >= 2007, Survey_year <= 2011)
subset_2012_2017 <- E8 %>% filter(Survey_year >= 2012, Survey_year <= 2017)

# Create boxplot for 1998–2002 (Panel A)
plot_A <- ggplot(subset_1998_2002, aes(x = SIZ, y = PREV_cr, fill = SIZ)) +
  geom_boxplot(outlier.size = 3, color = "black") +
  scale_y_continuous(name = "Crude Microfilarial Prevalence (%)", breaks = seq(0, 100, by = 10), limits = c(0, 60)) +
  scale_fill_manual(values = c("No" = "lightblue", "Yes" = "salmon")) +
  labs(x = "SIZ status", subtitle = "A") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),axis.title = element_text(size = 15),
        plot.subtitle = element_text(size = 15, face = "bold"))

# Create boxplot for 2007–2011 (Panel B)
plot_B <- ggplot(subset_2007_2011, aes(x = SIZ, y = PREV_cr, fill = SIZ)) +
  geom_boxplot(outlier.size = 3, color = "black") +
  scale_y_continuous(name = "Crude Microfilarial Prevalence (%)", breaks = seq(0, 100, by = 10), limits = c(0, 60)) +
  scale_fill_manual(values = c("No" = "lightblue", "Yes" = "salmon")) +
  labs(x = "SIZ status", subtitle = "B") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),axis.title = element_text(size = 15),
        plot.subtitle = element_text(size = 15, face = "bold"))

# Create boxplot for 2012–2017 (Panel C)
plot_C <- ggplot(subset_2012_2017, aes(x = SIZ, y = PREV_cr, fill = SIZ)) +
  geom_boxplot(outlier.size = 3, color = "black") +
  scale_y_continuous(name = "Crude Microfilarial Prevalence (%)", breaks = seq(0, 100, by = 10), limits = c(0, 60)) +
  scale_fill_manual(values = c("No" = "lightblue", "Yes" = "salmon")) +
  labs(x = "SIZ status", subtitle = "C") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),axis.title = element_text(size = 15),
        plot.subtitle = element_text(size = 15, face = "bold"))

# Combine the three plots side-by-side
combined_plot <- plot_A + plot_B + plot_C + plot_layout(ncol = 3)

# Display the combined plot
print(combined_plot)


#Figure 5 of Supplementary Material
#Where it stated "dat", dataset E8 can also be used
#### 12 Map with Longitude and latitude
#sf and mapview packages
library(sf)
library(mapview)
#CORRECT GBALIB MISSING VALUES
dat$GBALIB[dat$BASNOM=="NYIVE"] <- "TODZIE"
dat$GBALIB[dat$BASNOM=="OTI"] <- "OTI"
dat$GBALIB[dat$BASNOM=="HAHO" | dat$BASNOM=="MENOU"] <- "VOLTA LAC EAST"
dat$GBALIB[dat$BASNOM=="ASUKAWKAW" | dat$BASNOM=="ASSOUKOKO"] <- "VOLTA LAC"
dat$GBALIB[dat$VILNOM=="LOUGOU"] <- "WHITE VOLTA"
dat$GBALIB[dat$VILNOM=="LAKATA-KONDJI"] <- "VOLTA LAC EAST"
dat$GBALIB[dat$BASNOM=="BIANKOURI/WHITE VOLTA"] <- "WHITE VOLTA"
dat$GBALIB[dat$BASNOM=="AMOU/MONO" |dat$BASNOM=="AFAN/MONO"] <- "MONO"
dat$GBALIB[dat$VILNOM=="KOUKOUMBOU"] <- "OTI"
dat$GBALIB[dat$VILNOM=="Tsokple"] <- "MONO"
dat$GBALIB[dat$VILNOM=="Tchalanide "|dat$VILNOM=="Sogbadai "|dat$VILNOM=="Bouzalo "] <- "VOLTA LAC"
dat$GBALIB[dat$BASNOM=="MO"|dat$BASNOM=="KATCHA/MO"] <- "VOLTA LAC"
dat$GBALIB[dat$VILNOM=="Kpati Cope" | dat$VILNOM=="Atinkpassa " |dat$VILNOM=="Amouta "] <- "MONO"
dat$GBALIB[dat$VILNOM=="AGBARADA" | dat$VILNOM=="ABOUDA"] <- "OTI"
#now plot map
dat1 <- data.frame(dat$GBALIB,dat$LON,dat$LAT)
dat2 <- na.omit(dat1)
mmap <- st_as_sf(dat2, coords = c("dat.LON", "dat.LAT"), crs = 4326)
mapview(mmap)
#Although the coordinates are not accurate, and many points are outside Togo when all of them are in the country, the district data is correct (important for the model), I checked many of the values manually.
#correct some lat/lon
dat$LON[dat$VILNOM=="AGLAGO-KOPE"] <- "0.77"  #https://www.azpostcodes.com/tgo/place-plateaux-agou-aglagokope/
dat$LAT[dat$VILNOM=="AGLAGO-KOPE"] <- "6.76"  #https://www.azpostcodes.com/tgo/place-plateaux-agou-aglagokope/
dat$LON[dat$VILNOM=="FAZAO"] <- "0.77"  #adjusted accoridng to the values od the dataset
dat$LON[dat$VILNOM=="BANDA"] <- "0.375"  #BASED ON BOULOHOU
dat$LAT[dat$VILNOM=="BANDA"] <- "8.375"  #BASED ON BOULOHOU
dat$LAT[dat$VILNOM=="KATCHALIKADI"] <- "8.7"  #BASED ON THE VILLAGE OTHER VALUES
dat$LON[dat$VILNOM=="KATCHALIKADI"] <- "0.77"  #BASED ON THE VILLAGE OTHER VALUES
dat$LON[dat$VILNOM=="KPABTE"] <- "0.31"  #BASED ON THE VILLAGE OTHER VALUES
dat$LON[dat$VILNOM=="FAWUKPE"] <- "0.55"  #BASED ON OTHER VILLAGES IN DEVEHO
dat$LON[dat$VILNOM=="AMMOUTA"] <- "0.72"  #BASED ON THE VILLAGE OTHER VALUES
dat$LAT[dat$VILNOM=="AMMOUTA"] <- "7.33"  #BASED ON THE VILLAGE OTHER VALUES
dat$LAT[dat$VILNOM=="WASSITE"] <- "10.00"  #BASED ON THE VILLAGE FROM NABOULOUGOU
dat$LAT[dat$VILNOM=="KEDJEBOUA"] <- "8.4"  #BASED ON KEDJEBI-LOOU
dat$LAT[dat$VILNOM=="KPAGBAZIBIYO"] <- "9.6"  #BASED ON kOZAH PREFECTURE

dat1 <- data.frame(dat$REGNOM,dat$LON,dat$LAT)
dat2 <- na.omit(dat1)
Togo <- st_as_sf(dat2, coords = c("dat.LON", "dat.LAT"), crs = 4326)
mapview(Togo) #map of the districts
# Rename the REGNOM column to "Regions of Togo"
names(dat2)[names(dat2) == "dat.REGNOM"] <- "Regions"
# Convert the dataframe to a spatial dataframe
Togo <- st_as_sf(dat2, coords = c("dat.LON", "dat.LAT"), crs = 4326)
# Create the map with the renamed column
mapview(Togo, zcol = "Regions")
#highlight the country's borders
# Install necessary packages if not already installed
library(geodata)
library(sf)
library(mapview)
# Convert your data to a spatial dataframe
dat1 <- data.frame(dat$REGNOM, dat$LON, dat$LAT)
dat2 <- na.omit(dat1)
names(dat2)[names(dat2) == "dat.REGNOM"] <- "Regions"  # Rename column to "Regions"
# Create spatial points for your data
Togo_points <- st_as_sf(dat2, coords = c("dat.LON", "dat.LAT"), crs = 4326)
# Download the boundary for Togo (level 0 for country border)
togo_boundary <- gadm(country = "TGO", level = 0, path = tempdir())
togo_boundary_sf <- st_as_sf(togo_boundary)  # Convert to sf object
# Plot the map with points and highlighted Togo borders
mapview(Togo_points, zcol = "Regions") + 
  mapview(togo_boundary_sf, color = "black", lwd = 2)  # Togo borders in black



#Figure 6 of Supplementary Material
#Where dat was used, E8 can be used instead
summary(dat$PREV_cr) #21 NA's
summary(dat$PREV_st) #198 NA's
plot(dat$PREV_cr,dat$PREV_st) #they seem to have a nice correlation
lm <- lm(dat$PREV_st~dat$PREV_cr) #the order matters for the plot line. 
summary(lm)
confint (lm)
plot(dat$PREV_cr,dat$PREV_st, ylab = "Crude Prevalence",xlab = "Standardised Prevalence", main = "Relationship between Crude Prevalence and Standardised Prevalence") #they seem to have a nice correlation
abline(lm) # to add the linear regression line to the graph
#standardized prevalence would be better (standardized for age and sex), but has more missing values than the crude prevalence.
#Formally see if they correlate well with the Pearson's Correlation Coefficient:
cor.test(dat$PREV_st,dat$PREV_cr, method="pearson", use="complete.obs")
#0.99 correlation coefficient is very acceptable. Crude prevalence should be the choice.
complete_obs_count <- sum(complete.cases(dat$PREV_st, dat$PREV_cr))
complete_obs_count #1612 datapoints


#Figure 7 of Supplementary Material
E888 <- E8
E888$Survey_year[E888$Survey_year >= 1975 & E888$Survey_year <= 1985] <- "1975-1985"
E888$Survey_year[E888$Survey_year >= 1986 & E888$Survey_year <= 1995] <- "1986-1995"
E888$Survey_year[E888$Survey_year >= 1996 & E888$Survey_year <= 2005] <- "1996-2005"
E888$Survey_year[E888$Survey_year >= 2006 & E888$Survey_year <= 2015] <- "2006-2015"
E888$Survey_year[E888$Survey_year >= 2016 & E888$Survey_year <= 2018] <- "2016-2018"

boxplot(d3*100~E888$Survey_year,
        main="Relationship between the propotion of the population surveyed and the survey year",
        xlab="year",
        ylab="Proportion of each village population surveyed by skin snip",
        col="orange",
        border="brown")

