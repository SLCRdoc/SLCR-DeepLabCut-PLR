#20260317 PLR analysis from DeepLabCut output
#This is an analysis script specific to data in KM Daly et al., Cell Reports 2026 
# and is meant to be run section by section with user input. 

library(readr)
library(dplyr)
library(tidyverse)
library(zoo)   # for rollmedian
library(ggpubr)


#Prompt for data paths
meta_path   <- readline(prompt = "Enter full path to metadata CSV file: ")
results_dir <- readline(prompt = "Enter path to directory with result CSVs: ") #Results are .csv output from DLC tracking pupil edges for all mice under a single light intensity
#NOTE: directory needs a \ at the end of path name

if (!file.exists(meta_path))  stop("Metadata file not found: ", meta_path)
if (!dir.exists(results_dir)) stop("Results directory not found: ", results_dir)

metadata <- read_csv(meta_path)

result_files <- list.files(
  path = results_dir,
  pattern = ".csv",
  full.names = FALSE
)

metadata
result_files

setwd(results_dir)


#Load in data. 

#Function: Extract video ID from DLC result names
extract_video_id <- function(path) {
  fname <- basename(path)
  stem  <- tools::file_path_sans_ext(fname)
  vid   <- str_extract(stem, "C\\d{4}|MAH\\d{5}")
  vid
}

#READ one DLC file, build proper column names
# hard coding column names
dlc_col_names <- c(
  "FrameNumber",
  "Nx", "Ny", "Nlikelihood",
  "Ex", "Ey", "Elikelihood",
  "Sx", "Sy", "Slikelihood",
  "Wx", "Wy", "Wlikelihood",
  "NEx", "NEy", "NElikelihood",
  "NWx", "NWy", "NWlikelihood",
  "SEx", "SEy", "SElikelihood",
  "SWx", "SWy", "SWlikelihood" )
  # extend if needed

#Function: to read one DLC file
read_dlc_file <- function(path) {
  raw <- readr::read_csv(path, col_names = FALSE, show_col_types = FALSE)
  # Drop the first 3 rows (DLC header rows)
  dat <- raw[-c(1, 2, 3), , drop = FALSE]
  # Keep only as many columns as we have names for
  n_cols <- min(ncol(dat), length(dlc_col_names))
  dat <- dat[, seq_len(n_cols), drop = FALSE]
  colnames(dat) <- dlc_col_names[seq_len(n_cols)]
  # Ensure numeric
  dat[] <- lapply(dat, as.numeric)
  
  dat
}

# Build nested dataframe: one row per Video, data = cleaned DLC tibble
nested_results <- tibble::tibble(
  file_path = result_files,
  Video     = extract_video_id(result_files)
) %>%
  dplyr::left_join(metadata, by = "Video") %>%
  dplyr::mutate(
    data = purrr::map(file_path, read_dlc_file)
  )

# Add a pixels per mm column to the nested tibble. NOTE Values are different between video batches!! calculated with FIJI and logged in PLR mastersheet xls
nested_results <- nested_results %>%
  mutate(
    pixels_per_mm = case_when(
      str_starts(Video, "C")  ~ 132, #Change value if necessary
      str_starts(Video, "MAH") ~ 92, #Change value if necessary
      TRUE                     ~ NA_real_
    )
  )

#Function: Add pair distances and apply a rolling median filter to data results in tibble
add_pair_distances <- function(df, window = 3) {
  df <- df %>%
    mutate(
      # logical masks: likelihoods > 0.5
      valid_NS    = Nlikelihood  > 0.5 & Slikelihood  > 0.5,
      valid_EW    = Elikelihood  > 0.5 & Wlikelihood  > 0.5,
      valid_NE_SW = NElikelihood > 0.5 & SWlikelihood > 0.5,
      valid_NW_SE = NWlikelihood > 0.5 & SElikelihood > 0.5,
      
      # raw distances with NA when likelihood too low
      dist_NS_raw    = if_else(valid_NS,
                               sqrt((Nx  - Sx)^2  + (Ny  - Sy)^2),
                               NA_real_),
      dist_EW_raw    = if_else(valid_EW,
                               sqrt((Ex  - Wx)^2  + (Ey  - Wy)^2),
                               NA_real_),
      dist_NE_SW_raw = if_else(valid_NE_SW,
                               sqrt((NEx - SWx)^2 + (NEy - SWy)^2),
                               NA_real_),
      dist_NW_SE_raw = if_else(valid_NW_SE,
                               sqrt((NWx - SEx)^2 + (NWy - SEy)^2),
                               NA_real_)
    )
  
  # helper that keeps NAs and only smooths runs of non‑NA values
  smooth_med <- function(x, k) {
    # rollmedian with na.pad = TRUE leaves NAs at the edges
    zoo::rollmedian(x, k = k, fill = NA, na.pad = TRUE)
  }
  
  df %>%
    mutate(
      dist_NS    = smooth_med(dist_NS_raw,    window),
      dist_EW    = smooth_med(dist_EW_raw,    window),
      dist_NE_SW = smooth_med(dist_NE_SW_raw, window),
      dist_NW_SE = smooth_med(dist_NW_SE_raw, window)
    ) %>%
    select(
      -valid_NS, -valid_EW, -valid_NE_SW, -valid_NW_SE,
      -dist_NS_raw, -dist_EW_raw, -dist_NE_SW_raw, -dist_NW_SE_raw
    )
}


nested_results <- nested_results %>%
  mutate(data = map(data, ~ add_pair_distances(.x, window = 3)))



# Flatten results for plotting (add mm conversion)
flat_results <- nested_results %>%
  mutate(
    light_onset_s = `Corrected light onset time (s)`,
    mID = mID,
    Sex = Sex,
    Treatment = Treatment,
    Virus = Virus,
    pixels_per_mm = pixels_per_mm
  ) %>%
  select(Video, mID, Sex, Treatment, Virus, light_onset_s, pixels_per_mm, data) %>%
  unnest(data) %>%
  mutate(
    time_rel = FrameNumber / 23.98 - light_onset_s, #manually change framerate if necessary. True rate for videos was 23.98 fps
    # Convert all distances to mm
    dist_NS_mm = dist_NS / pixels_per_mm,
    dist_EW_mm = dist_EW / pixels_per_mm,
    dist_NE_SW_mm = dist_NE_SW / pixels_per_mm,
    dist_NW_SE_mm = dist_NW_SE / pixels_per_mm
  ) %>%
  filter(time_rel >= -5, time_rel <= 30)


# VISUALIZE PUPIL DIAMETER COORDINATE PAIRS ON XY PLOT --------------------
# ### plot pairwise distances for all videos
# 
# make_video_plot <- function(one_vid, data = flat_results) {
#   plot_df <- data %>%
#     filter(Video == one_vid) %>%
#     select(time_rel, dist_NS_mm, dist_EW_mm, dist_NE_SW_mm) %>%
#     pivot_longer(
#       cols = c(dist_NS_mm, dist_EW_mm, dist_NE_SW_mm),
#       names_to = "pair",
#       values_to = "distance"
#     )
# 
#   # mean only if there are at least 2 non‑NA distances
#   avg_df <- plot_df %>%
#     group_by(time_rel) %>%
#     summarise(
#       n_pairs = sum(!is.na(distance)),
#       distance = ifelse(n_pairs >= 2, mean(distance, na.rm = TRUE), NA_real_),
#       .groups = "drop"
#     )
# 
#   # limit jumps in mean distance to max 0.25 mm per step
#   avg_df <- avg_df %>%
#     arrange(time_rel) %>%
#     mutate(
#       distance_clamped = {
#         d <- distance
#         # iterate through and clamp changes
#         if (length(d) > 1) {
#           for (i in 2:length(d)) {
#             if (!is.na(d[i]) && !is.na(d[i - 1])) {
#               delta <- d[i] - d[i - 1]
#               if (delta > 0.25)  d[i] <- d[i - 1] + 0.25
#               if (delta < -0.25) d[i] <- d[i - 1] - 0.25
#             }
#           }
#         }
#         d
#       }
#     )
# 
#   # baseline mean over -5 to 0 s (using pair-level data as before)
#   baseline_mean <- plot_df %>%
#     filter(time_rel >= -5, time_rel <= 0) %>%
#     summarise(mean_dist = mean(distance, na.rm = TRUE)) %>%
#     pull(mean_dist)
# 
#   ggplot() +
#     geom_line(
#       data = plot_df,
#       aes(x = time_rel, y = distance, color = pair),
#       alpha = 0.8
#     ) +
#     geom_line(
#       data = avg_df,
#       aes(x = time_rel, y = distance_clamped),
#       color = "black",
#       linewidth = 0.5
#     ) +
#     geom_hline(
#       yintercept = baseline_mean,
#       linetype = "dashed",
#       color = "gray30"
#     ) +
#     geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 0.8) +
#     theme_minimal() +
#     labs(
#       title = paste("Pair distances aligned to light onset for", one_vid),
#       x = "Time relative to light onset (s)",
#       y = "Pupil Diameter (mm)",
#       color = "Pair"
#     ) +
#     coord_cartesian(xlim = c(-5, 30), ylim = c(0, 3))
# }
# 
# 
# # Get all unique video IDs
# video_ids <- unique(flat_results$Video)
# 
# # Create a named list of plots, one per video
# plot_list <- map(video_ids, make_video_plot, data = flat_results)
# names(plot_list) <- video_ids
# 
# # (Optional) Save each plot to file
# walk(video_ids, ~ ggsave(
#   filename = paste0("DLCpairs2_", .x, ".jpeg"),
#   plot = plot_list[[.x]],
#   width = 6, height = 4, dpi = 300
# ))





# BASELINE PUPIL MEASUREMENTS, Remove Outliers ---------------------------------------------
# Compute baseline average (mean across all 4 pairs) for each video - NOW IN MM
baseline_summary <- flat_results %>%
  filter(time_rel >= -5, time_rel <= 0) %>%
  group_by(Video, Treatment, Sex, Virus) %>%
  summarise(
    baseline_dist = mean(c(dist_NS_mm, dist_EW_mm, dist_NE_SW_mm), na.rm = TRUE),
    n_frames = n(),
    .groups = "drop"
  )

# Plot with boxplots by Treatment, points colored by Virus and shaped by Sex
ggplot(baseline_summary, aes(x = Treatment, y = baseline_dist, fill = Treatment)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.fill = "white") +
  geom_point(
    aes(color = Virus, shape = Sex),
    position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
    alpha = 0.7, size = 2.5
  ) +
  scale_fill_manual(values = c("CNO" = "grey80", "Saline" = "white")) +
  scale_color_manual(values = c("q" = "#40bfff", "i" = "#ff7e79")) +
  scale_shape_manual(values = c("m" = 15, "f" = 19)) +
  theme_minimal() +
  labs(
    title = "Baseline Pupil Diameter (mean of 4 pairs, -5 to 0s relative to light onset)",
    x = "Treatment",
    y = "Average Diameter (mm)",  # Changed units
    fill = "Treatment",
    color = "Virus",
    shape = "Sex"
  ) +
  theme(legend.position = "bottom")

# list the pupil baselines per video
baseline_summary %>%
  select(Video, Treatment, Sex, Virus, baseline_dist, n_frames) %>%
  mutate(baseline_dist_mm = round(baseline_dist, 3)) %>%
  arrange(Video) %>%
  print(n = Inf)

# Calculate outliers and list video ID and baseline pupil size
baseline_summary %>%
  group_by(Treatment) %>%
  filter(
    !is.na(baseline_dist),  # drop NAs before comparing
    baseline_dist < quantile(baseline_dist, 0.25, na.rm = TRUE) - 1.5 * IQR(baseline_dist, na.rm = TRUE) |
      baseline_dist > quantile(baseline_dist, 0.75, na.rm = TRUE) + 1.5 * IQR(baseline_dist, na.rm = TRUE)
  ) %>%
  select(Video, baseline_dist) %>%
  print(n = Inf)


# Manually define outliers with preconstricted pupils or NANs to exclude from above baseline and individual TS
#outlier_videos <- c("C0082", "MAH00095") #1 lux
#outlier_videos <- c("C0015", "C0018") #10lux
#outlier_videos <- c("C0085", "MAH00075", "MAH00103") #100lux "c0082", "c0071",
#outlier_videos <- c("C0047", "C0050") #1000lux


# Filter nested_results to exclude outliers
nested_results_clean <- nested_results %>%
  filter(!Video %in% outlier_videos)

# Now recreate flat_results without outliers
flat_results <- nested_results_clean %>%
  mutate(
    light_onset_s = `Corrected light onset time (s)`,
    mID = mID,
    Sex = Sex,
    Treatment = Treatment,
    Virus = Virus,
    pixels_per_mm = pixels_per_mm
  ) %>%
  select(Video, Sex, mID, Treatment, Virus, light_onset_s, pixels_per_mm, data) %>%
  unnest(data) %>%
  mutate(
    time_rel = FrameNumber / 23.98 - light_onset_s,
    dist_NS_mm = dist_NS / pixels_per_mm,
    dist_EW_mm = dist_EW / pixels_per_mm,
    dist_NE_SW_mm = dist_NE_SW / pixels_per_mm,
    dist_NW_SE_mm = dist_NW_SE / pixels_per_mm
  ) %>%
  filter(time_rel >= -5, time_rel <= 30)







# TIMESERIES PUPIL DIAMETER PLOTS per mouse or avg virus ------------------
## Avg for each condition e.g. X-Gi;Saline, X-Gi;CNO, X-Gq;Saline, X-Gq;CNO for a single mouse 
# 
# mouseID <- "6231" #manually update this as necessary
# 
# avg_timeseries_virus <- flat_results %>%
#   filter(mID == mouseID) %>% # comparing within one mouse. comment out other filters!
#   #filter(Virus == 'q') %>% # Comparing AVG pupil diameter of all mice expressing same virus (q or i). #Comment out other filters! 
#   #NOTE if both filters above commented out, will plot all mice in folder
#   
#   filter(time_rel >= -5, time_rel <= 30) %>%
#   select(Video, Treatment, mID, Virus, time_rel, dist_NS_mm, dist_EW_mm, dist_NE_SW_mm, dist_NW_SE_mm) %>%
#   pivot_longer(
#     cols = c(dist_NS_mm, dist_EW_mm, dist_NE_SW_mm, dist_NW_SE_mm),
#     names_to = "pair", values_to = "distance"
#   ) %>%
#   mutate(time_bin = round(time_rel * 2) / 2) %>%  # smooth 0.5s bins
#   group_by(Virus, Treatment, time_bin) %>%
#   summarise(
#     mean_distance = mean(distance, na.rm = TRUE),
#     sd_distance = sd(distance, na.rm = TRUE),
#     n_videos = n_distinct(Video),
#     .groups = "drop"
#   ) %>%
#   rename(time_rel = time_bin) %>%
#   unite(group, Virus, Treatment, sep = "_")  # combine for cleaner legend
# # Plot ALL 4 lines on single graph (Saline trial are separate here)
# ggplot(avg_timeseries_virus, aes(x = time_rel, y = mean_distance, color = group)) +
#   geom_line(linewidth = 1.5, alpha = 1) +
#   geom_ribbon(aes(ymin = pmax(mean_distance - sd_distance, 0),
#                   ymax = mean_distance + sd_distance, fill = group),
#               alpha = 0.25, color = NA) +
#   geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 1) +
# 
#   # Custom colors: q_Saline=gray, q_CNO=blue, i_Saline=red, i_CNO=orange
#   scale_color_manual(
#     values = c("q_Saline" = "gray50",
#                "q_CNO" = "blue",
#                "i_Saline" = "black",
#                "i_CNO" = "red")
#   ) +
#   scale_fill_manual(
#     values = c("q_Saline" = alpha("gray50", 0.25),
#                "q_CNO" = alpha("blue", 0.25),
#                "i_Saline" = alpha("black", 0.25),
#                "i_CNO" = alpha("red", 0.25))
#   ) +
# 
#   theme_minimal(base_size = 14) +
#   theme(
#     legend.position = "bottom",
#     panel.grid.minor = element_blank(),
#     legend.title = element_blank()
#   ) +
#   labs(
#     title = "Pupil Diameter by Virus and Treatment",
#     subtitle = "Shaded = ±1 SD | Dashed = light onset (t=0)",
#     x = "Time relative to light onset (s)",
#     y = "Pupil diameter (mm)"
#   ) +
#   coord_cartesian(xlim = c(-5, 30), ylim = c(0, NA))




# FINAL FIGURE TIMESERIES PLOTS -------------------------------------------
# Average time series for all (Saline from both viruses combined into one line)
avg_timeseries_virus <- flat_results %>%
  filter(time_rel >= -5, time_rel <= 30) %>%
  select(Video, Treatment, Virus, time_rel, dist_NS_mm, dist_EW_mm, dist_NE_SW_mm) %>%
  pivot_longer(
    cols = c(dist_NS_mm, dist_EW_mm, dist_NE_SW_mm),
    names_to = "pair", values_to = "distance"
  ) %>%
  # Combine ALL Saline videos (q+i) into one group
  mutate(Treatment_simple = if_else(Treatment == "Saline", "Saline", paste0(Virus, "_CNO"))) %>%
  mutate(time_bin = round(time_rel * 10) / 10) %>% #10 chosen for 0.1s bin interval. Cleaner plot without losing sensitivity
  group_by(Treatment_simple, time_bin) %>%
  summarise(
    n_pairs = sum(!is.na(distance)),
    mean_distance = ifelse(n_pairs >= 2, mean(distance, na.rm = TRUE), NA_real_),
    sd_distance = ifelse(n_pairs >= 2, sd(distance, na.rm = TRUE), NA_real_),
    n_videos = n_distinct(Video),
    .groups = "drop"
  ) %>%
  rename(time_rel = time_bin)


# Plot 3 clean lines on one graph
p <- ggplot(avg_timeseries_virus, aes(x = time_rel, y = mean_distance, color = Treatment_simple)) +
  geom_line(linewidth = 1, alpha = 1) +
  geom_ribbon(aes(ymin = pmax(mean_distance - sd_distance, 0), 
                  ymax = mean_distance + sd_distance, fill = Treatment_simple), 
              alpha = 0.2, color = NA) +
  
  # Saline upper / lower bounds as dotted lines
  geom_line(
    data = avg_timeseries_virus %>% filter(Treatment_simple == "Saline"),
    aes(y = pmax(mean_distance - sd_distance, 0), colour = 'black'),
    linetype = "dotted",
    linewidth = 0.6,
    show.legend = FALSE
  ) +
  geom_line(
    data = avg_timeseries_virus %>% filter(Treatment_simple == "Saline"),
    aes(y = mean_distance + sd_distance),
    linetype = "dotted",
    linewidth = 0.6,
    show.legend = FALSE
  ) +
   #light onset line
  geom_vline(xintercept = 0, color = "black", linetype = "solid", linewidth = 1) +
  
  # Colors: gray Saline (all viruses), blue q_CNO, red i_CNO
  scale_color_manual(
    values = c("Saline" = "gray50", 
               "q_CNO" = "#40bfff", 
               "i_CNO" = "#ff7e79")
  ) +
  scale_fill_manual(
    values = c("Saline" = alpha("white", 0), 
               "q_CNO" = alpha("#40bfff", 0.02), 
               "i_CNO" = alpha("#ff7e79", 0.02))
  ) +
  
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  # labs(
  #   title = "1 Lux Pupil Diameter Time Course",
  #   subtitle = "Gray = Saline (Gq + Gi combined) | Shaded Error = ±1 SD | Dashed = light onset",
  #   x = "Time relative to light onset (s)",
  #   y = "Pupil diameter (mm)"
  # ) +
  coord_cartesian(xlim = c(-5, 30), ylim = c(0, 2.5))

p





