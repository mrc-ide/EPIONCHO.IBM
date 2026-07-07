# Path to the unzipped "data_for_plots_only.zip", located in "data/"
data_folder <- "data/data_for_plots_only/"

#### Figure 1
make_fig_1 <- function(actual_data, projection_data) {
  mfp_plot <- ggplot() +
    geom_point(data=actual_data %>% filter(age_groups != 0.0), aes(x=age_groups, y=mf_prev*100, color="Observed Data")) +
    geom_errorbar(data=actual_data %>% filter(age_groups != 0.0), aes(x=age_groups, ymin=mf_wilson_lb*100, ymax=mf_wilson_ub*100), width=2.5, alpha=0.5) +
    geom_line(data=projection_data, aes(x=age_groups, y=mf_prev*100)) +
    geom_ribbon(data=projection_data, aes(x=age_groups, ymin=mf_q1*100, ymax=mf_q3*100), color="black", fill="black", linetype="dashed", alpha=0.2) +
    scale_color_manual("Type", values=c("black", "black")) +
    xlab("Age (years)") +
    ylab("Microfilarial prevalence (%)") +
    scale_x_continuous(breaks=seq(0, 80, 10)) +
    scale_y_continuous(breaks=seq(0, 20, 5), limits=c(0, 20.25)) +
    theme_bw() +
    theme(
      panel.grid.minor.x = element_blank()
    ) +
    guides(color="none") +
    facet_wrap(~exposure)
  print(mfp_plot)
}

fig_1_observed_data <- read.csv(paste0(data_folder, "figure_1_observed_atasame_data.csv"))
fig_1_model_data <- read.csv(paste0(data_folder, "figure_1_processed_model_data.csv"))

### Figure 1, with `use_onchosim = TRUE`
make_fig_1(
  fig_1_observed_data,
  fig_1_model_data
)

#### Figure 4

make_fig_4 <- function(observed_data, model_data) {
  ov16_plot <- ggplot() +
    geom_point(
      data=fig_4_observed_data,
      aes(x=age_groups, y=ov16_prev_rdt*100), alpha=0.5, color="black", size=0.5
    ) +
    geom_errorbar(
      data=fig_4_observed_data,
      aes(x=age_groups, ymin=(ov16_wilson_lb*100), ymax=(ov16_wilson_ub*100)), width=2.5, alpha=0.5, colour="black"
    ) +
    geom_line(
      data=fig_4_model_data,
      aes(x=age_groups, y=ov16_prev*100, color=Hypothesis, linetype=Hypothesis), linewidth=0.6
    ) +
    geom_ribbon(
      data=fig_4_model_data,
      aes(x=age_groups, ymin=ov16_lower*100, ymax=ov16_upper*100, fill=Hypothesis), color=NA, alpha=0.1
    ) +
    scale_color_manual(name="Hypotheses", values=c("black", "black")) +
    scale_fill_manual(name="Hypotheses", values=c("black", "black")) +
    scale_linetype_manual(name="Hypotheses", values=c("H4iii"="solid", "H4ii"="dotted")) +
    scale_alpha_continuous(range = c(0.3, 1)) +
    theme_bw() +
    xlab("Age (years)") +
    ylab("Anti-Ov16 seroprevalence (%)") +
    facet_wrap(~exposure) +
    guides(alpha="none", color="none", fill="none") +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      plot.title = element_text(size=10, hjust = 0.5),
      axis.title.y = element_text(size=10),
      title = element_text(size=10),
      legend.text = element_text(size=6),
      legend.title = element_text(size=8),
      legend.key.size = unit(0.8, 'lines'),
      axis.text.y = element_text(size=6),
      axis.text.x = element_text(size=6),
      strip.placement = "outside",
      panel.grid.minor = element_blank(),
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      strip.text.y = element_text(angle=0),
      axis.line.y.right = element_blank()
    )
  print(ov16_plot)
}

fig_4_observed_data <- read.csv(paste0(data_folder, "figure_4_observed_atasame_data.csv"))
fig_4_model_data <- read.csv(paste0(data_folder, "figure_4_processed_model_data.csv"))

make_fig_4(
  fig_4_observed_data,
  fig_4_model_data
)

#### Figure 5

make_fig_5 <- function(combined_df_2, descriptor="", hyps=c("H4i", "H4iii")) {
  ylab <- "all ages"
  pattern_vals <- c("none", "none", "stripe")
  names(pattern_vals) <- c("Observed", hyps[1], hyps[2])
  fill_vals <- c("grey", "white", "white")
  names(fill_vals) <- c("Observed", hyps[1], hyps[2])
  error_bar_colors <- c("Observed Error Bar"="black", "Error Bar"="black")
  combined_df_2 <- combined_df_2 %>%
    mutate(
      prefecture = factor(
        prefecture, levels=c("Oti/Kpendjal", "Keran", "Bassar")
      ),
      seroreversion = factor(
        seroreversion, levels=c("Observed", "H4ii", "H4iii")
      )
    )

  plot <- ggplot(data=combined_df_2) +
    geom_bar_pattern(aes(x=prefecture, y=ov16_prev*100, pattern=seroreversion, fill=seroreversion, group=seroreversion), pattern_fill="black", position=position_dodge(width=0.9), width=0.6, stat="identity", color="black", pattern_key_scale_factor=0.25) +
    geom_errorbar(aes(x=prefecture, ymin=ov16_prev_lb*100, ymax=ov16_prev_ub*100, group=seroreversion, color=error_bar_colors), width=0.3, position=position_dodge(width=0.9)) +
    scale_fill_manual("Model", values=fill_vals) +
    scale_pattern_manual("Model", values=pattern_vals) +
    scale_color_manual("", values=error_bar_colors) +
    xlab("Prefecture") +
    scale_y_continuous(
      paste("Anti-Ov16", ylab, "seroprevalence (%)"),
      limits = c(-5, max(combined_df_2$ov16_prev_ub*100)),
      breaks=seq(0, 60, 20)
    ) +
    guides(
      color="none",
      fill = guide_legend(title = ""),
      pattern = guide_legend(title = "")
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size=10, hjust = 0.5),
      axis.title.y = element_text(size=10),
      title = element_text(size=10),
      legend.text = element_text(size=6),
      legend.title = element_text(size=8),
      legend.key.size = unit(0.8, 'lines'),
      legend.position = "bottom",
      legend.direction = "horizontal",
      axis.text.y = element_text(size=6),
      axis.text.x = element_text(size=6, angle=70, vjust=0.5),
      strip.placement = "outside",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      strip.text.y = element_text(angle=0),
      axis.line.y.right = element_blank()
    ) +
    geom_segment(aes(x=0.55, xend=1.45, y=-5, yend=-5), color="black", size=0.1) +
    geom_segment(aes(x=0.55, xend=0.55, y=0, yend=-5), color="black", size=0.1) +
    geom_segment(aes(x=1.45, xend=1.45, y=0, yend=-5), color="black", size=0.1) +
    geom_segment(aes(x=1.55, xend=2.45, y=-5, yend=-5), color="black", size=0.1) +
    geom_segment(aes(x=1.55, xend=1.55, y=0, yend=-5), color="black", size=0.1) +
    geom_segment(aes(x=2.45, xend=2.45, y=0, yend=-5), color="black", size=0.1) +
    geom_segment(aes(x=2.55, xend=3.45, y=-5, yend=-5), color="black", size=0.1) +
    geom_segment(aes(x=2.55, xend=2.55, y=0, yend=-5), color="black", size=0.1) +
    geom_segment(aes(x=3.45, xend=3.45, y=0, yend=-5), color="black", size=0.1)
  print(plot)
}
figure_5_data <- read.csv(paste0(data_folder, "figure_5_data.csv"))

make_fig_5(
  figure_5_data,
  "both_all_ages_best_vc",
  hyps=c("H4ii", "H4iii")
)


#### Figure 7

make_fig_7 <- function(originalData, tmp_projection_data, descriptor, facet_value="sens ~ spec", facet_1, facet_2, best_fit_sens_spec) {
  color_vals <- c("H4i"="black", "H4ii"="black", "Midpoint"="black", "H4iii"="black")
  linetype_vals <- c("H4i"="solid", "H4ii"="dotted", "Midpoint"="dashed", "H4iii"="solid")
  p1 <- ggplot() +
    geom_line(
      aes(x=age_groups, y=ov16_prev*100,
          linetype=sero_type,
          color=sero_type
      ),
      size=0.5, data=tmp_projection_data
    ) +
    geom_ribbon(
      aes(x=age_groups, ymin=ov16_q1*100, ymax=ov16_q3*100,
          fill=sero_type,
          linetype=sero_type
      ),
      alpha=0.1,
      size=0.1,
      data=tmp_projection_data
    ) +
    geom_point(
      aes(x=age_groups, y=ov16_seroprev*100),
      color="black", alpha=0.5, data=originalData, size=0.5
    ) +
    geom_errorbar(
      aes(x=age_groups, ymin=ov16_wilson_lb*100, ymax=ov16_wilson_ub*100),
      width=2.5, alpha=0.5, color="black", data=originalData
    ) +
    xlab("Age (years)") + ylab("Anti-Ov16 seroprevalence (%)") +
    scale_color_manual("", values=color_vals) +
    scale_fill_manual("", values=color_vals) +
    scale_linetype_manual("", values=linetype_vals) +
    theme_bw() +
    scale_y_continuous(
      limits=c(0, 101),
      breaks=seq(0, 100, 25),
      sec.axis = sec_axis(~ . , name = facet_1, breaks = NULL, labels = NULL)
    ) +
    scale_x_continuous(
      limits=c(0, 81),
      breaks=seq(0, 100, 20)
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size=10, hjust = 0.5),
      axis.title.y = element_text(size=10),
      title = element_text(size=10),
      legend.text = element_text(size=6),
      legend.title = element_text(size=8),
      legend.key.size = unit(0.8, 'lines'),
      legend.position = "bottom",
      legend.direction = "horizontal",
      axis.text.y = element_text(size=6),
      axis.text.x = element_text(size=6),
      strip.placement = "outside",
      panel.grid.minor = element_blank(),
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      strip.text.y = element_text(angle=0),
      axis.line.y.right = element_blank()
    )
  print(p1)
}

komlan_data_fig_7 <- read.csv(paste0(data_folder, "figure_7_observed_data_komlan.csv"))
model_data_fig_7 <- read.csv(paste0(data_folder, "figure_7_processed_model_data.csv"))

make_fig_7(
  komlan_data_fig_7,
  model_data_fig_7,
  "best_vc_all", "~sens_spec",
  "",
  "Sensitivity, Specificity",
  best_fit_sens_spec = paste0(best_fit_sens, "%, ", round(best_fit_spec), "%")
)
