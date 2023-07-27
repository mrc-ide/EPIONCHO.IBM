library(dplyr)
library(ggplot2)
library(gridExtra)
library(fastR2)

getwd()
data <- read.csv("../Ov16 Data/MDA coverage per prefecture.csv")

unique(data$Prefecture)

colnames(data) <- c('Prefecture', 'Year', 'Endemicity', 'Total.Population', 'Population.at.Risk',
                    'Population.Treated', 'MDA.Delivered', 'MDA.Rounds', 'Cumulative.MDA.Years', 'Cumulative.MDA.Rounds',
                    'Effective.MDA.Coverage', 'Epidemilogical.Coverage', 'Data.Source', 'MDA.Geographical.Coverage')

countNonNAs <- function(x) {
  return(sum(is.na(x)))
}

countFunc <- function(x) {
  return(count(x))
}

nonNas <- data %>% mutate(
  Data.Source=case_when(
    Data.Source == "" ~ "No Source Provided",
    grepl('World Health Organization', Data.Source) ~ "WHO OCP",
    grepl('TOGO MoH', Data.Source) ~ "Togo MoH TIDC Database",
    TRUE ~ Data.Source)
) %>% group_by(Prefecture, Data.Source) %>% summarise_all(countNonNAs) %>% as.data.frame()

counts <- data %>% mutate(
  Data.Source=case_when(
    Data.Source == "" ~ "No Source Provided",
    grepl('World Health Organization', Data.Source) ~ "WHO OCP",
    grepl('TOGO MoH', Data.Source) ~ "Togo MoH TIDC Database",
    TRUE ~ Data.Source)
) %>% group_by(Prefecture, Data.Source) %>% summarise(count=n()) %>% as.data.frame() %>% arrange(desc(count))

# 1 = togoMoH, 2=all but togo, 3 = all
generatePlotsForPrefecture <- function(df, prefectures, colors, xval, yval, title="", ylabVal="", source=1, barGraph=TRUE, ybreaks=c(0,1)) {
  plots <- list()
  ylimits <- c(min(ybreaks), max(ybreaks))

  colorsToUse <- colors#sample(colors, length(prefectures))
  i <- 0
  for(prefecture in prefectures) {
    i <- i + 1
    tmpdf <- df %>% filter(Prefecture == prefecture & (grepl("TOGO MoH",Data.Source) & Data.Source != ""))
    if(source == 2) {
      tmpdf <- df %>% filter(Prefecture == prefecture & (!grepl("TOGO MoH",Data.Source) & Data.Source != ""))
    } else if(source == 3) {
      tmpdf <- df %>% filter(Prefecture == prefecture & (Data.Source != "")) %>% mutate(
        Epidemilogical.Coverage = ifelse((MDA.Delivered == 0) | is.na(MDA.Delivered), NA, ifelse(Epidemilogical.Coverage > 100, NA, Epidemilogical.Coverage)),
      ) %>% group_by(Prefecture, Year) %>% summarise(.groups="drop", Effective.MDA.Coverage = mean(Effective.MDA.Coverage, na.rm=TRUE), MDA.Delivered = mean(MDA.Delivered, na.rm=TRUE), Epidemilogical.Coverage=mean(Epidemilogical.Coverage, na.rm=TRUE), MDA.Rounds=mean(MDA.Rounds, na.rm=TRUE)) %>% as.data.frame()
    }
    tmpdf <- tmpdf %>% mutate(
      label=case_when(
        is.na(get(yval)) ~ 'NA',
        TRUE ~ ''
      )
    )
    plots[[i]] <- tmpdf %>% ggplot() +
      geom_line(aes(x=get(xval), y=get(yval)), color=colorsToUse[i]) +
      scale_y_continuous(name=ylabVal, breaks=ybreaks, limits=ylimits) +
      scale_x_continuous(name=xval, breaks=seq(1990, 2020, 5), limits=c(1988, 2020)) +
      theme(
        legend.position = "none",
        text=element_text(size=7),
        axis.text = element_text(size=7),
        axis.title= element_text(size=10)
      ) +
      ggtitle(paste(title, prefecture))
    if(barGraph) {
      plots[[i]] <- tmpdf %>% ggplot() +
        geom_bar(aes(x=get(xval), y=get(yval)), fill=colorsToUse[i], stat='identity') +
        scale_y_continuous(name=ylabVal, breaks=ybreaks, limits=ylimits) +
        scale_x_continuous(name=xval, breaks=seq(1990, 2020, 5), limits=c(1988, 2020)) +
        geom_text(aes(x=get(xval), y=0, label=label), position=position_dodge(width=0.9), size=2) +
        theme(
          legend.position = "none",
          text=element_text(size=7),
          axis.text = element_text(size=7),
          axis.title= element_text(size=10)
        ) +
        ggtitle(paste(title, prefecture))
    }
  }
  do.call("grid.arrange", c(plots))
  return(plots)
}

savePrefecturePlots <- function(prefecture, plots1, plots2, plots3, plots4, prefix="") {
  prefectures <- c('Bassar', 'Keran', 'KPENDJAL', 'OTI')
  prefecture_num <- which(prefecture == prefectures)
  arranged_plots <- grid.arrange(plots1[[prefecture_num]], plots2[[prefecture_num]], plots3[[prefecture_num]], plots4[[prefecture_num]])
  ggsave(paste("images/", prefix, "_", prefecture, '_plots.png', sep=""), arranged_plots, width=5.59, height=4.7)
}

colors <- c('darkgreen', 'darkblue', '#B05923', 'darkred')# grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]

# Togo MoH Source

mda_delivery_plots_raw <- generatePlotsForPrefecture(data, unique(data$Prefecture), colors, 'Year', 'MDA.Delivered', "Yearly MDA Delivery for", ylabVal="MDA Delivery")
mda_delivery_plots <- do.call("grid.arrange", c(mda_delivery_plots_raw))
#ggsave('images/mda_delivery_plots.png', mda_delivery_plots, width=1700, height=1200, units='px',dpi=300)

effective_mda_delivery_plots_raw <- generatePlotsForPrefecture(data, unique(data$Prefecture), colors, 'Year', 'Effective.MDA.Coverage', "Effective MDA Yearly Delivery for", ylabVal="Effective MDA Delivery")
effective_mda_delivery_plots <- do.call("grid.arrange", c(effective_mda_delivery_plots_raw))
#ggsave('images/effective_mda_delivery_plots.png', effective_mda_delivery_plots, width=1700, height=1200, units='px',dpi=300)

mda_rounds_plots_raw <- generatePlotsForPrefecture(data, unique(data$Prefecture), colors, 'Year', 'MDA.Rounds', "MDA Rounds for", ylabVal="Rounds", barGraph=FALSE, ybreaks=seq(0,3,1))
mda_rounds_plots <- do.call("grid.arrange", c(mda_rounds_plots_raw))
#ggsave('images/mda_rounds_plots.png', mda_rounds_plots, width=5.59, height=4.7)# width=3000, height=3000, units='px',dpi=400)

mda_coverage_plots_raw <- generatePlotsForPrefecture(data, unique(data$Prefecture), colors, 'Year', 'Epidemilogical.Coverage', "Epidemilogical MDA Coverage for", ylabVal="MDA Coverage (%)", barGraph=FALSE, ybreaks=seq(0,100, 10))
mda_coverage_plots <- do.call("grid.arrange", c(mda_coverage_plots_raw))
#ggsave('images/mda_coverage_plots.png', mda_coverage_plots, width=5.59, height=4.7)#width=3000, height=3000, units='px',dpi=400)

for(prefecture in unique(data$Prefecture)) {
  savePrefecturePlots(prefecture, mda_delivery_plots_raw,
                      effective_mda_delivery_plots_raw,
                      mda_rounds_plots_raw,
                      mda_coverage_plots_raw)
}

# Non Togo MoH Sources

non_moh_mda_delivery_plots_raw <- generatePlotsForPrefecture(data, unique(data$Prefecture), colors, 'Year', 'MDA.Delivered', "Yearly MDA Delivery for", ylabVal="MDA Delivery", source=2)
non_moh_mda_delivery_plots <- do.call("grid.arrange", c(non_moh_mda_delivery_plots_raw))
#ggsave('images/non_moh_mda_delivery_plots.png', non_moh_mda_delivery_plots, width=1700, height=1200, units='px',dpi=300)

non_moh_effective_mda_delivery_plots_raw <- generatePlotsForPrefecture(data, unique(data$Prefecture), colors, 'Year', 'Effective.MDA.Coverage', "Effective MDA Yearly Delivery for", ylabVal="Effective MDA Delivery", source=2)
non_moh_effective_mda_delivery_plots <- do.call("grid.arrange", c(non_moh_effective_mda_delivery_plots_raw))
#ggsave('images/non_moh_effective_mda_delivery_plots.png', non_moh_effective_mda_delivery_plots, width=1700, height=1200, units='px',dpi=300)

non_moh_mda_rounds_plots_raw <- generatePlotsForPrefecture(data, unique(data$Prefecture), colors, 'Year', 'MDA.Rounds', "MDA Rounds for", ylabVal="Rounds", barGraph=FALSE, ybreaks=seq(0,3,1), source=2)
non_moh_mda_rounds_plots <- do.call("grid.arrange", c(non_moh_mda_rounds_plots_raw))
#ggsave('images/non_moh_mda_rounds_plots.png', non_moh_mda_rounds_plots, width=3000, height=1500, units='px',dpi=400)

non_moh_mda_coverage_plots_raw <- generatePlotsForPrefecture(data, unique(data$Prefecture), colors, 'Year', 'Epidemilogical.Coverage', "Epidemilogical MDA Coverage for", ylabVal="MDA Coverage (%)", source=2, barGraph=FALSE, ybreaks=seq(0,100, 10))
non_moh_mda_coverage_plots <- do.call("grid.arrange", c(non_moh_mda_coverage_plots_raw))
#ggsave('images/non_moh_mda_coverage_plots.png', non_moh_mda_coverage_plots, width=3000, height=1500, units='px',dpi=400)

for(prefecture in unique(data$Prefecture)) {
  savePrefecturePlots(prefecture, non_moh_mda_delivery_plots_raw,
                      non_moh_effective_mda_delivery_plots_raw,
                      non_moh_mda_rounds_plots_raw,
                      non_moh_mda_coverage_plots_raw, prefix="non_moh_sources")
}

# All Sources

all_sources_mda_delivery_plots_raw <- generatePlotsForPrefecture(data, unique(data$Prefecture), colors, 'Year', 'MDA.Delivered', "Yearly MDA Delivery for", ylabVal="MDA Delivery", source=3)
all_sources_mda_delivery_plots <- do.call("grid.arrange", c(all_sources_mda_delivery_plots_raw))
#ggsave('images/all_sources_mda_delivery_plots.png', all_sources_mda_delivery_plots, width=1700, height=1200, units='px',dpi=300)

all_sources_effective_mda_delivery_plots_raw <- generatePlotsForPrefecture(data, unique(data$Prefecture), colors, 'Year', 'Effective.MDA.Coverage', "Effective MDA Yearly Delivery for", ylabVal="Effective MDA Delivery", source=3)
all_sources_effective_mda_delivery_plots <- do.call("grid.arrange", c(all_sources_effective_mda_delivery_plots_raw))
#ggsave('images/all_sources_effective_mda_delivery_plots.png', all_sources_effective_mda_delivery_plots, width=1700, height=1200, units='px',dpi=300)

all_sources_mda_rounds_plots_raw <- generatePlotsForPrefecture(data, unique(data$Prefecture), colors, 'Year', 'MDA.Rounds', "MDA Rounds for", ylabVal="Rounds", barGraph=FALSE, ybreaks=seq(0,3,1), source=3)
all_sources_mda_rounds_plots <- do.call("grid.arrange", c(all_sources_mda_rounds_plots_raw))
#ggsave('images/all_sources_mda_rounds_plots.png', all_sources_mda_rounds_plots, width=3000, height=1500, units='px',dpi=400)

all_sources_mda_coverage_plots_raw <- generatePlotsForPrefecture(data, unique(data$Prefecture), colors, 'Year', 'Epidemilogical.Coverage', "Epidemilogical MDA Coverage for", ylabVal="MDA Coverage (%)", source=3, barGraph=FALSE, ybreaks=seq(0,100, 10))
all_sources_mda_coverage_plots <- do.call("grid.arrange", c(all_sources_mda_coverage_plots_raw))
#ggsave('images/all_sources_mda_coverage_plots.png', all_sources_mda_coverage_plots, width=3000, height=1500, units='px',dpi=400)

for(prefecture in unique(data$Prefecture)) {
  savePrefecturePlots(prefecture, all_sources_mda_delivery_plots_raw,
                      all_sources_effective_mda_delivery_plots_raw,
                      all_sources_mda_rounds_plots_raw,
                      all_sources_mda_coverage_plots_raw, prefix="all_sources")
}


# Togo Serological Data

# all data from tables in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5849363/
togoMF <- function() {
  # 2 years each side except for 7.5 (2.5y gap) and 70.5 (9.5y gap)
  ageGroups <- c(7.5, 13, 18, 23, 28, 33, 38, 43, 48, 53, 58, 70.5)
  male_pop <- c(96, 106, 31, 32, 58, 45, 42, 45, 37, 32, 29, 42)
  female_pop <- c(76, 106, 64, 86, 137, 90, 63, 44, 53, 30, 32, 38)
  mf_pos_male <- c(1, 1, 0, 2, 9, 6, 6, 4, 3, 2, 1, 1)
  mf_pos_female <- c(3, 3, 1, 4, 7, 7, 8, 3, 3, 0, 3, 4)
  wilson_lb_m <- c()
  wilson_ub_m <- c()
  wilson_lb_f <- c()
  wilson_ub_f <- c()
  for(i in 1:length(ageGroups)) {
    wilson_ci_m <- wilson.ci(mf_pos_male[i], male_pop[i])
    wilson_ci_f <- wilson.ci(mf_pos_female[i], female_pop[i])

    wilson_lb_m <- c(wilson_lb_m, wilson_ci_m[1])
    wilson_ub_m <- c(wilson_ub_m, wilson_ci_m[2])
    wilson_lb_f <- c(wilson_lb_f, wilson_ci_f[1])
    wilson_ub_f <- c(wilson_ub_f, wilson_ci_f[2])
  }


  data <- data.frame(age_groups=ageGroups, male_pop=male_pop, female_pop=female_pop, mf_pos_male=mf_pos_male, mf_pos_female=mf_pos_female,
                     mf_wilson_lb_m=wilson_lb_m, mf_wilson_ub_m=wilson_ub_m,
                     mf_wilson_lb_f=wilson_lb_f, mf_wilson_ub_f=wilson_ub_f) %>%
    mutate(mf_prev_male = mf_pos_male/male_pop,
           mf_prev_female=mf_pos_female/female_pop,
           mf_wilson_lb_m=ifelse(mf_wilson_lb_m < 0, 0, mf_wilson_lb_m),
           mf_wilson_lb_f=ifelse(mf_wilson_lb_f < 0, 0, mf_wilson_lb_f))
  return(data)
}

togoOv16 <- function() {
  # 2 years each side except for 7.5 (2.5y gap) and 75.5 (4.5y gap)
  ageGroups <- c(7.5, 13, 18, 23, 28, 33, 38, 43, 48, 53, 58, 63, 68, 75.5)
  participants <- c(87, 76, 41, 44, 75, 58, 40, 27, 43, 24, 30, 13, 16, 2)
  ov16_pos <- c(13, 12, 7, 18, 31, 22, 20, 11, 15, 15, 17, 6, 9, 2)

  ageGroups <- c(0, 7.5, 13, 18, 23, 28, 33, 38, 43, 48, 53, 58, 63, 73)
  participants <- c(1, 87, 76, 41, 44, 75, 58, 40, 27, 43, 24, 30, 13, 18)
  ov16_pos <- c(0, 13, 12, 7, 18, 31, 22, 20, 11, 15, 15, 17, 6, 11)
  wilson_lb <- c()
  wilson_ub <- c()
  for(i in 1:length(ov16_pos)) {
    if(i==1) {
      wilson_lb <- c(wilson_lb, 0)
      wilson_ub <- c(wilson_ub, 0)
      next
    }
    wilson_ci <- wilson.ci(ov16_pos[i], participants[i])

    wilson_lb <- c(wilson_lb, wilson_ci[1])
    wilson_ub <- c(wilson_ub, wilson_ci[2])
  }
  data <- data.frame(age_groups=ageGroups, participants=participants, ov16_seropos=ov16_pos,
                     ov16_wilson_lb=wilson_lb, ov16_wilson_ub=wilson_ub) %>%
    mutate(ov16_seroprev = ov16_pos/participants)
  return(data)
}

togoVillage <- function() {
  riverBasin <- c(rep('Oti', 3), rep('Keran', 4), rep('Mo', 4))
  # https://journals.plos.org/plosntds/article/figure?id=10.1371/journal.pntd.0006312.g002
  # https://www.citypopulation.de/en/togo/admin/
  prefecture <- c('OTI', 'OTI', 'OTI', #c('KPENDJAL', 'KPENDJAL', 'OTI',#
                  'Keran', 'Keran', 'Keran', 'Keran',
                  'Bassar', 'Bassar', 'Bassar', 'Bassar')
  village <- c('Pancerys', 'Boutchakou', 'Koukoumbou',
               'Goulbi', 'Tchitchira', 'Koutougou Solla', 'Kpantiiyagou',
               'Bawlensi', 'Mo-Village', 'Katcha-Konkomba', 'Saboundi')
  mfExamined <- c(140, 127, 57,
                  131, 146, 81, 181,
                  148, 151, 188, 105)
  mfPositive <- c(4, 1, 3,
                  15, 15, 11, 14,
                  0, 13, 5, 2)
  ov16ElisaExaminedAllAges <- c(92, 92, 70,
                         92, 92, 79, 92,
                         NA, 22, NA, 13)
  ov16ElisaPositive <- c(10, 15, 9,
                         42, 42, 45, 46,
                         NA, 9, NA, 1)
  # Age 5-10
  ov16ElisaChildrenExamined <- c(25, 15, 17,
                         1, 1, 6, 22,
                         NA, NA, NA, NA)
  ov16ElisaChildrenPositive <- c(1, 1, 1,
                                 0, 0, 4, 6,
                                 NA, NA, NA, NA)
  df <- data.frame(river_basin=riverBasin, prefecture=prefecture, village=village, mf_pop=mfExamined,
                   mf_pos=mfPositive, ov16_pop_all_ages=ov16ElisaExaminedAllAges,
                   ov16_pos_all_ages=ov16ElisaPositive, ov16_pop_children=ov16ElisaChildrenExamined,
                   ov16_pos_children=ov16ElisaChildrenPositive)
}

togoMFData <- togoMF()

togo_mf_age_plot <- togoMFData %>% ggplot() + geom_point(aes(x=age_groups+0.5, y=mf_prev_male*100, color="Male")) +
  geom_point(aes(x=age_groups-0.5, y=mf_prev_female*100, color="Female")) +
  geom_errorbar(aes(x=age_groups-0.5, ymin=mf_wilson_lb_f*100, ymax=mf_wilson_ub_f*100, color="Female"), alpha=0.7) +
  geom_errorbar(aes(x=age_groups+0.5, ymin=mf_wilson_lb_m*100, ymax=mf_wilson_ub_m*100, color="Male"), alpha=0.7) +
  scale_y_continuous(breaks=seq(0, 30, 5), limits=c(0, 30)) +
  scale_x_continuous(breaks=seq(0, 80, 5), limits=c(0, 80)) +
  xlab("Age") +
  ylab("Microfilarial Prevalence (%)") +
  scale_color_manual(name="Sex", values=c("Male"='black', Female='red')) +
  ggtitle('Komlan et al. 2018 - Microfilarial Prevalence (Wilson 95% CI)')
togo_mf_age_plot

ggsave('images/togo_mf_age_plot.png', togo_mf_age_plot, width=5000, height = 4000, units="px", dpi=600)


togoOv16Data <- togoOv16()

togo_ov16_age_plot <- togoOv16Data %>% ggplot() + geom_point(aes(x=age_groups, y=ov16_seroprev*100)) +
  geom_errorbar(aes(x=age_groups, ymin=ov16_wilson_lb*100, ymax=ov16_wilson_ub*100)) +
  #geom_line(aes(x=age_groups, y=ov16_seroprev*100), color='red') +
  scale_y_continuous(breaks=seq(0, 100, 10), limits=c(0, 100)) +
  scale_x_continuous(breaks=seq(0, 80, 5), limits=c(0, 80)) +
  xlab("Age") +
  ylab("Ov16 Seroprevalence (%)") +
  ggtitle('Komlan et al. 2018 - Ov16 Seroprevalence (Wilson 95% CI)')
togo_ov16_age_plot

ggsave('images/togo_ov16_age_plot.png', togo_ov16_age_plot, width=5000, height = 4000, units="px", dpi=600)

togoVillageData <- togoVillage()

groupedVillageData <- togoVillageData %>% group_by(river_basin, prefecture) %>% summarise(across(-c('village'), ~ sum(.x, na.rm=TRUE)), .groups='drop') %>% as.data.frame() %>%
  mutate(
    mf_prev = mf_pos/mf_pop,
    ov16_prev_all_ages = ov16_pos_all_ages/ov16_pop_all_ages,
    ov16_prev_children = ov16_pos_children/ov16_pop_children
  )

togoVillageData <- togoVillageData %>% mutate(
  mf_prev=mf_pos/mf_pop
)

togo_village_data_mf_plot <- groupedVillageData %>% ggplot() +
  geom_bar(aes(x=river_basin, y=mf_prev*100, fill=river_basin), stat='identity') +
  scale_y_continuous(breaks=seq(0, 15, 3), limits=c(0, 12)) +
  scale_fill_manual(name="River Basin", values=c('darkgreen', 'darkblue', '#B05923')) +
  xlab("River Basin") +
  ylab("Microfilarial Prevalence (%)") +
  ggtitle('Komlan et al. 2018 - All Ages Microfilarial Prevalence by River Basin')
togo_village_data_mf_plot
ggsave('images/togo_village_data_mf_plot.png', togo_village_data_mf_plot, width=5000, height = 4000, units="px", dpi=600)


togo_village_data_ov16_all_ages_plot <- groupedVillageData %>% ggplot() +
  geom_bar(aes(x=river_basin, y=ov16_prev_all_ages*100, fill=river_basin), stat='identity') +
  scale_y_continuous(breaks=seq(0, 50, 5), limits=c(0, 50)) +
  scale_fill_manual(name="River Basin", values=c('darkgreen', 'darkblue', '#B05923')) +
  xlab("River Basin") +
  ylab("Ov16 Seroprevalence (%)") +
  ggtitle('Komlan et al. 2018 - All Ages Ov16 Seroprevalence by River Basin')
togo_village_data_ov16_all_ages_plot
ggsave('images/togo_village_data_ov16_all_ages_plot.png', togo_village_data_ov16_all_ages_plot, width=5000, height = 4000, units="px", dpi=600)



togo_village_data_ov16_children_plot <- groupedVillageData %>% ggplot() +
  geom_bar(aes(x=river_basin, y=ov16_prev_children*100, fill=river_basin), stat='identity') +
  scale_y_continuous(breaks=seq(0, 40, 5), limits=c(0, 40)) +
  scale_fill_manual(name="River Basin", values=c('darkgreen', '#B05923')) +
  xlab("River Basin") +
  ylab("Ov16 Seroprevalence (%)") +
  ggtitle('Komlan et al. 2018 - Children (5-10y) Ov16 Seroprevalence by River Basin')
togo_village_data_ov16_children_plot
ggsave('images/togo_village_data_ov16_children_plot.png', togo_village_data_ov16_children_plot, width=5000, height = 4000, units="px", dpi=600)

### Logistic regression

modified_data <- data %>% mutate(
  binned_coverage=case_when(
    is.na(Epidemilogical.Coverage) ~ "No Coverage",
    Epidemilogical.Coverage < 60 ~ "<60",
    Epidemilogical.Coverage < 80 ~ "60-80",
    Epidemilogical.Coverage < 95 ~ "80-95",
    TRUE ~ "95+"),
  rounded_coverage=case_when(
    is.na(Epidemilogical.Coverage) ~ 0,
    Epidemilogical.Coverage < 60 ~ 60,
    Epidemilogical.Coverage < 100 ~ ceiling(Epidemilogical.Coverage/5)*5,
    TRUE ~ 0
    ),
  Data.Source=case_when(
    Data.Source == "" ~ "No Source Provided",
    grepl('World Health Organization', Data.Source) ~ "WHO OCP",
    grepl('TOGO MoH', Data.Source) ~ "Togo MoH TIDC Database",
    TRUE ~ Data.Source)
  ) %>% filter(Data.Source == "Togo MoH TIDC Database") %>% select(Prefecture, Year, Endemicity, Data.Source, binned_coverage, rounded_coverage, Epidemilogical.Coverage)

prefecture_river_map <- togoVillageData %>% select(prefecture, river_basin) %>% unique()
#groupedVillageData

merged_data <- merge(modified_data, prefecture_river_map, by.x="Prefecture", by.y="prefecture")


fit <- lm(rounded_coverage ~ Year + Prefecture, data=merged_data)

summary(fit)


years_to_predict <- sort(unique(merged_data$Year))
fit_data <- data.frame(Years=rep(years_to_predict, 3), Prefecture=NA, cov=NA)
i <- 1
for(prefecture in c('OTI', 'Keran', 'Bassar')) {
  coverages <- fit$coefficients[[1]] + fit$coefficients[[2]] * years_to_predict
  if(prefecture == "OTI") {
    coverages <- coverages + fit$coefficients[[3]]
  } else if (prefecture == "Keran") {
    coverages <- coverages + fit$coefficients[[4]]
  }
  fit_data[i:(i+length(years_to_predict)-1), 'Prefecture'] <- prefecture
  fit_data[i:(i+length(years_to_predict)-1), 'cov'] <- coverages
  i <- i + length(years_to_predict)
}

fit_data %>% ggplot() +
  geom_line(aes(x=Years, y=cov, color=Prefecture)) +
  geom_point(aes(x=Year, y=Epidemilogical.Coverage, color=Prefecture), data=merged_data) +
  scale_y_continuous(breaks=seq(0, 100, 10), limits=c(0, 100)) +
  scale_x_continuous(breaks=seq(1990, 2020, 5), limits=c(1988, 2020)) +
  geom_vline(aes(xintercept=2003), linetype="dashed") +
  geom_vline(aes(xintercept=1996), linetype="dashed") +
  geom_vline(aes(xintercept=2000), color="darkgreen", linetype="dashed") +
  annotate("text", label="No MDA for all Prefectures", x=1996, y=10, hjust=-0.01) +
  annotate("text", label="Start of Twice a Year MDA", x=2003, y=95, hjust=-0.01) +
  annotate("text", label="Start of Consistent Once a Year MDA for Oti", x=2000, y=50, hjust=-0.01)


mda_scenarios <- merged_data %>% mutate(
  Epidemilogical.Coverage=ifelse(Epidemilogical.Coverage == 0, NA, Epidemilogical.Coverage)
  ) %>% mutate(
  binned_year=case_when(
    Year < 1995 ~ 1990,#mean(Epidemilogical.Coverage),
    Year < 2003 ~ 1996,#mean(Epidemilogical.Coverage),
    Year < 2018 ~ 2003,#mean(Epidemilogical.Coverage),
    TRUE ~ 2018
  )) %>% group_by(Prefecture, binned_year) %>% summarise(coverages=mean(Epidemilogical.Coverage, na.rm=TRUE)) %>%
  as.data.frame() %>%
  mutate(
    binned_year = case_when(
      Prefecture == "Bassar" ~ binned_year + 0.2,
      Prefecture == "Keran" ~ binned_year + 0.1,
      TRUE ~ binned_year
    ),
    coverages=case_when(
      is.nan(coverages) ~ 0,
      TRUE ~ coverages
    )
  )


mda_scenarios_plot <- mda_scenarios %>% ggplot() +
  geom_step(aes(x=binned_year, y=coverages, color=Prefecture)) +
  geom_point(aes(x=Year, y=Epidemilogical.Coverage, color=Prefecture), data=merged_data) +
  scale_y_continuous(breaks=seq(0, 100, 10), limits=c(0, 100)) +
  scale_x_continuous(breaks=seq(1990, 2020, 5), limits=c(1988, 2020)) +
  geom_vline(aes(xintercept=2003), linetype="dashed") +
  geom_vline(aes(xintercept=1996), linetype="dashed") +
  geom_vline(aes(xintercept=2000), color="darkgreen", linetype="dashed") +
  annotate("text", label="No MDA for all Prefectures in 1996", x=1996, y=10, hjust=-0.01) +
  annotate("text", label="Start of Twice a Year MDA", x=2003, y=95, hjust=-0.01) +
  annotate("text", label="Start of Consistent Once a Year MDA for Oti", x=2000, y=50, hjust=-0.01) +
  xlab("Year") +
  ylab("MDA Coverage (%)") +
  ggtitle("MDA Coverages in 3 Prefectures") +
  labs(clip = "off", hjust=0, caption="MDA Coverages calculated as a mean for each prefecture from 1988 - 1995, 1995-2002, and 2003+. NA values were removed from calculation.")
mda_scenarios_plot

ggsave("images/mda_scenarios_plot.png", mda_scenarios_plot, width=5000, height = 4000, units="px", dpi=500)
