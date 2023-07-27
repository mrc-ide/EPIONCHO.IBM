library(dplyr)
library(ggplot2)

readGabonDataRDTElisa <- function() {
  fileName <- "../Ov16 Data/Ov16 Gabon ELISA.csv"
  data <- read.csv(fileName) %>% select(1:14) %>% drop_na()
  colnames(data) <- c('id', 'ov16_status_rdt', 'mf_status', 'sample_id', 'inferred_conc', 'region', 'district_name', 'village_name', 'ov16_elisa_plate_num', 'sample_id_ind_id_mismatch', 'ov16_result_elisa', 'ov16_status_elisa', 'age', 'age_group_orig')
  # view duplicate individual id vals, remove them as we don't know which one is correct
  dupIds <- as.numeric(names(which(table(data$id) > 1)))
  data[data$id %in% dupIds,]
  data <- data %>% filter(!(id %in% dupIds)) %>% mutate(
    age_groups = case_when(
      age <= 30 ~ round(age/5)*5,
      age <= 75 ~ round(age/5)*5,
      TRUE ~ 80
    ),
    ov16_status_elisa = case_when(
      ov16_result_elisa == 'Invalid' ~ as.numeric(NaN),
      TRUE ~ as.numeric(ov16_status_elisa)
    )
  ) %>% drop_na()
  print(paste("Total Individuals:", length(unique(data$id))))
  print(paste("Overall Gabon MF Prevalence:", mean(data$mf_status)))
  print(paste("Overall Gabon Ov16 ELISA SeroPrevalence:", mean(data$ov16_status_elisa, na.rm=TRUE)))
  print(paste("Overall Gabon Ov16 RDT SeroPrevalence:", mean(data$ov16_status_rdt)))

  return(data)
}

data <- readGabonDataRDTElisa()

ggplot() + geom_histogram(aes(x=inferred_conc, fill=factor(ov16_status_elisa)), data=data[data$mf_status == 0 & data$ov16_status_elisa == 0,]) +
  geom_histogram(aes(x=inferred_conc, fill=factor(ov16_status_elisa)), data=data[data$mf_status == 0 & data$ov16_status_elisa == 1,], alpha=0.5) +
  ylim(0, 100) +
  scale_x_continuous(breaks=seq(0, 150, 10), limits=c(0, 150)) +
  ggtitle("MF neg")

ggplot() + geom_histogram(aes(x=inferred_conc, fill=factor(ov16_status_elisa)), data=data[data$mf_status == 1 & data$ov16_status_elisa == 0,]) +
  geom_histogram(aes(x=inferred_conc, fill=factor(ov16_status_elisa)), data=data[data$mf_status == 1 & data$ov16_status_elisa == 1,], alpha=0.5) +
  ylim(0, 100) +
  scale_x_continuous(breaks=seq(0, 150, 10), limits=c(0, 150)) +
  ggtitle("Mf pos")
