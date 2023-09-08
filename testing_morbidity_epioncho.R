# ========================================================================================================================== #
#                                 TESTING MORBIDITY (SKIN DISEASE) - 22/05/23                                                #
# ========================================================================================================================== #

devtools::load_all()

#length of simulation in years
timesteps = 200
timesteps = 130

#should treatment be given, when and how often
give.treat.in = 0
treat.strt = 1; treat.stp = 16
trt.int = 1

#annual biting rate, which determines infection prevalence (60% microfilarae prevalence)
ABR.in = 1082 # 60% mf prev
ABR.in = 2297 # 70% mf prev
ABR.in = 6969 # 80% mf prev
ABR.in = 607 # 50% mf prev

output_equilibrium_morbidity <- ep.equi.sim(time.its = timesteps,
                                  ABR = ABR.in,
                                  N.in = 440,
                                  treat.int = trt.int,
                                  treat.prob = 0.65,
                                  give.treat = give.treat.in,
                                  treat.start = treat.strt,
                                  treat.stop = treat.stp,
                                  pnc = 0.05,
                                  min.mont.age = 5,
                                  vector.control.strt = NA, #
                                  delta.hz.in = 0.186,
                                  delta.hinf.in = 0.003,
                                  c.h.in = 0.005,
                                  gam.dis.in = 0.3,
                                  run_equilibrium = TRUE,
                                  equilibrium,
                                  print_progress = TRUE,
                                  morbidity_module2 = "YES")

names(output_equilibrium_morbidity)

#output_equilibrium_morbidity_mf60prev <- output_equilibrium_morbidity
#output_equilibrium_morbidity_mf70prev <- output_equilibrium_morbidity
#output_equilibrium_morbidity_mf80prev <- output_equilibrium_morbidity
#output_equilibrium_morbidity_mf50prev <- output_equilibrium_morbidity
#output_equilibrium_morbidity_chesnaisABR <- output_equilibrium_morbidity

# saveRDS(output_equilibrium_morbidity_mf50prev , "output_equilibrium_morbidity_mf50prev.rds")
# saveRDS(output_equilibrium_morbidity_mf60prev , "output_equilibrium_morbidity_mf60prev.rds")
# saveRDS(output_equilibrium_morbidity_mf70prev , "output_equilibrium_morbidity_mf70prev.rds")
# saveRDS(output_equilibrium_morbidity_mf80prev , "output_equilibrium_morbidity_mf80prev.rds")

# output_equilibrium_morbidity <- output_equilibrium_morbidity_mf70prev

tme <- seq(1, 200*366-1)/366
tme <- seq(1, 130*366-1)/366
tme <- seq(1, 80*366-1)/366

plot(tme, output_equilibrium_morbidity$mf_prev, type = 'l', xlab = 'time (years)', ylab = 'microfilarial prevalence', ylim = c(0, 1))

plot(tme, output_equilibrium_morbidity$mf_intens, type = 'l', xlab = 'time (years)', ylab = 'mean microfilarial intensity')

plot(tme, output_equilibrium_morbidity$L3, type = 'l', xlab = 'time (years)', ylab = 'mean L3 in fly population')

# ============ #
#tme2 <- seq(1, 130*366)/366
#tme2 <- seq(1, 200*366)/366

#tme2 <- seq(1, 130*12)/366
tme2 <- seq(1, length(output_equilibrium_morbidity$severe_itch_prev))/12

plot(tme2, output_equilibrium_morbidity$severe_itch_prev, type = 'l', xlab = 'time (years)', ylab = 'Severe itch prevalence', ylim = c(0, 1))
tail(output_equilibrium_morbidity$severe_itch_prev)
mean(tail(output_equilibrium_morbidity$severe_itch_prev))
plot(tme2, output_equilibrium_morbidity$RSD_prev, type = 'l', xlab = 'time (years)', ylab = 'Reactive skin disease prevalence', ylim = c(0, 1))
tail(output_equilibrium_morbidity$RSD_prev)
mean(tail(output_equilibrium_morbidity$RSD_prev))
plot(tme2, output_equilibrium_morbidity$atrophy_prev, type = 'l', xlab = 'time (years)', ylab = 'Skin atrophy prevalence', ylim = c(0, 1))
tail(output_equilibrium_morbidity$atrophy_prev)
mean(tail(output_equilibrium_morbidity$atrophy_prev))
plot(tme2, output_equilibrium_morbidity$hanging_groin_prev, type = 'l', xlab = 'time (years)', ylab = 'Hanging groin prevalence', ylim = c(0, 1))
tail(output_equilibrium_morbidity$hanging_groin_prev)
mean(tail(output_equilibrium_morbidity$hanging_groin_prev))
plot(tme2, output_equilibrium_morbidity$mild_depigmentation_prev, type = 'l', xlab = 'time (years)', ylab = 'mild depigmentation prevalence', ylim = c(0, 1))
tail(output_equilibrium_morbidity$mild_depigmentation_prev)
mean(tail(output_equilibrium_morbidity$mild_depigmentation_prev))
plot(tme2, output_equilibrium_morbidity$severe_depigmentation_prev, type = 'l', xlab = 'time (years)', ylab = 'severe depigmentation prevalence', ylim = c(0, 1))
tail(output_equilibrium_morbidity$severe_depigmentation_prev)
mean(tail(output_equilibrium_morbidity$severe_depigmentation_prev))


plot(tme, output_equilibrium_morbidity$atrophy_prev, type = 'l', xlab = 'time (years)', ylab = 'Skin atrophy prevalence', ylim = c(0, 0.1), xlim = c(70,130))
plot(tme, output_equilibrium_morbidity$hanging_groin_prev, type = 'l', xlab = 'time (years)', ylab = 'Hanging groin prevalence', ylim = c(0, 0.1), xlim = c(70,130))
plot(tme, output_equilibrium_morbidity$mild_depigmentation_prev, type = 'l', xlab = 'time (years)', ylab = 'mild depigmentation', ylim = c(0, 0.1), xlim = c(70,130))

output_equilibrium_morbidity_60mfprev <- output_equilibrium_morbidity
output_equilibrium_morbidity_70mfprev <- output_equilibrium_morbidity


output_equilibrium_morbidity_60mfprev_v2 <- output_equilibrium_morbidity
output_equilibrium_morbidity_70mfprev_v2 <- output_equilibrium_morbidity

#output_equilibrium_morbidity <- output_equilibrium_morbidity_60mfprev_v2

# =============================================================================================================================== #
#                              Analysing / testing age-prevalence of skin disease                                                 #

#output_equilibrium_morbidity <- output_equilibrium_morbidity_70mfprev
# =============================================================================================================================== #
#                              Analysing / testing age-prevalence of skin disease                                                 #

#output_equilibrium_morbidity <- output_equilibrium_morbidity_70mfprev

names(output_equilibrium_morbidity)

ageprev_out <- output_equilibrium_morbidity$morbidity_ageprev_out

# SI modelled age-prev #
SI_0_1 <- mean(tail(ageprev_out$SI_prev0_1))
SI_2_4 <- mean(tail(ageprev_out$SI_prev2_4))
SI_5_9 <- mean(tail(ageprev_out$SI_prev5_9))
SI_10_19 <- mean(tail(ageprev_out$SI_prev10_19))
SI_20_29 <- mean(tail(ageprev_out$SI_prev20_29))
SI_30_49 <- mean(tail(ageprev_out$SI_prev30_49))
SI_50_80 <- mean(tail(ageprev_out$SI_prev50_80))

SI_ageprev_vec <- c(SI_0_1, SI_2_4, SI_5_9, SI_10_19, SI_20_29, SI_30_49, SI_50_80)
age_vec <- c(1, 3, 7.5, 15, 25, 40, 65)
SI_ageprev_df <- data.frame(age = age_vec, prop = SI_ageprev_vec)
SI_ageprev_df$prev <- SI_ageprev_df$prop * 100
SI_ageprev_df$sequela <- as.factor("severe itching")
SI_ageprev_df$agecat <- c("0 - 1", "2 - 4", "5 - 9", "10 - 19", "20 - 29", "30 - 49", "50 - 80")
SI_ageprev_df$agecat2 <- ordered(SI_ageprev_df$agecat, levels = c("0 - 1", "2 - 4", "5 - 9", "10 - 19", "20 - 29", "30 - 49", "50 - 80"))

require(ggplot2)
ggplot(data = SI_ageprev_df, aes(x = agecat2, y =prev, group = sequela))+
  geom_point()+
  geom_line()

# RSD modelled age-prev #
RSD_0_1 <- mean(tail(ageprev_out$RSD_prev0_1))
RSD_2_4 <- mean(tail(ageprev_out$RSD_prev2_4))
RSD_5_9 <- mean(tail(ageprev_out$RSD_prev5_9))
RSD_10_19 <- mean(tail(ageprev_out$RSD_prev10_19))
RSD_20_29 <- mean(tail(ageprev_out$RSD_prev20_29))
RSD_30_49 <- mean(tail(ageprev_out$RSD_prev30_49))
RSD_50_80 <- mean(tail(ageprev_out$RSD_prev50_80))

RSD_ageprev_vec <- c(RSD_0_1, RSD_2_4, RSD_5_9, RSD_10_19, RSD_20_29, RSD_30_49, RSD_50_80)
age_vec <- c(1, 3, 7.5, 15, 25, 40, 65)
RSD_ageprev_df <- data.frame(age = age_vec, prop = RSD_ageprev_vec)
RSD_ageprev_df$prev <- RSD_ageprev_df$prop * 100
RSD_ageprev_df$sequela <- as.factor("RSD")
RSD_ageprev_df$agecat <- c("0 - 1", "2 - 4", "5 - 9", "10 - 19", "20 - 29", "30 - 49", "50 - 80")
RSD_ageprev_df$agecat2 <- ordered(RSD_ageprev_df$agecat, levels = c("0 - 1", "2 - 4", "5 - 9", "10 - 19", "20 - 29", "30 - 49", "50 - 80"))

ggplot(data = RSD_ageprev_df, aes(x = agecat2, y =prev, group = sequela))+
  geom_point()+
  geom_line()

# Atrophy modelled age-prev #
Atrp_0_1 <- mean(tail(ageprev_out$Atrp_prev0_1))
Atrp_2_4 <- mean(tail(ageprev_out$Atrp_prev2_4))
Atrp_5_9 <- mean(tail(ageprev_out$Atrp_prev5_9))
Atrp_10_19 <- mean(tail(ageprev_out$Atrp_prev10_19))
Atrp_20_29 <- mean(tail(ageprev_out$Atrp_prev20_29))
Atrp_30_49 <- mean(tail(ageprev_out$Atrp_prev30_49))
Atrp_50_80 <- mean(tail(ageprev_out$Atrp_prev50_80))

Atrp_ageprev_vec <- c(Atrp_0_1, Atrp_2_4, Atrp_5_9, Atrp_10_19, Atrp_20_29, Atrp_30_49, Atrp_50_80)
age_vec <- c(1, 3, 7.5, 15, 25, 40, 65)
Atrp_ageprev_df <- data.frame(age = age_vec, prop = Atrp_ageprev_vec)
Atrp_ageprev_df$prev <- Atrp_ageprev_df$prop * 100
Atrp_ageprev_df$sequela <- as.factor("atrophy")
Atrp_ageprev_df$agecat <- c("0 - 1", "2 - 4", "5 - 9", "10 - 19", "20 - 29", "30 - 49", "50 - 80")
Atrp_ageprev_df$agecat2 <- ordered(Atrp_ageprev_df$agecat, levels = c("0 - 1", "2 - 4", "5 - 9", "10 - 19", "20 - 29", "30 - 49", "50 - 80"))

ggplot(data = Atrp_ageprev_df, aes(x = agecat2, y =prev, group = sequela))+
  geom_point()+
  geom_line()

# hanging groin modelled age-prev #
HG_0_1 <- mean(tail(ageprev_out$HG_prev0_1))
HG_2_4 <- mean(tail(ageprev_out$HG_prev2_4))
HG_5_9 <- mean(tail(ageprev_out$HG_prev5_9))
HG_10_19 <- mean(tail(ageprev_out$HG_prev10_19))
HG_20_29 <- mean(tail(ageprev_out$HG_prev20_29))
HG_30_49 <- mean(tail(ageprev_out$HG_prev30_49))
HG_50_80 <- mean(tail(ageprev_out$HG_prev50_80))

HG_ageprev_vec <- c(HG_0_1, HG_2_4, HG_5_9, HG_10_19, HG_20_29, HG_30_49, HG_50_80)
age_vec <- c(1, 3, 7.5, 15, 25, 40, 65)
HG_ageprev_df <- data.frame(age = age_vec, prop = HG_ageprev_vec)
HG_ageprev_df$prev <- HG_ageprev_df$prop * 100
HG_ageprev_df$sequela <- as.factor("hanging groin")
HG_ageprev_df$agecat <- c("0 - 1", "2 - 4", "5 - 9", "10 - 19", "20 - 29", "30 - 49", "50 - 80")
HG_ageprev_df$agecat2 <- ordered(HG_ageprev_df$agecat, levels = c("0 - 1", "2 - 4", "5 - 9", "10 - 19", "20 - 29", "30 - 49", "50 - 80"))

ggplot(data = HG_ageprev_df, aes(x = agecat2, y =prev, group = sequela))+
  geom_point()+
  geom_line()

# depigmentation modelled age-prev #
depigm_0_1 <- mean(tail(ageprev_out$depigm_prev0_1))
depigm_2_4 <- mean(tail(ageprev_out$depigm_prev2_4))
depigm_5_9 <- mean(tail(ageprev_out$depigm_prev5_9))
depigm_10_19 <- mean(tail(ageprev_out$depigm_prev10_19))
depigm_20_29 <- mean(tail(ageprev_out$depigm_prev20_29))
depigm_30_49 <- mean(tail(ageprev_out$depigm_prev30_49))
depigm_50_80 <- mean(tail(ageprev_out$depigm_prev50_80))

depigm_ageprev_vec <- c(depigm_0_1, depigm_2_4, depigm_5_9, depigm_10_19, depigm_20_29, depigm_30_49, depigm_50_80)
age_vec <- c(1, 3, 7.5, 15, 25, 40, 65)
depigm_ageprev_df <- data.frame(age = age_vec, prop = depigm_ageprev_vec)
depigm_ageprev_df$prev <- depigm_ageprev_df$prop * 100
depigm_ageprev_df$sequela <- as.factor("depigmentation")
depigm_ageprev_df$agecat <- c("0 - 1", "2 - 4", "5 - 9", "10 - 19", "20 - 29", "30 - 49", "50 - 80")
depigm_ageprev_df$agecat2 <- ordered(depigm_ageprev_df$agecat, levels = c("0 - 1", "2 - 4", "5 - 9", "10 - 19", "20 - 29", "30 - 49", "50 - 80"))

ggplot(data = depigm_ageprev_df, aes(x = agecat2, y =prev, group = sequela))+
  geom_point()+
  geom_line()

# combine #
modelled_age_sequela_prev_df <- rbind(SI_ageprev_df, RSD_ageprev_df, Atrp_ageprev_df, HG_ageprev_df, depigm_ageprev_df)

saveRDS(modelled_age_sequela_prev_df, "modelled_age_sequela_prev_df.rds")

saveRDS(modelled_age_sequela_prev_df, "modelled_age_sequela_prev_df_test6.rds")
saveRDS(modelled_age_sequela_prev_df, "modelled_age_sequela_prev_df_70mfprev.rds")
saveRDS(modelled_age_sequela_prev_df, "modelled_age_sequela_prev_df_80mfprev.rds")
saveRDS(modelled_age_sequela_prev_df, "modelled_age_sequela_prev_df_chesnais.rds")

# ===================== #
# age - mf prevalence   #

all.mats.temp_out <- output_equilibrium_morbidity$all_equilibrium_outputs[[1]] # to get mf out

mf.start <- 7
mf.end <- 6 + 21

mf_all_vec <- rowSums(all.mats.temp_out[, mf.start : mf.end])

mf_age_df <- data.frame(age = all.mats.temp_out[,2], mf_count = mf_all_vec)

mf_age_df$mf_pos <- ifelse(mf_age_df$mf_count > 0, 1, 0)

require(dplyr)
mf_age_df <- mf_age_df %>%
  mutate(
    #create categories
    age_cat = dplyr::case_when(
      age >= 0 & age < 2 ~ "0 - 1",
      age >= 2 & age < 5 ~ "2 - 4",
      age >= 5 & age < 10 ~ "5 - 9",
      age >= 10 & age < 20 ~ "10 - 19",
      age >= 20 & age < 30 ~ "20 - 29",
      age >= 30 & age < 50 ~ "30 - 49",
      age >= 50 & age <= 80 ~ "50 - 80",
      age > 80 ~ "> 80"
    ),
    # convert to factor
    age_cat = factor(
      age_cat,
      level = c("0 - 1","2 - 4", "5 - 9", "10 - 19", "20 - 29", "30 - 49", "50 - 80", "> 80")
    )
  )

modelled_mfprev_age_prop <- mf_age_df %>%
  group_by(age_cat) %>%
  summarise(x = sum(mf_pos),
            n = n(),
            prop = sum(mf_pos)/n())

# mfprev_df <- data.frame(age_cat = mfprev_age_prop$age_cat, x = mfprev_age_prop$x, n = mfprev_age_prop$n, prop = mfprev_age_prop$prop)

saveRDS(modelled_mfprev_age_prop, "modelled_age_mfprev_df.rds")

saveRDS(modelled_mfprev_age_prop, "modelled_age_mfprev_df_test6.rds")
saveRDS(modelled_mfprev_age_prop, "modelled_age_mfprev_df_70mfprev.rds")
saveRDS(modelled_mfprev_age_prop, "modelled_age_mfprev_df_80mfprev.rds")
saveRDS(modelled_mfprev_age_prop, "modelled_age_mfprev_df_chesnais.rds")


# ========================================================================================================== #
#                    Loop through to specific iteration in debugonce                                         #

# Define the total number of iterations
total_iterations <- 100

# Define the specific iteration to stop at
target_iteration <- 59

# Loop through the iterations
for (i in 1:total_iterations) {
  if (i == target_iteration) {
    # Set a debug flag on a specific function
    debugonce(ep.equi.sim)

    # Break the loop and start debugging
    break
  }

  # Perform your iteration logic here
  print(paste("Iteration:", i))

  output_equilibrium_morbidity <- ep.equi.sim(time.its = timesteps,
                                              ABR = ABR.in,
                                              N.in = 440,
                                              treat.int = trt.int,
                                              treat.prob = 0.65,
                                              give.treat = give.treat.in,
                                              treat.start = treat.strt,
                                              treat.stop = treat.stp,
                                              pnc = 0.05,
                                              min.mont.age = 5,
                                              vector.control.strt = NA, #
                                              delta.hz.in = 0.186,
                                              delta.hinf.in = 0.003,
                                              c.h.in = 0.005,
                                              gam.dis.in = 0.3,
                                              run_equilibrium = TRUE,
                                              equilibrium,
                                              print_progress = TRUE,
                                              morbidity_module2 = "YES")
}










