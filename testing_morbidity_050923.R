# ========================================================================================================================== #
#                   TESTING MORBIDITY (SKIN DISEASE) --- Taking EQ inputs -- 30/06/23                                                #
# ========================================================================================================================== #

devtools::load_all()

#length of simulation in years
# timesteps = 200
# timesteps = 130
# timesteps = 80

timesteps = 100

#should treatment be given, when and how often
give.treat.in = 0
treat.strt = 1; treat.stp = 16
trt.int = 1

#annual biting rate, which determines infection prevalence (60% microfilarae prevalence)
ABR.in = 607 # 50% mf prev
ABR.in = 1082 # 60% mf prev
ABR.in = 2297 # 70% mf prev
ABR.in = 6969 # 80% mf prev
ABR.in = 41922 # chesnais ABR

# output_equilibrium_morbidity <- ep.equi.sim(time.its = timesteps,
#                                             ABR = ABR.in,
#                                             N.in = 440,
#                                             treat.int = trt.int,
#                                             treat.prob = 0.65,
#                                             give.treat = give.treat.in,
#                                             treat.start = treat.strt,
#                                             treat.stop = treat.stp,
#                                             pnc = 0.05,
#                                             min.mont.age = 5,
#                                             vector.control.strt = NA, #
#                                             delta.hz.in = 0.186,
#                                             delta.hinf.in = 0.003,
#                                             c.h.in = 0.005,
#                                             gam.dis.in = 0.3,
#                                             run_equilibrium = TRUE,
#                                             equilibrium,
#                                             print_progress = TRUE,
#                                             morbidity_module = "YES")

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
                                            morbidity_module = "YES")
names(output_equilibrium_morbidity)

#output_equilibrium_morbidity_mf50prev <- output_equilibrium_morbidity
#output_equilibrium_morbidity_mf70prev <- output_equilibrium_morbidity
#output_equilibrium_morbidity_mf80prev <- output_equilibrium_morbidity
#output_equilibrium_morbidity_chesnaisABR <- output_equilibrium_morbidity
# tme <- seq(1, 130*366-1)/366
# tme <- seq(1, 80*366-1)/366
# tme <- seq(1, 200*366-1)/366

# ========== #
# mf outputs #
tme <- seq(1, 100*366-1)/366

plot(tme, output_equilibrium_morbidity$mf_prev, type = 'l', xlab = 'time (years)', ylab = 'microfilarial prevalence', ylim = c(0, 1))
tail(output_equilibrium_morbidity$mf_prev)
mean(tail(output_equilibrium_morbidity$mf_prev))

plot(tme, output_equilibrium_morbidity$mf_intens, type = 'l', xlab = 'time (years)', ylab = 'mean microfilarial intensity')

plot(tme, output_equilibrium_morbidity$L3, type = 'l', xlab = 'time (years)', ylab = 'mean L3 in fly population')

# ============ #
# skin disease #
#tme2 <- seq(1, 130*366)/366

# tme2 <- seq(1, 200*366)/366
# tme2 <- seq(1, 130*366)/366
tme2 <- seq(1, 100*366)/366

output_equilibrium_morbidity_test1 <- output_equilibrium_morbidity

plot(tme2, output_equilibrium_morbidity$severe_itch_prev, type = 'l', xlab = 'time (years)', ylab = 'Severe itch prevalence', ylim = c(0, 1))
tail(output_equilibrium_morbidity$severe_itch_prev)
mean(tail(output_equilibrium_morbidity$severe_itch_prev))

plot(tme2, output_equilibrium_morbidity$RSD_prev, type = 'l', xlab = 'time (years)', ylab = 'Reactive skin disease prevalence', ylim = c(0, 1))
tail(output_equilibrium_morbidity$RSD_prev)
mean(tail(output_equilibrium_morbidity$RSD_prev))

plot(tme2, output_equilibrium_morbidity$atrophy_prev, type = 'l', xlab = 'time (years)', ylab = 'Skin atrophy prevalence', ylim = c(0, 1))
tail(output_equilibrium_morbidity$atrophy_prev)
mean(tail(output_equilibrium_morbidity$atrophy_prev))

plot(tme2, output_equilibrium_morbidity$atrophy_prev, type = 'l', xlab = 'time (years)', ylab = 'atrophy prev',
     ylim = c(0, 0.05), xlim = c(0, 5))

plot(tme2, output_equilibrium_morbidity$hanging_groin_prev, type = 'l', xlab = 'time (years)', ylab = 'Hanging groin prevalence', ylim = c(0, 1))
tail(output_equilibrium_morbidity$hanging_groin_prev)
mean(tail(output_equilibrium_morbidity$hanging_groin_prev))

plot(tme2, output_equilibrium_morbidity$hanging_groin_prev, type = 'l', xlab = 'time (years)', ylab = 'HG prev',
     ylim = c(0, 0.05), xlim = c(0, 5))

plot(tme2, output_equilibrium_morbidity$depigmentation_prev, type = 'l', xlab = 'time (years)', ylab = 'depigmentation prevalence', ylim = c(0, 1))
tail(output_equilibrium_morbidity$depigmentation_prev)
mean(tail(output_equilibrium_morbidity$depigmentation_prev))

plot(tme2, output_equilibrium_morbidity$depigmentation_prev, type = 'l', xlab = 'time (years)', ylab = 'depigm prev',
     ylim = c(0, 0.05), xlim = c(0, 5))


# ========================= #
#     Eye disease           #

#tme3 <- seq(1, 130*366)/366
#tme3 <- seq(1, 132*366 - 1/366)/366

tme3 <- seq(1, 102*366 - 1/366)/366

plot(tme3, output_equilibrium_morbidity$blindness_prev, type = 'l', xlab = 'time (years)', ylab = 'blindness prev (lagged mat)', ylim = c(0, 1))
tail(output_equilibrium_morbidity$blindness_prev)
mean(tail(output_equilibrium_morbidity$blindness_prev))

plot(tme3, output_equilibrium_morbidity$blindness_prev, type = 'l', xlab = 'time (years)', ylab = 'blindness prev (lagged mat)',
     ylim = c(0, 0.05), xlim = c(0, 5)) # check 2 yr lag at start

plot(tme3, output_equilibrium_morbidity$visual_impairment_prev, type = 'l', xlab = 'time (years)', ylab = 'visual impairement prev', ylim = c(0, 1))
tail(output_equilibrium_morbidity$visual_impairment_prev)
mean(tail(output_equilibrium_morbidity$visual_impairment_prev))

plot(tme3, output_equilibrium_morbidity$visual_impairment_prev, type = 'l', xlab = 'time (years)', ylab = 'visual impairement prev',
     ylim = c(0, 0.05), xlim = c(0, 5)) # check 2 yr lag at start

# chop off final 2 years of lagged vec as two years longer than lagged mat

length_req_2yrs <- length(seq(0, 2 - 2/366, by = 1/366))
blind_lagged_vec <- head(output_equilibrium_morbidity$blindness_prev1, -length_req_2yrs)

plot(tme3, blind_lagged_vec, type = 'l', xlab = 'time (years)', ylab = 'blindness prev (lagged vec)', ylim = c(0, 1))
tail(blind_lagged_vec)
mean(tail(blind_lagged_vec))

plot(tme3, blind_lagged_vec, type = 'l', xlab = 'time (years)', ylab = 'blindness prev (lagged vec)',
     ylim = c(0, 0.05), xlim = c(0, 5)) # check 2 yr lag at start

# ============================== #
#   check blindness age-prev     #

names(output_equilibrium_morbidity)

ageprev_out <- output_equilibrium_morbidity$eye_morbidity_ageprev_out

# blind modelled age-prev #
blind_prev_0_1 <- mean(tail(ageprev_out$blind_prev0_1))
blind_prev_2_4 <- mean(tail(ageprev_out$blind_prev2_4))
blind_prev_5_9 <- mean(tail(ageprev_out$blind_prev5_9))
blind_prev_10_19 <- mean(tail(ageprev_out$blind_prev10_19))
blind_prev_20_29 <- mean(tail(ageprev_out$blind_prev20_29))
blind_prev_30_49 <- mean(tail(ageprev_out$blind_prev30_49))
blind_prev_50_80 <- mean(tail(ageprev_out$blind_prev50_80))

blind_prev_ageprev_vec <- c(blind_prev_0_1, blind_prev_2_4, blind_prev_5_9, blind_prev_10_19, blind_prev_20_29, blind_prev_30_49, blind_prev_50_80)
age_vec <- c(1, 3, 7.5, 15, 25, 40, 65)
blind_prev_ageprev_df <- data.frame(age = age_vec, prop = blind_prev_ageprev_vec)
blind_prev_ageprev_df$prev <- blind_prev_ageprev_df$prop * 100
blind_prev_ageprev_df$sequela <- as.factor("blindness")
blind_prev_ageprev_df$agecat <- c("0 - 1", "2 - 4", "5 - 9", "10 - 19", "20 - 29", "30 - 49", "50 - 80")
blind_prev_ageprev_df$agecat2 <- ordered(blind_prev_ageprev_df$agecat, levels = c("0 - 1", "2 - 4", "5 - 9", "10 - 19", "20 - 29", "30 - 49", "50 - 80"))

require(ggplot2)
ggplot(data = blind_prev_ageprev_df, aes(x = agecat2, y =prev, group = sequela))+
  geom_point()+
  geom_line()

# visual impairement modelled age-prev #
visual_imp_prev_0_1 <- mean(tail(ageprev_out$visual_imp_prev0_1))
visual_imp_prev_2_4 <- mean(tail(ageprev_out$visual_imp_prev2_4))
visual_imp_prev_5_9 <- mean(tail(ageprev_out$visual_imp_prev5_9))
visual_imp_prev_10_19 <- mean(tail(ageprev_out$visual_imp_prev10_19))
visual_imp_prev_20_29 <- mean(tail(ageprev_out$visual_imp_prev20_29))
visual_imp_prev_30_49 <- mean(tail(ageprev_out$visual_imp_prev30_49))
visual_imp_prev_50_80 <- mean(tail(ageprev_out$visual_imp_prev50_80))

visual_imp_prev_ageprev_vec <- c(visual_imp_prev_0_1, visual_imp_prev_2_4, visual_imp_prev_5_9,
                                 visual_imp_prev_10_19, visual_imp_prev_20_29, visual_imp_prev_30_49, visual_imp_prev_50_80)
age_vec <- c(1, 3, 7.5, 15, 25, 40, 65)
visual_imp_prev_ageprev_df <- data.frame(age = age_vec, prop = visual_imp_prev_ageprev_vec)
visual_imp_prev_ageprev_df$prev <- visual_imp_prev_ageprev_df$prop * 100
visual_imp_prev_ageprev_df$sequela <- as.factor("visual impairment")
visual_imp_prev_ageprev_df$agecat <- c("0 - 1", "2 - 4", "5 - 9", "10 - 19", "20 - 29", "30 - 49", "50 - 80")
visual_imp_prev_ageprev_df$agecat2 <- ordered(visual_imp_prev_ageprev_df$agecat, levels = c("0 - 1", "2 - 4", "5 - 9", "10 - 19", "20 - 29", "30 - 49", "50 - 80"))

require(ggplot2)
ggplot(data = visual_imp_prev_ageprev_df, aes(x = agecat2, y =prev, group = sequela))+
  geom_point()+
  geom_line()


eye_dis_df <- rbind(blind_prev_ageprev_df, visual_imp_prev_ageprev_df)

ggplot(data = eye_dis_df, aes(x = agecat2, y =prev, group = sequela))+
  geom_point()+
  geom_line()+
  facet_wrap(~sequela)

# ==========================================================================================================#
#        Running with equilibria inputs - will not work as blindness needs to be updated to take EQ inputs  #
# ==========================================================================================================#

devtools::load_all()

timesteps = 30
treat.strt = 1; treat.stp = 26 #if model run stops in year 30, the last round is at the beginning of year 25
gv.trt = 0
trt.int = 1

output_equilibrium <- output_equilibrium_morbidity

morbidity_prev_eq <- list(output_equilibrium$severe_itch_prev, output_equilibrium$RSD_prev,
                       output_equilibrium$atrophy_prev, output_equilibrium$hanging_groin_prev,
                       output_equilibrium$hanging_groin_prev)

output_equilibrium_morbidity2 <- ep.equi.sim(time.its = timesteps,
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
                                             run_equilibrium = FALSE,
                                             equilibrium = output_equilibrium[[4]],
                                             morbidity_eq = output_equilibrium[[13]],
                                             morbidity_ageprev_eq = output_equilibrium[[10]],
                                             morbidity_prev_eq = morbidity_prev_eq,
                                             print_progress = TRUE,
                                             morbidity_module2 = "YES",
                                             mf_ageprev_eq = output_equilibrium[[14]],
                                             mf_agintens_eq = output_equilibrium[[15]])


# morbidity_ageprev_eq_test <- output_equilibrium[[10]]
# morb.test <- output_equilibrium$morbidity.outputs[[1]]
# morb.test2 <- output_equilibrium[[13]]


tme <- seq(1, 30*366)/366

plot(tme, output_equilibrium_morbidity2$mf_prev, type = 'l', xlab = 'time (years)', ylab = 'microfilarial prevalence', ylim = c(0, 1))
tail(output_equilibrium_morbidity2$mf_prev)
mean(tail(output_equilibrium_morbidity2$mf_prev))

plot(tme, output_equilibrium_morbidity$mf_intens, type = 'l', xlab = 'time (years)', ylab = 'mean microfilarial intensity')

plot(tme, output_equilibrium_morbidity$L3, type = 'l', xlab = 'time (years)', ylab = 'mean L3 in fly population')

# ============ #
#tme2 <- seq(1, 130*366)/366

tme2 <- seq(1, 160*366-1)/366 # 130 yrs + 30 yrs

plot(tme2, output_equilibrium_morbidity2$severe_itch_prev, type = 'l', xlab = 'time (years)', ylab = 'Severe itch prevalence', ylim = c(0, 1))
tail(output_equilibrium_morbidity2$severe_itch_prev)
mean(tail(output_equilibrium_morbidity2$severe_itch_prev))

plot(tme2, output_equilibrium_morbidity2$RSD_prev, type = 'l', xlab = 'time (years)', ylab = 'Reactive skin disease prevalence', ylim = c(0, 1))
tail(output_equilibrium_morbidity2$RSD_prev)
mean(tail(output_equilibrium_morbidity2$RSD_prev))

plot(tme2, output_equilibrium_morbidity2$atrophy_prev, type = 'l', xlab = 'time (years)', ylab = 'Skin atrophy prevalence', ylim = c(0, 1))
tail(output_equilibrium_morbidity$atrophy_prev)
mean(tail(output_equilibrium_morbidity$atrophy_prev))

plot(tme2, output_equilibrium_morbidity2$hanging_groin_prev, type = 'l', xlab = 'time (years)', ylab = 'Hanging groin prevalence', ylim = c(0, 1))
tail(output_equilibrium_morbidity2$hanging_groin_prev)
mean(tail(output_equilibrium_morbidity2$hanging_groin_prev))

plot(tme2, output_equilibrium_morbidity2$depigmentation_prev, type = 'l', xlab = 'time (years)', ylab = 'depigmentation prevalence', ylim = c(0, 1))
tail(output_equilibrium_morbidity2$depigmentation_prev)
mean(tail(output_equilibrium_morbidity2$depigmentation_prev))

