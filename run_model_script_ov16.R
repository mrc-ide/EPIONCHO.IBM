library(devtools)
library(dplyr)
library(ggplot2)
devtools::load_all()

# First running to equilibrium
# ==================================================== #

#iter <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
#iter <- 1000

#set.seed(iter + (iter*3758))

DT.in <- 1/366
timesteps = 200
give.treat.in = 0
treat.strt = 1; treat.stp = 16
trt.int = 1


#ABR_in <- 41922 #annual biting rate
brz <- c(1500, 2000, 3000)

abr_vec <- rep(brz, 400)


ABR.in <- abr_vec[runif(1, 1, length(abr_vec))]
ABR.in <- 1500

output_equilibrium <-  ep.equi.sim(time.its = timesteps,
                       ABR = ABR.in,
                       treat.int = trt.int,
                       treat.prob = 0.80,
                       give.treat = give.treat.in,
                       treat.start = treat.strt,
                       treat.stop = treat.stp,
                       pnc = 0.01,
                       min.mont.age = 5,
                       vector.control.strt = NA,
                       delta.hz.in = 0.186,
                       delta.hinf.in = 0.003,
                       c.h.in = 0.005,
                       gam.dis.in = 0.3,
                       run_equilibrium=TRUE,
                       print_progress=TRUE)


# Now add treatment

treat.len <- 18 #treatment duration in years
treat.strt.yrs <- 10
yrs.post.treat <- 30

treat.strt  = round(treat.strt.yrs / (DT.in )); treat.stp = treat.strt + round(treat.len / (DT.in ))

treat.strt = treat.strt.yrs; treat.stp = treat.strt + treat.len
timesteps = treat.stp + round(yrs.post.treat / (DT.in )) #final duration
timesteps = treat.stp + yrs.post.treat #final duration

give.treat.in = 1 # set to 0 to turn off treatment
trt.int = 1 #treatment interval (years, 0.5 gives biannual)

numRuns <- 100
allOutputs <- data.frame(matrix(ncol=9))
colnames(allOutputs) <- c("age", "sex", "ov16_pos", "mf_prev", "age_pre", "sex_pre", "ov16_pos_pre", "mf_prev_pre", "run_num")

#---- Multi Run ------#

for(i in 1:numRuns) {
  print(paste("Run Number:", i))
  output <- ep.equi.sim(time.its = timesteps,
                        ABR = ABR.in,
                        treat.int = trt.int,
                        treat.prob = 0.80,
                        give.treat = give.treat.in,
                        treat.start = treat.strt,
                        treat.stop = treat.stp,
                        treat.timing = NA,
                        pnc = 0.01,
                        min.mont.age = 5,
                        vector.control.strt = NA,
                        delta.hz.in = 0.186,
                        delta.hinf.in = 0.003,
                        c.h.in = 0.005,
                        gam.dis.in = 0.3,
                        run_equilibrium = FALSE,
                        equilibrium = output_equilibrium[[4]],
                        print_progress=TRUE,
                        calc_ov16 = TRUE)

  age <- output$all_infection_burdens[,2]

  age_pre <- output$all_infection_burdens_pre_treatment[,2]

  sex <- ifelse(output$all_infection_burdens[,3]==1, "Male", "Female")

  sex_pre <- ifelse(output$all_infection_burdens_pre_treatment[,3]==1, "Male", "Female")

  mf_prev <- output$mf_indv_prevalence

  mf_prev_pre <- output$mf_indv_prevalence_pre_treatment

  ov16_seropos <- output$ov16_seropositive

  ov16_seropos_pre <- output$ov16_seropositive_pre_treatment

  tmpNumRows <- length(age)
  startIndex <- 1+tmpNumRows*(i-1)
  endIndex <- tmpNumRows*i
  allOutputs[startIndex:endIndex,-9] <- list(age, sex, ov16_seropos, mf_prev, age_pre, sex_pre, ov16_seropos_pre, mf_prev_pre)
  allOutputs[startIndex:endIndex,9] <- i
}

write.csv(allOutputs, "data.csv")


#---- Single Run -----#

output <- ep.equi.sim(time.its = timesteps,
                       ABR = ABR.in,
                       treat.int = trt.int,
                       treat.prob = 0.80,
                       give.treat = give.treat.in,
                       treat.start = treat.strt,
                       treat.stop = treat.stp,
                       treat.timing = NA,
                       pnc = 0.01,
                       min.mont.age = 5,
                       vector.control.strt = NA,
                       delta.hz.in = 0.186,
                       delta.hinf.in = 0.003,
                       c.h.in = 0.005,
                       gam.dis.in = 0.3,
                       run_equilibrium = FALSE,
                       equilibrium = output_equilibrium[[4]],
                       print_progress=TRUE,
                       calc_ov16 = TRUE)

age <- output$all_infection_burdens[,2]

age_pre <- output$all_infection_burdens_pre_treatment[,2]

sex <- ifelse(output$all_infection_burdens[,3]==1, "Male", "Female")

sex_pre <- ifelse(output$all_infection_burdens_pre_treatment[,3]==1, "Male", "Female")

mf_prev <- output[[2]] # todo

ov16_prev <- output$ov16_seropositive

ov16_prev_pre <- output$ov16_seropositive_pre_treatment

i <- 1
tmpNumRows <- length(age)
startIndex <- 1+tmpNumRows*(i-1)
endIndex <- tmpNumRows*i
allOutputs[startIndex:endIndex,] <- list(age, sex, ov16_prev)

# viz


allOutputs <- allOutputs %>% mutate(new_age = round(age/10)*10,
                                    age_groups = case_when(
                                      #age < 5 & as.integer(age) %% 2 == 1 ~ as.integer(age),
                                      #age < 5 & as.integer(age) %% 2 == 0 ~ age + 1,
                                      age < 1 ~ 1,
                                      age < 3 ~ 3,
                                      age < 5 ~ 5,
                                      age < 10 ~ 10,
                                      age < 15 ~ 15,
                                      age < 20 ~ 20,
                                      age < 30 ~ 30,
                                      age < 40 ~ 40,
                                      age < 50 ~ 50,
                                      age < 60 ~ 60,
                                      age < 70 ~ 70,
                                      TRUE ~ 80
                                    ),
                                    age_groups_pre = case_when(
                                      age_pre < 1 ~ 1,
                                      age_pre < 3 ~ 3,
                                      age_pre < 5 ~ 5,
                                      age_pre < 10 ~ 10,
                                      age_pre < 15 ~ 15,
                                      age_pre < 20 ~ 20,
                                      age_pre < 30 ~ 30,
                                      age_pre < 40 ~ 40,
                                      age_pre < 50 ~ 50,
                                      age_pre < 60 ~ 60,
                                      age_pre < 70 ~ 70,
                                      TRUE ~ 80
                                    ))

tmpDf <- allOutputs %>% dplyr::group_by(age_groups, sex) %>% dplyr::summarise(ov16_prev=mean(ov16_pos), mf_prev=mean(mf_prev)) %>% as.data.frame() #%>% tidyr::pivot_longer(c(ov16_prev, mf_prev), names_to="treatment", values_to="ov16_prev") %>% as.data.frame()

tmpDf2 <- allOutputs %>% dplyr::group_by(age_groups_pre, sex_pre) %>% dplyr::summarise(ov16_prev_pre=mean(ov16_pos_pre), mf_prev_pre=mean(mf_prev_pre)) %>% as.data.frame() #%>% tidyr::pivot_longer(c(mf_prev, mf_prev_pre), names_to="treatment", values_to="mf_prev") %>% as.data.frame()

ggplot() +
  geom_smooth(aes(x=age_groups, y=ov16_prev, color=sex, linetype='Post Treatment'), se=F, data=tmpDf) +
  geom_smooth(aes(x=age_groups_pre, y=ov16_prev_pre, color=sex_pre, linetype='Pre Treatment'), se=F, data=tmpDf2) +
  scale_linetype_manual(values=c("solid", "dashed"))

ggplot() +
  geom_smooth(aes(x=age_groups, y=mf_prev, color=sex, linetype='Post Treatment'), se=F, data=tmpDf) +
  geom_smooth(aes(x=age_groups_pre, y=mf_prev_pre, color=sex_pre, linetype='Pre Treatment'), se=F, data=tmpDf2) +
  scale_linetype_manual(values=c("solid", "dashed"))



allOutputs %>% dplyr::group_by(age_groups, sex) %>% dplyr::summarise(ov16_prev=mean(ov16_pos), ov16_prev_pre=mean(ov16_pos_pre)) %>% as.data.frame() %>%
  ggplot(aes(x=age_groups, y=ov16_prev, color=sex)) +
  geom_point() +
  geom_line()

ov16_repl_graph

ggsave("Second_Attempt_ov16.png", ov16_repl_graph)

# ov16_repl_graph <- allOutputs %>% dplyr::group_by(new_age, sex) %>% dplyr::summarise(ov16_prev=mean(ov16_pos)) %>% as.data.frame() %>%
#   ggplot(aes(x=new_age, y=ov16_prev, color=sex)) +
#   geom_point() +
#   geom_line()
# ov16_repl_graph
