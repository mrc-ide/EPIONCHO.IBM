output_prefix = Sys.getenv("OUTPUT_PREFIX")
output_folder = Sys.getenv("OUTPUT_FOLDER")
iter <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
set.seed(iter + (iter*3758))

iters_per_paramset = 200

DT_in <- 1/366

abrs_to_test = c(50, 70, 90, 100, 130, 160, 200, 300, 400, 500, 600, 800, 1000, 1250, 1500, 2000, 3000, 4000, 5000, 6000, 10000, 20000, 50000)



kEs = rep(c(rep(0.2, iters_per_paramset), rep(0.3, iters_per_paramset), rep(0.4, iters_per_paramset)), each=length(abrs_to_test))
kE <- kEs[iter]

# 200x50, 200x70, ... 200x50000, then starts again at 200x50... 3 sets of 200 for each ABR
abrs_to_test_iter = rep(rep(abrs_to_test, each=iters_per_paramset), length(unique(kEs)))

# no treatment to be given
give_treat = 0
treatment_interval = 1
timesteps_before_treatment = 100
treatment_length = 0
years_post_treatment = 0
treatment_start_year = timesteps_before_treatment
treatment_stop_year = treatment_start_year + treatment_length
total_timesteps = treatment_stop_year + years_post_treatment #final duration

output <- ep.equi.sim(
  time.its = total_timesteps,
  ABR = ABR_in,
  N.in = 400,
  treat.int = treatment_interval,
  treat.prob = 0.65,
  give.treat = give_treat,
  treat.start = treatment_start_year,
  treat.stop = treatment_stop_year,
  pnc = 0.05,
  min.mont.age = 0,
  vector.control.strt = NA,
  gam.dis.in = kE,
  run_equilibrium = TRUE,
  morbidity_module = "YES",
  equilibrium = NA,
  print_progress = TRUE,
  prob_serorevert_fast = 0.5,
  output_age_groups = list(c(0, 5), c(5, 10), c(10, 15), c(15, 20), c(20, 30), c(30, 40) c(40, 50), c(50, 60), c(60, 70), c(70, 81))
)

params <- list(ABR_in, kE)
names(params) <- c('ABR', 'Ke')
output <- append(output, params)
tmp_dir <- Sys.getenv("TMPDIR")
saveRDS(output, paste(tmp_dir, "/", output_folder, "/",
                      output_prefix , "_", iter, ".rds", sep = ""))