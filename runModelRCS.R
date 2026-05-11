output_prefix = Sys.getenv("OUTPUT_PREFIX")
output_folder = Sys.getenv("OUTPUT_FOLDER")
iter <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
set.seed(iter + (iter*3758))

num_iters_per_paramset = 200
DT_in <- 1/366

kEs = c(rep(0.2, num_iters/2), rep(0.3, num_iters/2))
kE <- kEs[iter]


if (kE == 0.2) {
  ABR_in <- 73
} else {
  ABR_in <- 176
}

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
  min.mont.age = 5,
  vector.control.strt = NA,
  gam.dis.in = kE,
  run_equilibrium = TRUE,
  morbidity_module = "YES",
  equilibrium = NA,
  print_progress = TRUE,
  prob_serorevert_fast = 0.5
)

params <- list(ABR_in, kE)
names(params) <- c('ABR', 'Ke')
output <- append(output, params)
tmp_dir <- Sys.getenv("TMPDIR")
saveRDS(output, paste(tmp_dir, "/", output_folder, "/",
                      output_prefix , "_", iter, ".rds", sep = ""))