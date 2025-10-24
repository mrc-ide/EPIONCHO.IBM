iter <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
set.seed(iter + (iter*3758))

#Biting Rates Vector (Hyper.abr1-4, Meso.abr, Hypo.abr)
#Hyper.abr1 to 4 ----------------------------------------------------------------------


ABR_values <- rep(c(290, 615, 2200, 60000), each=100)
tlength <- c(30)
treatz <- c(1) 

tot.parms <- expand.grid(ABR_values, treatz)

#coverage

os.cov <- function(all.dt, pncomp, covrg, coi) 
  
{                   
  pop.ages <- all.dt[,2]
  
  iny <- which(pop.ages < 5 | all.dt[,1] == 1)
  
  nc.age <- length(iny) / length(pop.ages)
  
  covrg <- covrg / (1 - nc.age)
  
  out.cov <- rep(covrg, length(pop.ages))
  
  out.cov[iny] <- 0
  
  f.cov <- rep(0, 500)
  
  r.nums <- runif(500, 0, 1)
  
  inds.c <- which(r.nums < out.cov)
  
  f.cov[inds.c] <- 1
  
  #print('final cov')
  #print(sum(f.cov) / 500)
  
  return(f.cov)
  
}



#function rotate matrix, used in mf function 
rotate <- function(x) t(apply(x, 2, rev))

change.micro <- function(dat, num.comps, mf.cpt, num.mf.comps, ws, DT, time.each.comp, mu.rates.mf, fec.rates, mf.move.rate,
                         up, kap, iteration, treat.vec, give.treat, treat.start)
  
{
  N <- length(dat[,1])
  #indexes for fertile worms (to use in production of mf)
  fert.worms.start <-  ws + num.comps*2 
  fert.worms.end <-  (ws-1) + num.comps*3
  
  #indexes to check if there are males (males start is just 'ws')
  mal.worms.end <- (ws-1) + num.comps
  mf.mu <- rep(mu.rates.mf[mf.cpt], N)
  fert.worms <- dat[, fert.worms.start:fert.worms.end]
  
  if(give.treat == 1 & iteration >= treat.start)
  {
    tao <- ((iteration-1)*DT) - treat.vec #tao is zero if treatment has been given at this timestep
    
    mu.mf.prime <- ((tao + up) ^ (-kap))
    
    mu.mf.prime[which(is.na(mu.mf.prime) == TRUE)] <- 0
    
    mf.mu <- mf.mu + mu.mf.prime
    
  }
  
  if(mf.cpt == 1)
  {
    mp <- rep(0, N)
    
    inds.fec <- which(rowSums(dat[, ws : mal.worms.end]) > 0); mp[inds.fec] <- 1     #need to check there is more than one male
    
    k1 <- derivmf.one(fert.worms = fert.worms, mf.in = dat[, 6 + mf.cpt], ep.in = fec.rates, mf.mort = mf.mu, mf.move = mf.move.rate, mp = mp, k.in = 0)  #fert worms and epin are vectors
    k2 <- derivmf.one(fert.worms = fert.worms, mf.in = dat[, 6 + mf.cpt], ep.in = fec.rates, mf.mort = mf.mu, mf.move = mf.move.rate, mp = mp, k.in = DT*k1/2) 
    k3 <- derivmf.one(fert.worms = fert.worms, mf.in = dat[, 6 + mf.cpt], ep.in = fec.rates, mf.mort = mf.mu, mf.move = mf.move.rate, mp = mp, k.in = DT*k2/2)  
    k4 <- derivmf.one(fert.worms = fert.worms, mf.in = dat[, 6 + mf.cpt], ep.in = fec.rates, mf.mort = mf.mu, mf.move = mf.move.rate, mp = mp, k.in = DT*k3) 
    
    out <- dat[, 6 + mf.cpt] + DT/6*(k1+2*k2+2*k3+k4)
    
  }
  
  
  if(mf.cpt > 1)
  {
    k1 <- derivmf.rest(mf.in = dat[, 6 + mf.cpt], mf.mort = mf.mu, mf.move = mf.move.rate, mf.comp.minus.one = dat[, 6 + mf.cpt - 1], k.in = 0) 
    k2 <- derivmf.rest(mf.in = dat[, 6 + mf.cpt], mf.mort = mf.mu, mf.move = mf.move.rate, mf.comp.minus.one = dat[, 6 + mf.cpt - 1], k.in = DT*k1/2) 
    k3 <- derivmf.rest(mf.in = dat[, 6 + mf.cpt], mf.mort = mf.mu, mf.move = mf.move.rate, mf.comp.minus.one = dat[, 6 + mf.cpt - 1], k.in = DT*k2/2) 
    k4 <- derivmf.rest(mf.in = dat[, 6 + mf.cpt], mf.mort = mf.mu, mf.move = mf.move.rate, mf.comp.minus.one = dat[, 6 + mf.cpt - 1], k.in = DT*k3) 
    
    out <- dat[, 6 + mf.cpt] + DT/6*(k1+2*k2+2*k3+k4)
    
  }  
  
  
  return(out)
}

derivmf.one <- function(fert.worms, mf.in, ep.in, mf.mort, mf.move, mp, k.in)  #fert worms and epin are vectors
{
  
  new.in <- (rotate(fert.worms)*ep.in) #need to rotate matrix to each column is multiplied by respective fecundity rate, not each row
  new.in <- rotate(rotate(rotate(new.in)))
  new.in <- rowSums(new.in)
  
  #out <- mp * new.in - mf.mort*(mf.in + k.in) - mf.move * (mf.in + k.in) 
  
  
  mort.temp <- mf.mort*(mf.in + k.in)
  move.temp <- mf.move * (mf.in + k.in) 
  
  mort.temp[which(mort.temp < 0)] <- 0
  move.temp [which(move.temp < 0)] <- 0
  
  if(length(which(mf.mort*(mf.in + k.in) < 0)) >0) {print('WARNING MF NEGATIVE1')}
  if(length(which(mf.move * (mf.in + k.in) < 0)) >0) {print('WARNING MF NEGATIVE2')}
  
  
  out <- mp * new.in - mort.temp - move.temp
  
  # if(length(which(out < 0)) >0) {print('WARNING MF out')}
  # 
  # out[which(out < 0)] <- 0
  
  return(out)
}


derivmf.rest <- function(mf.in, mf.mort, mf.move, mf.comp.minus.one, k.in) 
{
  #out <- mf.comp.minus.one*mf.move - mf.mort*(mf.in + k.in) - mf.move * (mf.in + k.in)
  
  move.last <- mf.comp.minus.one*mf.move
  mort.temp <- mf.mort * (mf.in + k.in)
  move.temp <- mf.move * (mf.in + k.in)
  
  move.last[which(move.last < 0)] <- 0
  mort.temp[which(mort.temp < 0)] <- 0
  move.temp [which(move.temp < 0)] <- 0
  
  if(length(which(mf.mort*(mf.in + k.in) < 0)) >0) {print('WARNING MF NEGATIVE3')}
  if(length(which(mf.move * (mf.in + k.in) < 0)) >0) {print('WARNING MF NEGATIVE4')}
  
  
  out <- move.last - mort.temp - move.temp
  
  # if(length(which(out < 0)) >0) {print('WARNING MF out')}
  # 
  # out[which(out < 0)] <- 0
  
  return(out)
}


#prop of L3 larvae developing into adult worms in one human, expos = total exposure for an individual

delta.h <- function(delta.hz, delta.hinf, c.h, L3, m , beta, expos)
  
{
  out <- (delta.hz + delta.hinf * c.h * m * beta *  L3 * expos) / (1 + c.h * m * beta * L3 * expos) 
  return(out)
}

#L1, L2, L3 dynamics

#proportion of mf per mg developing into infective larvae within the vector

delta.v <- function(delta.vo, c.v, mf, expos)
  
{
  
  out <- delta.vo / (1 + c.v * mf *expos)
  #out <- delta.vo / (1 + c.v * mf)
  
  return(out)
  
}

calc.L1 <- function(beta, mf, mf.delay.in, expos, delta.vo, c.v, nuone, mu.v, a.v, expos.delay)
  
{
  delta.vv <- delta.v(delta.vo, c.v, mf, expos)
  
  #out <- (delta.vv * beta * expos *  mf)  / (mu.v + a.v * mf + nuone) #no delay
  
  out <- (delta.vv * beta * expos *  mf)  / ((mu.v + a.v * mf*expos) + (nuone * exp (-(4/366) * (mu.v + (a.v * mf.delay.in*expos.delay)))))
  
  return(out)
}

calc.L2 <- function(nuone, L1.in, mu.v, nutwo, mf, a.v, expos) #### need to have delay expos
  
{
  #out <- (nuone * L1.in) / (mu.v + nutwo) #no delay
  
  out <- (L1.in * (nuone * exp (-(4/366) * (mu.v + (a.v * mf * expos))))) / (mu.v + nutwo)
  
  #if(out < 0) {print('negative L2'); out <- 0}
  
  return(out)
}

calc.L3 <- function(nutwo, L2.in, a.H, g, mu.v, sigma.L0)
  
{
  out <- (nutwo * L2.in) / ((a.H / g) + mu.v + sigma.L0) 
  
  #if(out < 0) {print('negative L3'); out <- 0}
  
  return(out)
}

#rate of acquisition of new infections in humans

Wplus1.rate <- function(delta.hz, delta.hinf, c.h, L3, m , beta, expos, DT)
  
{
  
  dh <- delta.h(delta.hz, delta.hinf, c.h, L3, m , beta, expos)
  
  out <- DT * m * beta * dh * expos * L3
  
  return(out)
  
}



#calculate number of mf in skin snip


mf.per.skin.snip <- function(ss.wt, num.ss, slope.kmf, int.kMf, data, nfw.start, fw.end,  ###check vectorization 
                             mf.start, mf.end, pop.size)
  
{
  
  all.mfobs <- c()
  
  kmf <- slope.kmf * (rowSums(data[,nfw.start:fw.end])) + int.kMf #rowSums(da... sums up adult worms for all individuals giving a vector of kmfs
  
  mfobs <- rnbinom(pop.size, size = kmf, mu = ss.wt * (rowSums(data[,mf.start:mf.end])))
  
  nans <- which(mfobs == 'NaN'); mfobs[nans] <- 0
  
  if(num.ss > 1)
    
  {
    
    tot.ss.mf <- matrix(, nrow = length(data[,1]), ncol = num.ss)
    tot.ss.mf[,1] <- mfobs
    
    for(j in 2 : (num.ss)) #could be vectorized
      
    {
      
      temp <- rnbinom(pop.size, size = kmf, mu = ss.wt * (rowSums(data[,mf.start:mf.end])))
      
      nans <- which(temp == 'NaN'); temp[nans] <- 0
      
      tot.ss.mf[,j] <- temp
      
    }
    
    mfobs <- rowSums(tot.ss.mf)  
    
  } 
  
  mfobs <- mfobs / (ss.wt * num.ss) 
  # mfobs <- mfobs / (num.ss) 
  
  list(mean(mfobs), mfobs)
  
}

#calculate CMFL from skin snip output

calc.CMFL <- function(ss.mf.pi, data, ss.wt)
  
{
  ind.ov.20 <-  which(data[,2] >= 20)
  mf.o.t <- log((ss.mf.pi[ind.ov.20] * ss.wt) + 1)
  CMFL <- exp(sum(mf.o.t)/length(ind.ov.20)) - 1  
  
  return(CMFL)
}


change.worm.per.ind<- function(delta.hz, delta.hinf, c.h, L3, m , beta, compartment, total.dat, num.comps,
                               w.f.l.c, lambda.zero, omeg, expos, ws, DT, mort.rates, time.each.comp, new.worms.m, new.worms.nf.fo,
                               rn, lam.m, phi, treat.stop, iteration, treat.int, treat.prob, cum.infer, treat.vec, give.treat, treat.start, N, onchosim.cov, times.of.treat)
  
{
  
  N <- length(treat.vec)
  
  lambda.zero.in <- rep(lambda.zero * DT, N) #loss of fertility
  omeg <- rep(omeg * DT, N) #becoming fertile
  
  #male worms
  
  cl <- (ws-1) + compartment #calculate which column to use depending on sex, type (fertile or infertile) and compartment
  
  cur.Wm <- total.dat[, cl] #take current number of worms from matrix
  
  worm.dead.males <- rbinom(N, cur.Wm, rep(mort.rates[compartment], N))
  worm.loss.males <- rbinom(N, (cur.Wm - worm.dead.males), rep((DT / time.each.comp), N))
  
  
  if(compartment == 1)
    
  {
    
    male.tot.worms <- cur.Wm + new.worms.m - worm.loss.males - worm.dead.males
  }
  
  if(compartment > 1)
    
  {
    male.tot.worms <- cur.Wm + w.f.l.c[[2]] - worm.loss.males - worm.dead.males
  }
  
  #female worms
  
  clnf <- (ws - 1) + num.comps + compartment #column for infertile females, num.comps skips over males
  
  clf <- (ws - 1) + 2*num.comps + compartment #column for fertile females, 2*num.comps skips over males and infertile females
  
  cur.Wm.nf <- total.dat[, clnf] #take current worm number from matrix infertile females 
  
  cur.Wm.f <- total.dat[, clf] #take current worm number from matrix fertile females
  
  mort.fems <- rep(mort.rates[compartment], N)
  
  #####################################################  
  ######### 
  #treatment
  #########
  
  #check if a complier, older than five, time for treatment, still within treatment timeframe (if standard coverage and compliance)
  #approach assumes individuals which are moved from fertile to non fertile class due to treatment re enter fertile class at standard rate 
  
  #treat.times <- seq(treat.start, treat.stop)
  
  
  if(give.treat == 1 & iteration >= treat.start) 
  {
    
    # if((((iteration-1)*DT) %% treat.int < DT) & iteration <= treat.stop) #is it a treatment time
    
    #if((((((iteration-1)*DT) %% treat.int) + 0.00001) < DT) & iteration <= treat.stop)  
    
    if((sum(times.of.treat == iteration) == 1) & iteration <= treat.stop)
      
    {
      print('TREATMENT GIVEN')
      #inds.to.treat <- which(total.dat[,1]  == 0 & total.dat[,2] > 5 & rn <= treat.prob)  #find individuals which are compliers, older than 5 and under the coverage
      
      inds.to.treat <- which(onchosim.cov == 1) 
      
      treat.vec[inds.to.treat]  <-  (iteration-1) * DT #alter time since treatment 
      if(iteration > treat.start) {mort.fems[inds.to.treat] <- mort.fems[inds.to.treat] + (cum.infer); print('cum.infer')} #alter mortality 
    }
    
    
    tao <- ((iteration-1)*DT) - treat.vec #CHECK #vector of toas, some will be NA
    
    lam.m.temp <- rep(0, N); lam.m.temp[which(is.na(treat.vec) != TRUE)] <- lam.m
    
    f.to.nf.rate <- DT * (lam.m.temp * exp(-phi * tao)) 
    
    f.to.nf.rate[which(is.na(treat.vec) == TRUE)] <- 0 #these entries in f.to.nf.rate will be NA, lambda.zero.in cannot be NA
    
    lambda.zero.in <- lambda.zero.in + f.to.nf.rate #update 'standard' fertile to non fertile rate to account for treatment 
    
  }
  ############################################################
  
  #.fi = 'from inside': worms moving from a fertile or infertile compartment
  #.fo = 'from outside': completely new adult worms 
  
  
  worm.dead.nf <- rbinom(N, cur.Wm.nf, mort.fems) #movement to next compartment
  
  worm.dead.f <- rbinom(N, cur.Wm.f, mort.fems)
  
  worm.loss.age.nf <- rbinom(N, (cur.Wm.nf - worm.dead.nf), rep((DT / time.each.comp), N))
  
  worm.loss.age.f <- rbinom(N, (cur.Wm.f - worm.dead.f), rep((DT / time.each.comp), N))
  
  # worm.loss.age.nf <- rbinom(N, (cur.Wm.nf), rep((DT / time.each.comp), N)) ####### CHANGED
  # # 
  # worm.loss.age.f <- rbinom(N, (cur.Wm.f), rep((DT / time.each.comp), N))
  
  
  #calculate worms moving between fertile and non fertile and deaths and aging 
  
  #females from fertile to infertile
  
  new.worms.nf.fi <- rep(0, N)
  
  trans.fc <- which((cur.Wm.f - worm.dead.f - worm.loss.age.f) > 0)
  #print(trans.fc)
  
  if(length(trans.fc) > 0) ##CHANGED
  {
    new.worms.nf.fi[trans.fc] <- rbinom(length(trans.fc), (cur.Wm.f[trans.fc] - worm.dead.f[trans.fc] - worm.loss.age.f[trans.fc]), lambda.zero.in[trans.fc])
  }
  
  
  #new.worms.nf.fi<- rbinom(N, cur.Wm.f, lambda.zero.in) 
  #else {new.worms.nf.fi <- 0}
  
  #females worms from infertile to fertile #this happens independent of males, but production of mf depends on males
  
  new.worms.f.fi <- rep(0, N)
  
  trans.fc <-  which((cur.Wm.nf - worm.dead.nf - worm.loss.age.nf) > 0)
  if(length(trans.fc) > 0)
  {
    new.worms.f.fi[trans.fc] <- rbinom(length(trans.fc), (cur.Wm.nf[trans.fc] - worm.dead.nf[trans.fc] - worm.loss.age.nf[trans.fc]), omeg[trans.fc])#females moving from infertile to fertile
  }
  
  
  #new.worms.f.fi<- rbinom(N, cur.Wm.nf, omeg)#females moving from infertile to fertile
  
  
  if(compartment == 1)
    
  {
    nf.out <- cur.Wm.nf + new.worms.nf.fo + new.worms.nf.fi - worm.loss.age.nf - new.worms.f.fi - worm.dead.nf#final number of infertile worms
    
    f.out <- cur.Wm.f + new.worms.f.fi - worm.loss.age.f - new.worms.nf.fi - worm.dead.f#final number of fertile worms
  }       
  
  if(compartment > 1)
    
  {
    nf.out <- cur.Wm.nf + new.worms.nf.fi - worm.loss.age.nf - new.worms.f.fi + w.f.l.c[[5]] - worm.dead.nf#w.f.l.c = worms from previous compartment
    
    f.out <- cur.Wm.f + new.worms.f.fi - worm.loss.age.f - new.worms.nf.fi + w.f.l.c[[6]] - worm.dead.f
  }   
  
  # if(male.tot.worms < 0) {male.tot.worms <- 0}
  # if(nf.out < 0) {nf.out <- 0}
  # if(f.out < 0) {f.out <- 0}
  
  list(male.tot.worms,
       worm.loss.males,
       nf.out,
       f.out,
       worm.loss.age.nf,
       worm.loss.age.f, treat.vec)  
}


weibull.mortality <- function(DT, par1, par2, age.cats)
  
{
  out <- DT  * (par1 ^ par2) * par2 * (age.cats ^ (par2-1))
  
  return(out)
}

prevalence.for.age <- function(age, ss.in, main.dat)
  
{
  
  inds <- which(main.dat[,2] >= age)
  
  out <- length(which(ss.in[[2]][inds] > 0)) / length(inds) 
  
  return(out)
}

p.flies.infected <- function(lthree) 
  
{
  kl <- (lthree * 0.1938) + 0.0017
  ana.p <- 1 - ((1 + (lthree / kl)) ^ - kl)
  return(ana.p)
}



#need to correct function inputs when running at end of code

os.age.exp <- function(full.dat)
  
{
  
  out <-  c()
  ages <- full.dat[,2]
  sexes <- full.dat[,3]
  
  ind.m <- which(sexes == 1)
  ind.f <- which(sexes == 0)
  
  #exposure for males across ages 
  out[ind.m] <- 0.05 * ages[ind.m]
  inds <- which(ages >= 20 & sexes == 1)
  if(length(inds) > 0) {out[inds] <- 1}
  
  #exposure for females across ages
  out[ind.f] <- 0.04 * ages[ind.f]
  inds <- which(ages >= 20 & sexes == 0)
  if(length(inds) > 0) {out[inds] <- 0.8}
  
  return(out)
  
}

dietz.exp <- function(age, A, beta)
{
  out <- 1 / (1 + exp(- beta * (age - A)))
  return(out)
}


ep.equi.sim <- function(delta.hz, delta.hinf, c.h,  m , beta, m.exp, f.exp, age.exp.m, age.exp.f, mu.v, int.mf, sigma.L0,
                        a.H, g, a.v, num.comps.worm, real.max.age, N, mean.age, time.its, int.L3, int.L2, int.L1, lambda.zero, omeg, delta.vo, c.v, 
                        num.mf.comps, DT, int.worms, ss.wt, num.ss, slope.kmf, int.kMf, sex.rat, nuone, nutwo, mean.in.ln, sig.in.ln,
                        time.each.comp.worms, time.each.comp.mf, mu.w1, mu.w2, fec.w.1, fec.w.2, mu.mf1, mu.mf2, it.thresh, 
                        sero.sensitivity.in, mf.move.rate, type.sero, l3.delay, dt.days,
                        input.var, lam.m, phi, treat.int, treat.int1, treat.prob, cum.infer,
                        up, kap, give.treat,treat.start, treat.stop, input.eq.mat, 
                        input.sensitivity, input.exposure, input.delay.L3, input.ser.vec, pnc, gam.dis, age.prev.in, input.eq, out.all, record.age.profs, 
                        num.age.bins.in, act.age.bins.in, threshold, smpl.sero, sens, specs, E0, q, test.thresh,
                        vector_control_start, vector_control_stop, 
                        treat.start2,  treat.start3, treat.start4)
  
  
{ 
  times.of.treat.in <- c(seq(treat.start, treat.stop - (treat.int / DT), treat.int / DT))
  #The above line defines the annnual biannual treatment
  
  stop.t <- 0
  years.req <- rep(NA, length(smpl.sero)) # if the threshold is never met, years will just be the maxium duration of treatment
  end.thresh <- vector(length = 0) #if the threshold is never met the length of this vector will be 0 and can be checked when loading
  #new.time <- 0 
  #columns to set to zero when an individual dies
  
  cols.to.zero <- seq(from = 1, to = (6 + num.mf.comps + 3*num.comps.worm))  #should this include 1 for treatment?
  cols.to.zero <- cols.to.zero[-c(1,5, 6)] #compliance, L2 and L3 do not become zero when an individual dies
  
  #used to perform operations on different worm and mf compartments 
  tot.worms <- num.comps.worm*3
  num.cols <- 6 + num.mf.comps + tot.worms 
  worms.start <- 7 + num.mf.comps
  
  
  nfw.start <- 7 + num.mf.comps + num.comps.worm #start of infertile worms
  fw.end <- num.cols #end of fertile worms 
  mf.start <- 7
  mf.end <- 6 + num.mf.comps
  
  #age dependent mortality and fecundity rates
  
  age.cats <- seq(0, 20, length = num.comps.worm)
  
  mort.rates.worms <- weibull.mortality(DT = DT, par1 = mu.w1, par2 = mu.w2, age.cats = age.cats)
  
  #mort.rates.worms <- rep(0.1*DT, num.comps.worm) #to test no senescence 
  
  fec.rates.worms <- 1.158305 * fec.w.1 / (fec.w.1 + (fec.w.2 ^ -age.cats) - 1) #no DT - Rk4
  
  #fec.rates.worms <- rep(1.15, num.comps.worm) #to test no senescence 
  
  age.cats.mf <- seq(0, 2.5, length = num.mf.comps)
  
  mort.rates.mf <- weibull.mortality(DT = 1, par1 = mu.mf1, par2 = mu.mf2, age.cats = age.cats.mf)
  
  #mort.rates.mf <- rep(0, length = num.mf.comps)  #onchosim mf mortality
  
  #mort.rates.mf[6] <- 183 #onchosim mf mortality
  
  #mort.rates.mf <- rep(1.2, num.mf.comps) #to test no senescence 
  
  
  if(input.var == 0) #are we inputting an equilibirum condition?
    
  {
    
    sero.prev.vec <- 0
    #sero sensitivity 
    
    num <- ceiling(sero.sensitivity.in * N)
    ind.sensit <- sample(seq(1, N), num)
    
    sensitivity <- rep(0, N)
    sensitivity[ind.sensit] <-  1 #infected individuals either test positive or not for life
    
    
    #create inital age distribution and simulate stable age distribution
    
    cur.age <- rep(0, N)
    
    for(i in 1 : 75000) #if at equilibrium you saved the age at which inds die and simulated further, you should get an exponential distribution
    {
      cur.age <- cur.age + DT
      
      death.vec <- rbinom(N, 1, (1/mean.age) * DT) 
      
      cur.age[which(death.vec == 1)] <- 0 #set individuals which die to age 0
      cur.age[which(cur.age >= real.max.age)] <- 0
    }
    
    
    #create list to store matrices (where required) and mean number of parasites per person
    
    all.mats <- vector("list", length = time.its + 1) 
    
    ex.vec <-rgamma(N, gam.dis, gam.dis)
    
    ###############################################
    #matrix for delay in l3 establishment in humans 
    num.delay.cols <- l3.delay * (28 / dt.days) 
    l.extras <- matrix(0, ncol= num.delay.cols, nrow= N)
    inds.l.mat <- seq(2,(length(l.extras[1,]))) #for moving columns along with time
    
    ################################################
    l1.delay <- rep(int.L1, N)
    
    ###############################################
    #matrix for tracking mf for l1 delay
    num.mfd.cols <- 4 / dt.days
    mf.delay <- matrix(int.mf, ncol= num.mfd.cols, nrow= N)
    inds.mfd.mats <- seq(2,(length(mf.delay[1,])))
    #mf.delay <- rep(int.mf, N) 
    
    ###############################################
    #matrix for exposure for L1 delay
    num.exp.cols <- 4 / dt.days
    exposure.delay <- matrix(ex.vec, ncol= num.exp.cols, nrow= N)
    inds.exp.mats <- seq(2,(length(exposure.delay[1,])))  
    
    #matrix for first time step
    all.mats[[1]] <- matrix(, nrow=N, ncol=num.cols)
    
    all.mats[[1]][,  (worms.start) : num.cols] <- int.worms
    
    all.mats[[1]][, 4] <- int.L1
    
    all.mats[[1]][, 5] <- int.L2
    
    all.mats[[1]][, 6] <- int.L3
    
    all.mats[[1]][, 7 : (7 + (num.mf.comps-1))] <- int.mf
    
    all.mats[[1]][,1] <- rep(0, N) #column used during treatment
    all.mats[[1]][,2] <- cur.age
    
    #assign sex to humans 
    
    sex <- rbinom(N, 1, sex.rat)
    
    all.mats[[1]][,3] <- sex
    
    #mf per snip and CMFL
    
    temp <- mf.per.skin.snip(ss.wt, num.ss, slope.kmf, int.kMf, data = all.mats[[1]], nfw.start, fw.end, mf.start, mf.end, pop.size = N)
    mean.mf.per.snip <- temp[[1]]
    CMFL<- calc.CMFL(ss.mf.pi = temp[[2]], data = all.mats[[1]], ss.wt)
    prev.by.ss <- prevalence.for.age(age = age.prev.in, ss.in = temp, main.dat = all.mats[[1]])
    
    sero.prev <- c()
    sero.prev.9 <- c() #0-9
    sero.prev.4 <- c() #0-4
    sero.prev.5_9 <- c()#5-9
    sero.prev.10_14 <- c() # 10-14
    sero.prev.15_19 <- c() # 15-19
    sero.prev.5_14 <- c() # 5-14
    
    #test.prev <- c()
    
    L3 <- mean(all.mats[[1]][, 6])
    
    L2 <- mean(all.mats[[1]][, 5])
    
    L1 <- mean(all.mats[[1]][, 4])
    
    non.comp <- ceiling(N * pnc)
    out.comp <- rep(0, N)
    s.comp <- sample(N, non.comp)
    out.comp[s.comp] <- 1
    all.mats[[1]][,1] <- out.comp
    
    treat.vec.in <- rep(NA, N) #for time since treatment calculations 
    
    pfert.worms.per.person <- rowSums(all.mats[[1]][,(mf.end + 1) : fw.end]) 
    
    worm.prev <- length(which(pfert.worms.per.person > 0)) / length(pfert.worms.per.person) 
    
    prop.flies.infected <- p.flies.infected(L3) 
    
  }
  
  ###
  #if inputing existing population
  ####
  
  
  
  else
    
  {
    
    all.mats <- list(); all.mats[[1]] <- input.eq[[2]]
    ex.vec <- input.eq[[3]]
    sensitivity <- input.eq[[8]]
    #sero.prev.vec <- input.eq[[12]]
    exposure.delay <- input.eq[[17]]
    
    #num.exp.cols <- 2 #THIS NEEDS TO MADE FLEXIBLE WHEN DELAY SITUATION IS 'SOLVED'
    inds.exp.mats <- seq(2,(length(exposure.delay[1,]))) 
    
    l.extras <- input.eq[[9]]
    inds.l.mat <- seq(2,(length(l.extras[1,])))
    
    l1.delay <- input.eq[[13]]
    #num.mfd.cols <- 2 #THIS NEEDS TO MADE FLEXIBLE WHEN DELAY SITUATION IS 'SOLVED'
    mf.delay <- input.eq[[14]]
    inds.mfd.mats <- seq(2,(length(mf.delay[1,])))
    
    
    L3 <- mean(all.mats[[1]][, 6])
    
    L2 <- mean(all.mats[[1]][, 5])
    
    L1 <- mean(all.mats[[1]][, 4])
    
    temp <- mf.per.skin.snip(ss.wt, num.ss, slope.kmf, int.kMf, data = all.mats[[1]], nfw.start, fw.end, mf.start, mf.end, pop.size = N)
    mean.mf.per.snip <- temp[[1]]
    CMFL<- calc.CMFL(ss.mf.pi = temp[[2]], data = all.mats[[1]], ss.wt)
    
    inds.five <- which(all.mats[[1]][,2] >= 5)
    
    prev.by.ss <- length(which(temp[[2]][inds.five] > 0)) / length(inds.five) 
    
    treat.vec.in <- rep(NA, N)
    
    if(type.sero == 1)
    {
      sero.prev.vec <- input.eq[[12]]
      sero.prev <- sum(sero.prev.vec) / N
      sero.prev.9  <- sum(sero.prev.vec[which(all.mats[[1]][,2] < 10)]) / length(which(all.mats[[1]][,2] < 10))
    }
    
    #treat.vec.in <- input.eq[[8]]#for time since treatment calculations 
    
    #if(input.var == 'treatment') {iter.last.treat <- input.eq[[10]]}
    
    
    pfert.worms.per.person <- rowSums(all.mats[[1]][,(mf.end + 1) : fw.end]) 
    
    worm.prev <- length(which(pfert.worms.per.person > 0)) / length(pfert.worms.per.person) 
    
    prop.flies.infected <- p.flies.infected(L3) 
    
  }
  
  #nw.rate <- 1 ###may not be correct when an already treated population is an input
  sero.prev.by.age.m <- list(); sero.prev.by.age.f <- list()
  prev.mf.m <- list(); prev.mf.f <- list()
  mf.intens.m <- list(); mf.intens.f <- list()
  
  coi.in <- runif(500, 0, 1)
  
  i <- 1
  
  while(i < time.its) #over time
    
  {
    print(paste(round(i * DT, digits = 2), 'yrs', sep = ' '))
    if(i > 1) 
    {
      means <- vector(length=2)
      
      L3 <- mean(all.mats[[i-1]][, 6])
      
      if(num.mf.comps == 1){means[1] <- mean(as.numeric((all.mats[[i-1]][, 7])))} 
      else means[1] <- mean(as.numeric(rowSums(all.mats[[i-1]][, 7:(worms.start-1)]))) #means for mf
      
      means[2] <- mean(as.numeric(rowSums(all.mats[[i-1]][, worms.start : ((worms.start  + num.comps.worm * 3) - 1 )]))) #means for adult worms
      
      all.mats[[i-1]] <- c(means, L3) #replace matrix at previous timestep with means for worms and L3
    }
    
    all.mats.cur <- all.mats[[i]] #create temporary matrix for t rather than using list element to increase speed
    all.mats.temp <- all.mats.cur
    
    
    if(i >= treat.start) {cov.in <- os.cov(all.dt = all.mats.cur, pncomp = pnc, covrg = C1, coi = coi.in)}
    if(i >= treat.start2) {cov.in <- os.cov(all.dt = all.mats.cur, pncomp = pnc, covrg = C2, coi = coi.in)}
    if(i >= treat.start3) {cov.in <- os.cov(all.dt = all.mats.cur, pncomp = pnc, covrg = C3, coi = coi.in)}
    if(i >= treat.start){print(sum(cov.in))}    
    
    if(i < vector_control_start) {ABR <- ABR_current}
    if(i == vector_control_start) {ABR <- ABR_current*vc.eff}
    if(i > vector_control_stop) {ABR <- ABR_current}
    
    if(i > 0)  {m <-  ABR * ((1/104) / prop.hums)}
    
    
    #if(i < vector_control_start) {m <- ABR * ((1/104) / prop.hums)}
    #if(i > vector_control_start) {m <- (ABR*0.01) * ((1/104) / prop.hums)}
    #if(i > vector_control_stop)  {m <-  (ABR) * ((1/104) / prop.hums)}
    print(m); print(ABR)
    
    
    #if(i == vector_control_start) {print(vc_eff_current)}
    #if(i > vector_control_stop) {vc_eff_current <- 1}
    
    #sex and age dependent exposure
    
    mls <- which(all.mats.cur[,3] == 1)
    fmls <- which(all.mats.cur[,3] == 0)
    
    #s.a.exp <- c()
    s.a.exp <- rep(0, N)
    
    #s.a.exp <- os.age.exp(all.mats.cur)
    #s.a.exp[mls] <- dietz.exp(age = all.mats.cur[mls, 2], A = 10, beta = 0.55)
    s.a.exp[mls] <- m.exp * exp(-age.exp.m * (all.mats.cur[mls, 2]))
    
    gam.m <- 1 / mean(s.a.exp[mls]) #normalize so mean = 1
    s.a.exp[mls] <- s.a.exp[mls] * gam.m
    
    s.a.exp[fmls] <- f.exp * exp(-age.exp.f * (all.mats.cur[fmls, 2]))
    #s.a.exp[fmls] <- dietz.exp(age = all.mats.cur[fmls, 2], A = 10, beta = 0.55)
    
    gam.f <- 1 / mean(s.a.exp[fmls]) #normalize so mean = 1
    s.a.exp[fmls] <- s.a.exp[fmls] * gam.f
    
    # print(mean(s.a.exp[fmls]))
    # print(mean(s.a.exp[mls]))
    
    ex.vec <- ex.vec * (1 / mean(ex.vec)) #normalize so mean = 1
    
    tot.ex.ai <- s.a.exp * ex.vec
    tot.ex.ai <- tot.ex.ai * (1 / mean(tot.ex.ai)) #normalize so mean = 1
    
    #increase age (for next time step)
    
    all.mats.temp[,2] <- (all.mats.cur[,2]) + DT #increase age for all individuals
    
    death.vec <- rbinom(N, 1, (1/mean.age) * DT) #select individuals to die
    
    to.die <- which(death.vec == 1)
    
    at.ab.max <- which(all.mats.temp[,2] >= real.max.age)
    
    to.die <- c(to.die, at.ab.max)
    
    to.die <- unique(to.die) #may have repeated indivudals i.e selected by binom and >80
    
    ##################
    #delay calculations 
    ##################
    
    new.worms.m <- c()
    new.worms.nf <- c()
    
    new.worms.m <- rbinom(N, size = l.extras[,length(l.extras[1,])], prob = 0.5) #draw males and females from last column of delay matrix
    new.worms.nf <- l.extras[,length(l.extras[1,])] - new.worms.m
    
    #move individuals along
    
    l.extras[,inds.l.mat] <- l.extras[,(inds.l.mat-1)]
    
    #new establishing L3 vectorized for all individuals
    
    L3.in <- mean(all.mats.cur[, 6])
    
    nw.rate <- Wplus1.rate(delta.hz, delta.hinf, c.h, L3 = L3.in, m ,
                           beta, expos = tot.ex.ai, DT)
    
    
    # if(i == treat.start) {print(L3); print(L3.in)}
    
    new.worms <- rpois(N, nw.rate) #total new establishing L3 for each individual 
    
    l.extras[,1] <- new.worms
    
    
    rann <- runif(N, 0, 1) #random number for each individual to see if treatment is given (depending on compliance)
    
    for(k in 1 : num.comps.worm) #go through each adult worm compartment
      
    {
      
      if(k == 1) {from.last <- rep(0, N)} #create vector for worms coming from previous compartment (needs to be 0 when k ==1)
      
      
      # if(i == 1) {treat.vec.in <- treat.vec} #debugging 
      
      res <- change.worm.per.ind(delta.hz = delta.hz, delta.hinf = delta.hinf, c.h = c.h, L3 = L3.in, m = m , beta = beta, compartment = k, 
                                 total.dat = all.mats.cur, num.comps = num.comps.worm,
                                 w.f.l.c = from.last, lambda.zero = lambda.zero, omeg = omeg, expos = tot.ex.ai, 
                                 ws = worms.start, DT = DT, mort.rates = mort.rates.worms, time.each.comp = time.each.comp.worms, new.worms.m = new.worms.m, 
                                 new.worms.nf.fo = new.worms.nf, rn = rann, lam.m = lam.m, phi = phi, treat.stop = treat.stop, iteration = i, treat.int = treat.int, treat.prob = treat.prob, 
                                 cum.infer = cum.infer, treat.vec = treat.vec.in, 
                                 give.treat = give.treat, treat.start = treat.start, N = N, onchosim.cov = cov.in, times.of.treat = times.of.treat.in)
      
      
      
      from.last <- res #assign output to use at next iteration, indexes 2, 5, 6 (worms moving through compartments)
      
      #update male worms in matrix for compartment k
      
      all.mats.temp[, (6 + num.mf.comps + k)] <- res[[1]]
      
      #update females worms in matrix
      
      all.mats.temp[, (6 + num.mf.comps + num.comps.worm + k)] <- res[[3]] #infertile, num.comps.worm skips over males
      all.mats.temp[, (6 + num.mf.comps + 2*num.comps.worm + k)] <- res[[4]] #fertile, num.comps.worm skips over males and infertile females
      
      
    }
    
    if(give.treat == 1 & i >= treat.start) {treat.vec.in <- res[[7]]}
    
    for(mf.c in 1 : num.mf.comps)   
      
    {
      
      res.mf <- change.micro(dat = all.mats.cur, num.comps =num.comps.worm, mf.cpt = mf.c, 
                             num.mf.comps = num.mf.comps, ws=worms.start, DT=DT, time.each.comp = time.each.comp.mf, 
                             mu.rates.mf = mort.rates.mf, fec.rates = fec.rates.worms, mf.move.rate = mf.move.rate, up = up, kap = kap, iteration = i, treat.vec = treat.vec.in, give.treat = give.treat, treat.start = treat.start)
      
      all.mats.temp[, 6 + mf.c] <- res.mf
    }
    
    
    ####l1 and mf delay
    
    exp.delay.temp <- exposure.delay[, length(exposure.delay[1,])]
    mf.delay.temp <- mf.delay[, length(mf.delay[1,])]
    l1.delay.temp <- l1.delay
    
    exposure.delay[, inds.exp.mats] <- exposure.delay[, (inds.exp.mats -1)]
    mf.delay[, inds.mfd.mats] <- mf.delay[, (inds.mfd.mats - 1)] 
    
    #new L1 L2 and L3
    
    mf.temp <- rowSums(all.mats.cur[, 7 : (6 + num.mf.comps)]) #sum mf over compartments, mf start in column 7
    
    all.mats.temp[, 4] <- calc.L1(beta, mf = mf.temp, mf.delay.in = mf.delay.temp, expos = tot.ex.ai, delta.vo, c.v, nuone, mu.v, a.v, expos.delay = exp.delay.temp)
    all.mats.temp[, 5] <- calc.L2(nuone, L1.in = l1.delay.temp, mu.v, nutwo, mf = mf.delay.temp, a.v, expos = exp.delay.temp)
    all.mats.temp[, 6] <- calc.L3(nutwo, L2.in = all.mats.cur[, 5], a.H, g, mu.v, sigma.L0)
    
    #add new parasites 
    
    l1.delay <- all.mats.temp[, 4]
    mf.delay[, 1] <- rowSums(all.mats.cur[, 7 : (6 + num.mf.comps)])
    exposure.delay[, 1] <- tot.ex.ai
    
    
    ######IVM kills l3 l4
    
    # if(i >= treat.start & (sum(times.of.treat.in == i) == 1) & i <= treat.stop) 
    #   
    #   #if((((((i-1)*DT) %% treat.int) + 0.00001) < DT) & i >= treat.start & i <= treat.stop)  
    # {
    #   #print('yes')
    #   inds.prop <- which(cov.in == 1)
    #   l.extras[inds.prop, 1:14] <- 0
    # } 
    
    #if(i >= treat.start & sum(times.of.treat.in == i) == 1) {print(l.extras[, 1:14])}
    ####
    #######
    
    
    #new individual exposure for newborns, clear rows for new borns
    
    if(length(to.die) > 0)
    {
      ex.vec[to.die] <- rgamma(length(to.die), gam.dis, gam.dis)
      
      coi.in[to.die] <- runif(length(to.die), 0, 1)
      
      l.extras[to.die, ] <- 0 #establishing adult worms 
      
      
      mf.delay[to.die, 1] <- 0 #individual dies so no contribution to l1s at this timestep
      #l1.delay[to.die,1] <- 0
      l1.delay[to.die] <- 0
      
      treat.vec.in[to.die] <- NA
      
      all.mats.temp[to.die, cols.to.zero] <- 0 #set age, sex and parasites to 0 (includes L1, but not L2 L3)
      all.mats.temp[to.die, 3] <- rbinom(length(to.die), 1, 0.5) #draw sex
    }
    
    
    all.mats[[i + 1]] <- all.mats.temp #update list of matrices
    
    #update mf per snip and CMFL
    
    temp <- mf.per.skin.snip(ss.wt, num.ss, slope.kmf, int.kMf, data = all.mats.temp, nfw.start, fw.end, mf.start, mf.end, pop.size = N)
    
    mean.mf.per.snip[i + 1] <- temp[[1]]
    
    CMFL[i + 1] <- calc.CMFL(ss.mf.pi = temp[[2]], data = all.mats.temp, ss.wt)
    
    prev.by.ss[i + 1] <- prevalence.for.age(age = age.prev.in, ss.in = temp, main.dat = all.mats.temp)
    
    
    L2[i + 1] <- mean(all.mats.temp[, 5])
    
    L1[i + 1] <- mean(all.mats.temp[, 4])
    
    
    pfert.worms.per.person <- rowSums(all.mats.temp[,(mf.end + 1) : fw.end]) 
    
    worm.prev[i + 1] <- length(which(pfert.worms.per.person > 0)) / length(pfert.worms.per.person) 
    
    prop.flies.infected[i + 1]  <- p.flies.infected(mean(all.mats.temp[,6])) 
    
    
    #if sero prev =< thresh , treat.stop <- timestep - 1 (controlled for DT)
    #if treat.stop >= max.treat.years, some.var <- 1, else some.var <- 0 (add some.var to return)
    
    ##age structured infection outputs
    
    ##################################################
    ###SEROLOGY & MF AGE INTENSITY & INTENSITY PROFILES
    
    mf.temp <- rowSums(all.mats.temp[, 7 : (6 + num.mf.comps)])
    mls <- which(all.mats.temp[,3] == 1)
    fmls <- which(all.mats.temp[,3] == 0)
    
    if(i == start.sero){
      if(input.var == 0){
        sero.prev.vec.main <- rep(0, N)
        sero.prev.vec <- rep(0, N)
        
        ind.one <- which(rowSums(all.mats.temp[, worms.start:num.cols]) > 0 & sensitivity == 1) # this is only updated when a new individual is infected, or an infected individuals dies
        
        
        
        sero.prev.vec.main[ind.one] <- 1 #who has at least one adult worm and can mount a response 
        
        inds.infected <- which(sero.prev.vec.main > 0)
        inds.not.infected <- which(sero.prev.vec.main == 0)
        
        temp.sen <- rep(0, N)
        temp.spec <- rep(0, N)
        
        gen.vec <- seq(1, N)
        
        inds.sen.inf <- sample(inds.infected, length(inds.infected) * sens)
        sen.inf <- inds.sen.inf
        
        inds.spec.nf<- sample(inds.not.infected, length(inds.not.infected) * (1 - specs))
        
        spec.nf <- inds.spec.nf
        
        temp.sen[sen.inf] <-1 #0s here indicate false and true negatives
        temp.spec[spec.nf] <- 1 #0s here indicate false and true positives
        
        ind.two<- which(temp.sen == 1) 
        
        ind.three<- which(temp.spec == 1) 
        
        
        sero.prev.vec[unique(c(ind.two, ind.three))] <- 1
      }
      
      else{ ####THIS NEEDS TO BE UPDATED
        ind <- which(rowSums(all.mats.temp[, worms.start:num.cols]) > 0 & sero.prev.vec == 0 & sensitivity == 1) #which inds have worms but didn't at the last time step
        sero.prev.vec[ind] <- 1 #new (adult worm) infections become seropositive
        if(length(to.die) > 0) {sero.prev.vec[to.die] <- 0} #newborns are negative
      }
    }
    
    if(i > start.sero){
      
      ####
      #######assumes individuals not able to mount are response are not infected 
      
      ind <- which(rowSums(all.mats.temp[, worms.start:num.cols]) > 0 & sero.prev.vec.main == 0 & sensitivity == 1) #which inds have worms but didn't at the last time step
      sero.prev.vec.main[ind] <- 1 #new (adult worm) infections become seropositive
      
      sero.prev.vec.main[to.die] <- 0 #newborns are negative
      
      sero.prev.vec <- rep(0, N)
      #assume tests are conducted at each time step, so resample individuals mounting a response based of sensitivity and sensitivity
      
      inds.infected <- which(sero.prev.vec.main > 0)
      inds.not.infected <- which(sero.prev.vec.main == 0)
      
      temp.sen <- rep(0, N)
      temp.spec <- rep(0, N)
      
      gen.vec <- seq(1, N)
      
      inds.sen.inf <- sample(inds.infected, length(inds.infected) * sens)
      sen.inf <- inds.sen.inf
      
      inds.spec.nf<- sample(inds.not.infected, length(inds.not.infected) * (1 - specs))
      
      spec.nf <- inds.spec.nf
      
      temp.sen[sen.inf] <-1 #0s here indicate false and true negatives
      temp.spec[spec.nf] <- 1 #0s here indicate false and true positives
      
      ind.two<- which(temp.sen == 1) 
      
      ind.three<- which(temp.spec == 1) 
      # 
      # print(length(ind.two))
      # print(length(ind.three))
      
      sero.prev.vec[unique(c(ind.two, ind.three))] <- 1
      
      #####
      #######
      # 
      # print(sum(temp.sen) / N)
      # print(1 - sum(temp.spec) / N)
      
    }
    
    if(i >= start.sero){
      sero.prev[i + 1] <- sum(sero.prev.vec) / N #prevalence in whole population 
      
      #print(sum(sero.prev.vec))
      
      sero.prev.9 [i + 1] <- sum(sero.prev.vec[which(all.mats.temp[,2] < 10)]) / length(which(all.mats.temp[,2] < 10)) #prevalence in under 10s
      
      
      sero.prev.4[i + 1] <- sum(sero.prev.vec[which(all.mats.temp[,2] < 5)]) / length(which(all.mats.temp[,2] < 5)) #0-4
      sero.prev.5_9[i + 1] <- sum(sero.prev.vec[which(all.mats.temp[,2] < 10 & all.mats.temp[,2] >= 5)]) / length(which(all.mats.temp[,2] < 10 & all.mats.temp[,2] >= 5))#5-9
      sero.prev.10_14[i + 1] <- sum(sero.prev.vec[which(all.mats.temp[,2] < 15 & all.mats.temp[,2] >=10)]) / length(which(all.mats.temp[,2] < 15 & all.mats.temp[,2] >=10)) # 10-14
      sero.prev.15_19[i + 1] <- sum(sero.prev.vec[which(all.mats.temp[,2] < 20 & all.mats.temp[,2] >= 15)]) / length(which(all.mats.temp[,2] < 20 & all.mats.temp[,2] >= 15)) # 15-19
      sero.prev.5_14[i + 1] <- sum(sero.prev.vec[which(all.mats.temp[,2] < 15 & all.mats.temp[,2] >= 5)]) / length(which(all.mats.temp[,2] < 15 & all.mats.temp[,2] >= 5)) # 5-14
      
      
    }
    
    if(record.age.profs == 'all') {  
      if(num.age.bins.in != 0) {bz.m <- bin(all.mats.temp[mls,2], nbins = num.age.bins.in); bz.f <- bin(all.mats.temp[fmls,2], nbins = num.age.bins.in)}
      if(length(act.age.bins.in) != 1) {bz.m <- cut(all.mats.temp[mls,2], act.age.bins.in, include.lowest = TRUE); bz.f <- cut(all.mats.temp[fmls,2], act.age.bins.in, include.lowest = TRUE)}
      
      #######NEED TO CREATE THESE LISTS AT START OF FUNCTION  
      sero.prev.by.age.m[[i]] <- tapply(sero.prev.vec[mls], bz.m, function(x) sum(x) / length(x))
      sero.prev.by.age.f[[i]] <- tapply(sero.prev.vec[fmls], bz.f, function(x) sum(x) / length(x))
      
      infected <- rep(0, N)
      
      inf <-  which(temp[[2]] > 0)
      
      infected[inf] <- 1
      
      prev.mf.m[[i]] <- tapply(infected[mls], bz.m, function(x) sum(x) / length(x))
      prev.mf.f[[i]] <- tapply(infected[fmls], bz.f, function(x) sum(x) / length(x))
      
      mf.intens.m[[i]] <- tapply(temp[[2]][mls], bz.m, mean)
      mf.intens.f[[i]] <- tapply(temp[[2]][fmls], bz.f, mean)
    }
    
    if(record.age.profs == 'last' & i == (treat.start - 1)) {
      
      if(num.age.bins.in != 0) {bz.m <- bin(all.mats.temp[mls,2], nbins = num.age.bins.in); bz.f <- bin(all.mats.temp[fmls,2], nbins = num.age.bins.in)}
      
      if(length(act.age.bins.in) != 1) {bz.m <- cut(all.mats.temp[mls,2], act.age.bins.in, include.lowest = TRUE); bz.f <- cut(all.mats.temp[fmls,2], act.age.bins.in, include.lowest = TRUE)}
      
      sero.prev.by.age.m <- tapply(sero.prev.vec[mls], bz.m, function(x) sum(x) / length(x))
      sero.prev.by.age.f <- tapply(sero.prev.vec[fmls], bz.f, function(x) sum(x) / length(x))
      
      infected <- rep(0, N)
      
      inf <-  which(temp[[2]] > 0)
      
      infected[inf] <- 1
      
      prev.mf.m <- tapply(infected[mls], bz.m, function(x) sum(x) / length(x))
      prev.mf.f <- tapply(infected[fmls], bz.f, function(x) sum(x) / length(x))
      
      mf.intens.m <- tapply(temp[[2]][mls], bz.m, mean)
      mf.intens.f <- tapply(temp[[2]][fmls], bz.f, mean)
    }
    
    if(record.age.profs == 'last' & give.treat == 1 & i == (treat.stop + (treat.int / DT)))
      
    {
      
      if(num.age.bins.in != 0) {bz.m <- bin(all.mats.temp[mls,2], nbins = num.age.bins.in); bz.f <- bin(all.mats.temp[fmls,2], nbins = num.age.bins.in)}
      
      if(length(act.age.bins.in) != 1) {bz.m <- cut(all.mats.temp[mls,2], act.age.bins.in, include.lowest = TRUE); bz.f <- cut(all.mats.temp[fmls,2], act.age.bins.in, include.lowest = TRUE)}
      
      sero.prev.by.age.m.treat <- tapply(sero.prev.vec[mls], bz.m, function(x) sum(x) / length(x))
      sero.prev.by.age.f.treat <- tapply(sero.prev.vec[fmls], bz.f, function(x) sum(x) / length(x))
      
      infected <- rep(0, N)
      
      inf <-  which(temp[[2]] > 0)
      
      infected[inf] <- 1
      
      prev.mf.m.treat <- tapply(infected[mls], bz.m, function(x) sum(x) / length(x))
      prev.mf.f.treat <- tapply(infected[fmls], bz.f, function(x) sum(x) / length(x))
      
      mf.intens.m.treat <- tapply(temp[[2]][mls], bz.m, mean)
      mf.intens.f.treat <- tapply(temp[[2]][fmls], bz.f, mean)
      
    }  
    
    
    #sero.prev.5_14[i + 1] <- sum(sero.prev.vec[which(all.mats.temp[,2] < 15 & all.mats.temp[,2] >= 5)]) / length(which(all.mats.temp[,2] < 15 & all.mats.temp[,2] >= 5)) # 5-14
    
    # if(i > treat.start & (((i-1)*DT) %% treat.int < DT) & i <= treat.stop)
    # {
    #   for(l in 1 : length(smpl.sero))  #this could be vectorized
    #     
    #   {
    #     
    #     if(sum(is.na(years.req[l])) == 0) next 
    #     all.inds <- which(all.mats.temp[,2] < 9)
    #     # print(length(all.inds))
    #     # print(length(all.inds) / smpl.sero[l])
    #     inds.sample.sero <- sample(all.inds, length(all.inds) * smpl.sero[l]) 
    #     
    #     #if element in yrs and threshold is not NA next
    #     #if length = 0 is NA break
    #     
    #     #calculate prevalence 
    #     
    #     #note inds.sample.sero might not necessarily be <9
    #     prev.9 <- sum(sero.prev.vec[inds.sample.sero]) / length(inds.sample.sero) #prevalence in under 10s
    #     
    #     
    #     if(prev.9 <= threshold & is.na(years.req[l]) == TRUE) 
    #       
    #     {
    #       
    #       years.req[l] <-  (i - treat.start) * DT#save treatment year, 
    #       # treat.stop <- i #to stop treatment  
    #       end.thresh[l] <- prev.9 #to check it has stopped at the right time 
    #       # time.its <- i + 5/DT#alter total timesteps so it just runs for 50 year mores to check elimination 
    #       #when load need to check if treatment year is less
    #       #new.time <- 1
    #     }
    #     
    #     
    #     
    #   }
    # }
    
    # print(stop.t)
    # print(sum(is.na(years.req)))
    
    #use while loop instead so i can be altered
    
    #if(sum(is.na(years.req)) == 0 & stop.t == 0 & test.thresh == 1) {time.its <- i + 50/DT; treat.stop <- i; stop.t <-  1; print('yes'); stop("checking thresh")} 
    #END SECTION
    
    i <- i + 1
    
    #if(new.time == 1 & i == time.its) break     
    
  }
  
  i <- i - 1
  
  final.mat <- all.mats[[i + 1]] #save matrix from final time step before replacing with means
  
  #same process as at beginning of main for loop
  
  means <- vector(length=2)
  
  L3 <- sum(as.numeric(all.mats[[i]][, 6])) / N
  
  if(num.mf.comps == 1){means[1] <- mean(as.numeric((all.mats[[i]][, 7])))} 
  else means[1] <- mean(as.numeric(rowSums(all.mats[[i]][, 7:(worms.start-1)]))) #sum up mf over compartments for each individual
  
  means[2] <- mean(as.numeric(rowSums(all.mats[[i]][, worms.start : ((worms.start  + num.comps.worm * 3) - 1 )])))
  
  all.mats[[i]] <- c(means, L3)
  
  means <- vector(length=2)
  
  L3 <- sum(as.numeric(all.mats[[i + 1]][, 6])) / N
  
  if(num.mf.comps == 1){means[1] <- mean(as.numeric((all.mats[[i + 1]][, 7])))}
  else means[1] <- mean(as.numeric(rowSums(all.mats[[i + 1]][, 7:(worms.start-1)])))
  
  means[2] <- mean(as.numeric(rowSums(all.mats[[i + 1]][, worms.start : ((worms.start  + num.comps.worm * 3) - 1 )])))
  
  all.mats[[i + 1]] <- c(means, L3)
  
  
  #ouput
  
  
  if(out.all == 1 & type.sero != 1) #enough to put output into a new simulation
    
  {
    return(list(all.mats, final.mat, ex.vec, mean.mf.per.snip, CMFL, prev.by.ss, temp[[2]], sensitivity, l.extras, L2, L1, sero.prev.vec, l1.delay, mf.delay, treat.vec.in, i, exposure.delay, worm.prev, prop.flies.infected))
  }  
  
  if(out.all == 1 & type.sero == 1) #enough to put output into a new simulation and plot age profiles for mf and seroprevalence
    
  {
    # return(list(all.mats, final.mat, ex.vec, mean.mf.per.snip, CMFL, prev.by.ss, temp[[2]], sensitivity, l.extras, L2, L1, sero.prev.vec, l1.delay, mf.delay, treat.vec.in, i, exposure.delay,
    #             sero.prev.by.age.m, sero.prev.by.age.f, prev.mf.m, prev.mf.f, sero.prev, sero.prev.9 , tot.ex.ai, worm.prev, prop.flies.infected,
    #             sero.prev.4, sero.prev.5_9, sero.prev.10_14, sero.prev.15_19, sero.prev.5_14, mf.intens.m, mf.intens.f,
    #             sero.prev.by.age.m.treat, sero.prev.by.age.f.treat, prev.mf.m.treat, prev.mf.f.treat, mf.intens.m.treat, mf.intens.f.treat))
    
    return(list(all.mats, final.mat, ex.vec, mean.mf.per.snip, CMFL, prev.by.ss, temp[[2]], sensitivity, l.extras, L2, L1, sero.prev.vec, l1.delay, mf.delay, treat.vec.in, i, exposure.delay,
                sero.prev.by.age.m, sero.prev.by.age.f, prev.mf.m, prev.mf.f, sero.prev, sero.prev.9 , tot.ex.ai, worm.prev, prop.flies.infected,
                sero.prev.4, sero.prev.5_9, sero.prev.10_14, sero.prev.15_19, sero.prev.5_14, mf.intens.m, mf.intens.f, end.thresh, years.req))
  }  
  
  if(out.all == 0) #basic output
    
  {
    return(list(all.mats, final.mat, ex.vec, mean.mf.per.snip, CMFL, prev.by.ss, temp[[2]], worm.prev, prop.flies.infected))
  }
  
}

#fix functionality for return age profiles after treatment 

#OEPA ELISA: sen: 45 - 60, spec: 99.98 (determination of sensitivity unclear)
#CDC / OEPA ELISE: sen: 88.2, spec: 99.7 (one study)
#PATH ELISA:sen 84.3, spec: 94 (validation study) 
#SD bioline RDT:sen 89.1 / 91, spec 97 / 98 (two studies)

#RDT POC is instant result
#RDT in lab (different from ELISA)

#must mount IGG4 to OV16

spec.in <- 1
sens.in <- 1

DT.in <- 1/366
dt.days.in <- DT.in*366

treat.len <- tlength

start.sero = 1 #when to start recording sero prev
pop.size.in <- 500


gamin <- 0.3

reuse <- 1

give.treat.in = 1

delta.hz.in <- 0.1864987
delta.hinf.in <- 0.002772749
ch.in <-  0.004900419

# ch.in <-  0
# delta.hz.in <- 0.0031

prop.hums <- 0.63

threshold.in <- -1

#input parms for exposure 

#cam
q.in <- 0
E0.in <- 0
m.exp.in <- 1.08
f.exp.in <-  0.9
age.exp.m.in <-  0.007
age.exp.f.in <-  -0.023

#### treat.start is times for ivermectin
treat.start =  round(151 / (DT.in )) #1991-1995
treat.start2 = round(156 / (DT.in )) #1996-2001
treat.start3 = round(162 / (DT.in )) #2002 onwards
treat.start4 = round(174 / (DT.in )) #2014 onwards
treat.stop = treat.start3 + round(28 / (DT.in )) # to end in 2030 (as we do not have information after the pandemic started)
##Regarding the above line (treat.stop) - Time for least treatment relative to the final coverage shift


#### vector_control_start is times for vc to start (set it relative to IVM starting)
vector_control_start = round(149 / (DT.in )); vector_control_stop = round(162 / (DT.in ))
#Vector control from 1989 to 2002

#### total simulation length
timesteps = round(240 / (DT.in )) #Until 2080

## Coverages associted with Coverages 1 --> Coverages 2 --> Coverages 3 
C1 <- 0.65
C2 <-  0.75
C3 <- 0.8

## Sytematic Non-Comp, Linked to last coverage value (larger C3, smaller SNC)
SNC <- 0.01 #1% if final coverage (C3) is 80%, 2.5% if 75%, and 5% if 65%

## 1- porpotional reduction in biting rate --> 99% reunction = 0.01
vc.eff <- 0.1  # 0.1 = 90% vector control efficiency; 0.25 = 75%; 0.4 = 60%

#treat.int.in <- 1
#ABR_current <- 600
#ABR <- ABR_current

treat.int.in <-tot.parms[iter, 2]
treat.int.in1 <- 0.5
ABR_current <- tot.parms[iter, 1]
ABR <- ABR_current 


out.tot<- ep.equi.sim(delta.hz  = delta.hz.in, delta.hinf = delta.hinf.in , c.h = ch.in,  m = ABR * ((1/104) / prop.hums) , beta = prop.hums / (1/104),
                      m.exp = m.exp.in, f.exp = f.exp.in, age.exp.m = age.exp.m.in, age.exp.f = age.exp.f.in, mu.v = 26, int.mf = 0, sigma.L0 = 52,
                      a.H = 0.8, g = 0.0096, a.v = 0.39, num.comps.worm = 21, real.max.age = 80, N = pop.size.in, mean.age = 50, time.its = timesteps,
                      int.L3 = 0.03, int.L2 = 0.03, int.L1 = 0.03, lambda.zero = 0.33, omeg = 0.59,
                      delta.vo = 0.0166, c.v = 0.0205, num.mf.comps = 21, DT = DT.in, int.worms=1, ss.wt = 2, num.ss = 2,
                      slope.kmf = 0.0478, int.kMf = 0.313, sex.rat = 0.5, nuone = 201.6189, nutwo = 207.7384,
                      time.each.comp.worms = 1, time.each.comp.mf = 0.125, mu.w1 = 0.09953, mu.w2 = 6.00569, fec.w.1 = 70, fec.w.2 = 0.72,
                      mu.mf1 = 1.089, mu.mf2 = 1.428, it.thresh=start.sero, sero.sensitivity.in = 0.8, mf.move.rate = 8.133333,
                      type.sero = 1, l3.delay = 10, dt.days = dt.days.in, input.var = 0,  treat.int = treat.int.in, treat.int1 = treat.int.in1, treat.prob = 0.5, cum.infer= 0.345,
                      lam.m = 32.4, phi = 19.6, up = 0.0096, kap = 1.25, give.treat = give.treat.in, treat.start = treat.start, treat.stop = treat.stop, pnc = SNC,  gam.dis = gamin, age.prev.in = 5, input.eq = 0, out.all = 1,
                      record.age.profs = 'last', num.age.bins.in = 0, act.age.bins.in = c(seq(0, 30, 0.5), seq(35, 80, 5)), threshold = threshold.in,
                      smpl.sero = c(1, 0.9, 0.8, 0.7, 0.6, 0.5), sens = sens.in, specs = spec.in, E0 = E0.in, q = q.in, test.thresh = 0,
                      vector_control_start = vector_control_start, vector_control_stop = vector_control_stop , 
                      treat.start2 = treat.start2,  treat.start3 = treat.start3, treat.start4=treat.start4)




  save.image(paste("outputs/",iter,".RData", sep=""))