
# ===== report_est =======
report_est <- function(Y, A, D, muhat, e, a, opts = list()){
  # horizon
  TT = length(Y)
  # delay probability
  if (is.null(opts$pD)){
    pD = rep(1, TT) # do not adjust for delay probability
  }
  else{
    pD = opts$pD
  }
  # computing adaptive weights
  h = sqrt(e) #* sqrt(pD)
  
  # report estimators
  A_flag = (A == a)
  D_flag = (D <= (TT-1:TT))
  ## hajek
  gamma = A_flag * D_flag/(e * pD)
  est_aipw = sum(h * (Y - muhat) * gamma)/sum(h * gamma) +
    sum(h * muhat) / sum(h)
  var_aipw = sum(h^2 * (Y - est_aipw)^2 * gamma^2) / sum(h * gamma)^2
  
  ## ignoring arm dependency
  est_aipw_iad = sum(h * (Y - muhat) * gamma)/sum(h * D_flag/pD) +
    sum(h * muhat) / sum(h)
  phat = sum(h * D_flag/pD)/sum(h)
  var_aipw_iad = sum(h^2 * ((Y - muhat) * gamma/phat + muhat - est_aipw_iad)^2) / sum(h)^2
  
  ## non-hajek
  est_aipw_nh = sum(h * (Y - muhat) * gamma)/sum(h) +
    sum(h * muhat) / sum(h)
  var_aipw_nh = sum(h^2 * ((Y - muhat) * gamma + muhat - est_aipw_nh)^2) / sum(h)^2
  
  ## no adaptive weighting
  est_aipw_naw = sum(1 * (Y - muhat) * gamma)/sum(1 * gamma) +
    sum(1 * muhat) / TT
  var_aipw_naw = sum(1^2 * (Y - est_aipw)^2 * gamma^2) / sum(1 * gamma)^2
  
  est = c(est_aipw, est_aipw_iad, est_aipw_nh, est_aipw_naw)
  var_est = c(var_aipw, var_aipw_iad, var_aipw_nh, var_aipw_naw)
  
  data.frame(est = est, 
             var_est = var_est,
             row.names = c("Hajek", "Unadjusted", "NoHajek", "NoWeight"))
}


# ====== get delay ======
get_delay <- function(delay_dist, pmt, partial = 1.0){
  censored <- rbinom(1, 1, 1-partial)
  D <- switch(
    delay_dist,
    NODELAY = 0,
    POISSON = rpois(1, pmt),
    NEGBINOM = rnbinom(1, pmt$size, pmt$prob),
    PARETO  = floor(rpareto(1, 1, pmt) - 1),
    0
  ) 
  
  ifelse(censored, Inf, D)
}

# ====== eps greedy ======
eps_greedy <- function(TT, d, mu, alpha, opts = list()){
  Y = rep(NA, TT)
  A = rep(NA, TT)
  D = rep(0,  TT)
  e = matrix(NA, nrow = d, ncol = TT)
  n = rep(0, d)
  avg = rep(0, d)
  AVG = matrix(0, nrow = d, ncol = TT)
  bufferD = c()
  bufferA = c()
  bufferY = c()
  
  ### extract delay setup
  if (is.null(opts$delay_dist)){
    delay_dist = rep("NODELAY", d)
    pmt = rep(NA, d)
    partial = rep(1.0, d)
  }else
  {
    delay_dist = opts$delay_dist
    pmt = opts$pmt
    partial = opts$partial
  }
  
  ### burn-in stage: randomized pull
  burn_in = 0.1*TT
  for (t in 1:burn_in){
    #### Use Buffer to update n and avg
    if (length(bufferD) != 0){
      bufferD = bufferD - 1
      uniqueA = unique(bufferA[bufferD == 0])
      
      if (length(uniqueA) != 0){
        for (a in uniqueA){
          observed = (bufferD == 0 & bufferA == a)
          avg[a] = (avg[a] * (n[a]) + sum(bufferY[observed])) / 
            (n[a] + sum(observed))
          n[a] = n[a] + sum(observed)
        }
      }
      bufferA = bufferA[bufferD != 0]
      bufferY = bufferY[bufferD != 0]
      bufferD = bufferD[bufferD != 0]
    }
    
    #### update n and avg and Buffer using data from stage t
    
    e[, t] = rep(1/d, d)
    A[t] = which(rmultinom(1, 1, e[, t]) == 1)
    Y[t] = rnorm(1, mu[A[t]])
    D[t] = get_delay(delay_dist[A[t]], pmt[[A[t]]], partial[A[t]])
    if (D[t] == 0){
      n[A[t]] = n[A[t]] + 1
      avg[A[t]] = (avg[A[t]] * (n[A[t]]-1) + Y[t]) / n[A[t]]
    }else{
      bufferD = c(bufferD, D[t])
      bufferA = c(bufferA, A[t])
      bufferY = c(bufferY, Y[t])
    }
    
    #### record AVG
    AVG[, t] = avg
  }
  
  ### greedy stage
  for (t in ((burn_in + 1):TT)){
    #### Use Buffer to update n and avg
    if (length(bufferD) != 0){
      bufferD = bufferD - 1
      uniqueA = unique(bufferA[bufferD == 0])
      if (length(uniqueA) != 0){
        for (a in uniqueA){
          observed = (bufferD == 0 & bufferA == a)
          avg[a] = (avg[a] * (n[a]) + sum(bufferY[observed])) / 
            (n[a] + sum(observed))
          n[a] = n[a] + sum(observed)
        }
      }
      bufferA = bufferA[bufferD != 0]
      bufferY = bufferY[bufferD != 0]
      bufferD = bufferD[bufferD != 0]
    }
    
    #### update n and avg and Buffer using data from stage t
    
    e[, t] = t^(-alpha)/(d-1)
    amax = which.max(avg)
    e[amax, t] = 1 - t^(-alpha)
    
    A[t] = which(rmultinom(1, 1, e[, t]) == 1)
    Y[t] = rnorm(1, mu[A[t]])
    D[t] = get_delay(delay_dist[A[t]], pmt[[A[t]]], partial[[A[t]]])
    
    if (D[t] == 0){
      n[A[t]] = n[A[t]] + 1
      avg[A[t]] = (avg[A[t]] * (n[A[t]]-1) + Y[t]) / n[A[t]]
    }else{
      bufferD = c(bufferD, D[t])
      bufferA = c(bufferA, A[t])
      bufferY = c(bufferY, Y[t])
    }
    
    #### record AVG
    AVG[, t] = avg
  }
  
  list(Y = Y, A = A, D = D, e = e, AVG = AVG)
}

