
##-----------------------
## estimate theta(j,k,p)
##-----------------------
estima <- function(data, library) {

    n <- nrow(data) # sample size
    K = max(data$study) # max study index
    W = ncol(data) - 5 # max covariate index
    

    # Predict outcome from all other variables
    fita <- SuperLearner(Y = data$y,
                         X = select(data, a, m, paste0('w', 0:W), study),
                         SL.library = library,
                         verbose = FALSE,
                         method = "method.NNLS",
                         family = if (all(data$y %in% c(0, 1))) binomial() else gaussian())
    
    
    # Predict outcome under other studies, keeping m, x, l the same
    preda <- lapply(0:K, function(ss) {
        predict(fita,
                newdata = data %>% select(a, m, paste0('w', 0:W), study) %>%
                    mutate(study = ss))$pred[,1]
    })
    
    # Predict preda for all studies given treatment and covariates
    fitb <- lapply(0:K, function(ss) {
        SuperLearner(Y = preda[[1 + ss]],
                     X = select(data, a, paste0('w', 0:W), study),
                     SL.library = library,
                     verbose = FALSE,
                     method = "method.NNLS")
    })

    # Predict treatment from covariates, mediators, and study
    fitg <-  SuperLearner(Y = pull(data, a),
                          X = select(data, m, paste0('w', 0:W), study),
                          SL.library = library,
                          verbose = FALSE,
                          method = "method.NNLS",
                          family = binomial())

    # Predict study from covariates and mediators
    fite <- lapply(0:K, function(ss) {
        SuperLearner(Y = (pull(data, study) == ss) * 1,
                     X = select(data, m, paste0('w', 0:W)),
                     SL.library = library,
                     verbose = FALSE,
                     method = "method.NNLS",
                     family = binomial())
    })
    
    # Predict study from covariates
    fitr <- lapply(0:K, function(ss) {
        SuperLearner(Y = (pull(data, study) == ss) * 1,
                     X = select(data, paste0('w', 0:W)),
                     SL.library = library,
                     verbose = FALSE,
                     method = "method.NNLS",
                     family = binomial())
    })

    # Predict treatment from study and covariates
    t1 <- SuperLearner(Y = pull(data, a),
                       X = select(data, paste0('w', 0:W), study),
                       SL.library = library,
                       verbose = FALSE,
                       method = "method.NNLS",
                       family = binomial())
    
    # Get proportion for each study
    h <- function(ss)mean(data$study == ss)

    onestep <- function(params) {

        k  <- params[1] # s_y
        p  <- params[2] # s_m
        j  <- params[3] # s_w

        b <- lapply(0:1, function(as){
            predict(fitb[[1 + k]],
                    newdata = data %>% select(a, paste0('w', 0:W), study) %>%
                        mutate(a = as, study = p))$pred[,1]
        })

        ba <- data$a * b[[2]] + (1 - data$a) * b[[1]]
        a <- preda[[1 + k]]

        diffb <- b[[2]] - b[[1]]
        theta <- mean(diffb[data$study == j])

        gp <- predict(fitg, newdata = data %>% select(m, paste0('w', 0:W), study) %>%
                                mutate(study = p))$pred[,1]
        gk <- predict(fitg, newdata = data %>% select(m, paste0('w', 0:W), study) %>%
                                mutate(study = k))$pred[,1]

        grat <- (data$a * gp + (1 - data$a) * (1 - gp)) / (data$a * gk + (1 - data$a) * (1 - gk))

        ep <- predict(fite[[1 + p]], newdata = data %>% select(m, paste0('w', 0:W)))$pred[,1]
        ek <- predict(fite[[1 + k]], newdata = data %>% select(m, paste0('w', 0:W)))$pred[,1]

        erat <- ep / ek

        rj <- predict(fitr[[1 + j]], newdata = data %>% select(paste0('w', 0:W)))$pred[,1]
        rp <- predict(fitr[[1 + p]], newdata = data %>% select(paste0('w', 0:W)))$pred[,1]

        rrat <- rj / rp

        tp  <-  predict(t1, newdata = data %>% select(paste0('w', 0:W), study) %>%
                             mutate(study = p))$pred[,1]

        tasp <- data$a * tp + (1 - data$a) * (1 - tp)

        akrat <- (2 * data$a - 1) * (data$study == k) / (tasp * h(j))
        aprat <- (2 * data$a - 1) * (data$study == p) / (tasp * h(j))

        jrat <- (data$study == j) / h(j)

        eif1 <- grat * erat * rrat * akrat * (data$y - a)
        eif2 <- rrat * aprat * (a - ba)
        eif3 <- jrat * (diffb - theta)
        estima  <-  theta + mean(eif1 + eif2 + eif3)
        eif <- eif1 + eif2 + eif3 - estima
        return(list(estima = estima, eif = eif))
    }

    vals <- expand.grid(k = 0:K, p = 0:K, j = 0:K)
    thetas <- apply(vals, 1, onestep)
    estimates <- sapply(thetas, function(x)x[[1]])
    eifs <- sapply(thetas, function(x)x[[2]])

    return(list(vals = vals, n = n, estimates = estimates, eifs = eifs))

}

# Calculate variance of vector
myvar <- function(x)mean((x - mean(x))^2)

# Calculate variance parameter
varest <- function(ests, eifs){
    kk <- length(ests)
    cons <- (kk - 1) / kk * (1 / (choose(kk, 2) * 2))
    var <- eif <- 0
    for(i in 1:kk) {
        for(j in 1:kk) {
            var <- var + cons * (j < i) * (ests[i] - ests[j])^2
            eif <- eif + 2 * cons * (j < i) * (ests[i] - ests[j]) * (eifs[, i] - eifs[,j])
        }
    }
    return(list(var = var, ests = ests, eif = eif, eifs = eifs))
}

# Extract entries for k, p
funcm <- function(valskp, out) {

    k <- valskp[1]
    p <- valskp[2]
    idx <- which(out$vals$k == k & out$vals$p == p)
    ests <- out$estimates[idx]
    eifs <- out$eifs[, idx]

    return(varest(ests, eifs))

}

# Extract entries for k
funem <- function(k, out) {

    idx <- which(out$vals$k == k)
    ests <- out$estimates[idx]
    eifs <- out$eifs[, idx]

    return(varest(ests, eifs))

}

# Estimate variance decomposition
fun_decomp <- function(out) {

    n <- out$n
    K <- length(out$vals)^(1/3)

    valskp <- expand.grid(k = 0:K, p = 0:K)
    tmp1 <- apply(valskp, 1, funcm, out)
    tcm <- mean(sapply(tmp1, function(x)x$var))
    eifcm <- rowMeans(sapply(tmp1, function(x)x$eif))
    setcm <- sqrt(myvar(eifcm) / n)

    kappa <- sapply(tmp1, function(x)mean(x$ests))
    kappaeif <- sapply(tmp1, function(x)rowMeans(x$eifs))
    tmp2 <- varest(kappa, kappaeif)
    teh <- tmp2$var
    eifeh <- tmp2$eif
    seteh <- sqrt(myvar(eifeh) / n)

    outem <- list(vals = valskp, estimates = kappa, eifs = kappaeif)
    tmp3 <- lapply(0:K, funem, outem)
    tmv <- mean(sapply(tmp3, function(x)x$var))
    eifmv <- rowMeans(sapply(tmp3, function(x)x$eif))
    setmv <- sqrt(myvar(eifmv) / n)

    gamma <- sapply(tmp3, function(x)mean(x$ests))
    gammaeif <- sapply(tmp3, function(x)rowMeans(x$eifs))
    tmp4 <- varest(gamma, gammaeif)
    tem <- tmp4$var
    eifem <- tmp4$eif
    setem <- sqrt(myvar(eifem) / n)

    tmp5 <- varest(out$estimates, out$eifs)
    tau <- tmp5$var
    eiftau <- tmp5$eif
    setau <- sqrt(myvar(eiftau) / n)

    variances  <-  c(tau = tau, tau_cm = tcm, tau_eh = teh, tau_em = tem, tau_mv = tmv)
    std.error  <-  c(tau = setau, tau_cm = setcm, tau_eh = seteh, tau_em = setem, tau_mv = setmv)

    perc <- variances / tau * 100
    std.error.pct <- 100 * c(NA,
                       sqrt(myvar(1 / tau * eifcm - tcm / tau^2 * eiftau) / n),
                       sqrt(myvar(1 / tau * eifeh - teh / tau^2 * eiftau) / n),
                       sqrt(myvar(1 / tau * eifem - tem / tau^2 * eiftau) / n),
                       sqrt(myvar(1 / tau * eifmv - tmv / tau^2 * eiftau) / n))

    return(data.frame(variances, std.error, perc, std.error.pct))

}

# Special case for two studies
binary_decomp <- function(out, conf.level = 0.95) {
  ests <- out$estimates
  eifs <- out$eifs
  z <- qnorm((1 - conf.level) / 2 + conf.level)
  n <- dim(eifs)[1]
  
  if (length(ests) != 8){
    stop('must be run with exactly two study levels')
  }
  
  calc_param <- function(index1, index2){
    est <- ests[index1] - ests[index2]
    eifs_diff <- eifs[, index1] - eifs[, index2]
    se <- sd(eifs_diff)/sqrt(n)
    lower <- est - z * se
    upper <- est + z * se
    return(list(est = est, se = se, lower = lower, upper = upper))
  }
  
  dcm <- calc_param(8, 4)
  deh <- calc_param(4, 1)
  dem <- calc_param(4, 3)
  dmv <- calc_param(3, 1)
  
  return(list(dcm = dcm, deh = deh, dem = dem, dmv = dmv))
}

# Clamp function for simulation
clamp = function(a, lower, upper) {
  pmax(pmin(a, upper), lower)
}

# Data generation for simulation
generate = function(fw, fs, fa, fm, fy, n){
  w = fw(n); s = fs(w); a = fa(w, s); m = fm(w, s, a); y = fy(w, s, a, m)
  return(data = data.frame(w0 = w, study = s, a = a, m = m, y = y))
}

# Simulate dcm for simulation check
sim_dcm = function(fw, fs, fa, fm, fy, n = 1e6){
  # calc term 1
  w = fw(n)
  s = fs(w)
  w = w[which(s == 1)]
  s = 1
  a = 1
  m = fm(w, s, a)
  y = fy(w, s, a, m)
  term1 = mean(y)
  
  # calc term 2
  w = fw(n)
  s = fs(w)
  w = w[which(s == 1)]
  s = 1
  a = 0
  m = fm(w, s, a)
  y = fy(w, s, a, m)
  term2 = mean(y)
  
  # calc term 3
  w = fw(n)
  s = fs(w)
  w = w[which(s == 0)]
  s = 1
  a = 1
  m = fm(w, s, a)
  y = fy(w, s, a, m)
  term3 = mean(y)
  
  # calc term 4
  w = fw(n)
  s = fs(w)
  w = w[which(s == 0)]
  s = 1
  a = 0
  m = fm(w, s, a)
  y = fy(w, s, a, m)
  term4 = mean(y)
  
  return(list(dcm = (term1 - term2) - (term3 - term4)))
}

# Simulate deh for simulation check
sim_deh <- function(fw, fs, fa, fm, fy, n = 1e6){
    # calc term 1
    w = fw(n)
    s = fs(w)
    w = w[which(s == 0)]
    s = 1
    a = 1
    m = fm(w, s, a)
    y = fy(w, s, a, m)
    term1 = mean(y)
    
    # calc term 2
    w = fw(n)
    s = fs(w)
    w = w[which(s == 0)]
    s = 1
    a = 0
    m = fm(w, s, a)
    y = fy(w, s, a, m)
    term2 = mean(y)
    
    # calc term 3
    w = fw(n)
    s = fs(w)
    w = w[which(s == 0)]
    s = 0
    a = 1
    m = fm(w, s, a)
    y = fy(w, s, a, m)
    term3 = mean(y)
    
    # calc term 4
    w = fw(n)
    s = fs(w)
    w = w[which(s == 0)]
    s = 0
    a = 0
    m = fm(w, s, a)
    y = fy(w, s, a, m)
    term4 = mean(y)
    
    return(list(deh = (term1 - term2) - (term3 - term4)))
}

# Simulate dem for simulation check
sim_dem = function(fw, fs, fa, fm, fy, n = 1e6){
  # calc term 1
  w = fw(n)
  s = fs(w)
  w = w[which(s == 0)]
  s = 1
  a = 1
  m = fm(w, s, a)
  s = 1
  y = fy(w, s, a, m)
  term1 = mean(y)
  
  # calc term 2
  w = fw(n)
  s = fs(w)
  w = w[which(s == 0)]
  s = 1
  a = 0
  m = fm(w, s, a)
  s = 1
  y = fy(w, s, a, m)
  term2 = mean(y)
  
  # calc term 3
  w = fw(n)
  s = fs(w)
  w = w[which(s == 0)]
  s = 1
  a = 1
  m = fm(w, s, a)
  s = 0
  y = fy(w, s, a, m)
  term3 = mean(y)
  
  # calc term 4
  w = fw(n)
  s = fs(w)
  w = w[which(s == 0)]
  s = 1
  a = 0
  m = fm(w, s, a)
  s = 0
  y = fy(w, s, a, m)
  term4 = mean(y)
  
  return(list(dem = (term1 - term2) - (term3 - term4)))
}

# Simulate dmv for simulation check
sim_dmv = function(fw, fs, fa, fm, fy, n = 1e6){
  # calc term 1
  w = fw(n)
  s = fs(w)
  w = w[which(s == 0)]
  s = 1
  a = 1
  m = fm(w, s, a)
  s = 0
  y = fy(w, s, a, m)
  term1 = mean(y)
  
  # calc term 2
  w = fw(n)
  s = fs(w)
  w = w[which(s == 0)]
  s = 1
  a = 0
  m = fm(w, s, a)
  s = 0
  y = fy(w, s, a, m)
  term2 = mean(y)
  
  # calc term 3
  w = fw(n)
  s = fs(w)
  w = w[which(s == 0)]
  s = 0
  a = 1
  m = fm(w, s, a)
  s = 0
  y = fy(w, s, a, m)
  term3 = mean(y)
  
  # calc term 4
  w = fw(n)
  s = fs(w)
  w = w[which(s == 0)]
  s = 0
  a = 0
  m = fm(w, s, a)
  s = 0
  y = fy(w, s, a, m)
  term4 = mean(y)
  
  return(list(dmv = (term1 - term2) - (term3 - term4)))
}

# Calculate dcm analytically for simulation
calc_dcm <- function(Q, B, C){
  term1 <- 1 + 2*(2*Q^3/3-Q^2+2/3)*C + 1/3*(1 + 2*Q^3/3-Q^2+2/3 + B)
  
  term2 <- 0
  
  term3 <- 1 + 2*(2*Q^3/3+Q^2*(1-Q)-Q^2-Q*(1-Q)^2+Q-2*(1-Q)^3/3+(1-Q)^2)*C + 1/3*(1 + 2*Q^3/3+Q^2*(1-Q)-Q^2-Q*(1-Q)^2+Q-2*(1-Q)^3/3+(1-Q)^2 + B)
  
  term4 <- 0
  
  return(list(dcm = (term1 - term2) - (term3 - term4)))
}

# Calculate deh analytically for simulation
calc_deh <- function(Q, B, C){
  term1 <- 1 + 2*(2*Q^3/3+Q^2*(1-Q)-Q^2-Q*(1-Q)^2+Q-2*(1-Q)^3/3+(1-Q)^2)*C + 1/3*(1 + 2*Q^3/3+Q^2*(1-Q)-Q^2-Q*(1-Q)^2+Q-2*(1-Q)^3/3+(1-Q)^2 + B)
  
  term2 <- 0
  
  term3 <- 1 + 1*(2*Q^3/3+Q^2*(1-Q)-Q^2-Q*(1-Q)^2+Q-2*(1-Q)^3/3+(1-Q)^2)*C + 1/3*(1 + 2*Q^3/3+Q^2*(1-Q)-Q^2-Q*(1-Q)^2+Q-2*(1-Q)^3/3+(1-Q)^2)
  
  term4 <- 0
  
  return(list(deh = (term1 - term2) - (term3 - term4)))
}

# Calculate dem analytically for simulation
calc_dem <- function(Q, B, C){
  term1 <- 1 + 2*(2*Q^3/3+Q^2*(1-Q)-Q^2-Q*(1-Q)^2+Q-2*(1-Q)^3/3+(1-Q)^2)*C + 1/3*(1 + 2*Q^3/3+Q^2*(1-Q)-Q^2-Q*(1-Q)^2+Q-2*(1-Q)^3/3+(1-Q)^2 + B)
  
  term2 <- 0
  
  term3 <- 1 + 1*(2*Q^3/3+Q^2*(1-Q)-Q^2-Q*(1-Q)^2+Q-2*(1-Q)^3/3+(1-Q)^2)*C + 1/3*(1 + 2*Q^3/3+Q^2*(1-Q)-Q^2-Q*(1-Q)^2+Q-2*(1-Q)^3/3+(1-Q)^2 + B)
  
  term4 <- 0
  
  return(list(dem = (term1 - term2) - (term3 - term4)))
}

# Calculate dmv analytically for simulation
calc_dmv <- function(Q, B, C){
  return(list(dmv = calc_deh(Q, B, C)$deh - calc_dem(Q, B, C)$dem))
}



