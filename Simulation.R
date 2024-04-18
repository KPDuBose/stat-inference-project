## Read required packages ####
# installpackages("RSpectra")
# devtools::install_github("variani/lme4qtl")
# installpackages("pbmcapply")

pkgs <- c('dplyr',
          'Matrix',
          'tictoc',
          'RSpectra',
          'parallel',
          'lme4qtl',
          'lme4',
          'data.table',
          'pbmcapply')
invisible(sapply(pkgs, library, character.only = TRUE))
`%+%` <- function(...) paste0(...)
##
## Set global parameters ####

### Name of simulation run/setting ####
name <- 'test'

### Family types ####
a_1 <- Matrix(data = c(1.0, 0.0, 0.5,
                       0.0, 1.0, 0.5,
                       0.5, 0.5, 1.0),
              nrow = 3,
              ncol = 3,
              byrow = TRUE)
a_2 <- Matrix(data = c(1.0, 0.0, 0.5, 0.5,
                       0.0, 1.0, 0.5, 0.5,
                       0.5, 0.5, 1.0, 0.5,
                       0.5, 0.5, 0.5, 1.0),
              nrow = 4,
              ncol = 4,
              byrow = TRUE)
a_3 <- Matrix(data = c(1.0, 0.0, 0.5, 0.5, 0.5,
                       0.0, 1.0, 0.5, 0.5, 0.5,
                       0.5, 0.5, 1.0, 0.5, 0.5,
                       0.5, 0.5, 0.5, 1.0, 0.5,
                       0.5, 0.5, 0.5, 0.5, 1.0),
              nrow = 5,
              ncol = 5,
              byrow = TRUE)
a_4 <- Matrix(data = c(1.0, 0.0, 0.5, 0.5, 0.5, 0.5,
                       0.0, 1.0, 0.5, 0.5, 0.5, 0.5,
                       0.5, 0.5, 1.0, 0.5, 0.5, 0.5,
                       0.5, 0.5, 0.5, 1.0, 0.5, 0.5,
                       0.5, 0.5, 0.5, 0.5, 1.0, 0.5,
                       0.5, 0.5, 0.5, 0.5, 0.5, 1.0),
              nrow = 6,
              ncol = 6,
              byrow = TRUE)
as <- lapply(list(a_1,a_2,a_3,a_4),
             function(x) t(chol(x)))

### Number of random normal fixed effects ####
p <- 1

### Fixed effect coefficients and names ####
beta <- rep(1, p)
xnames <- 'x' %+% 1:p

### Weights for family types ####
ws <- c(0.2,0.3,0.3,0.2)

### Number of families to choose: (10,50,100,200,400) ####
K <- 100

### Error variance ####
s2E <- 1

### Heritability: (0.0, 0.1, 0.2, 0.5) ####
h2 <- 0.0

### Alphas (significance level) to consider ####
alphas <- c(0.05, 0.01)

### Number of samples for exact distribution ####
samps_exact <- 10000

### How many cores to use to calc exact distribution ####
cores_exact <- 1

### Number of simulations ####
nsim <- 100

### Number of cores for simulation ####
cores_sim <- 1

### Statistics to calculate ####
statistics <- c('LRT', 'RLRT')# 'Q', 'Qtilde', 'Fstat')

### Fit linear model ####
##-- X does NOT include the intercept, one is automatically added
##-- corr MUST have row and column names for subject ID (same order as Y and X)
fit_lmm <- function(Y,
                    X,
                    corr = NULL,
                    REML) {
  if(!is.null(corr)) {
    d <- data.frame('y' = as.matrix(Y),
                    as.matrix(X),
                    id = colnames(corr))
    model <- relmatLmer(y ~ . - id + (1|id),
                        d,
                        relmat = list(id = corr),
                        REML = REML)
    beta <- getME(model, name = "fixef")
    beta <- beta[order(names(beta))]

    foo <- as.data.frame(VarCorr(model))[, c("var1", "vcov")] %>%
      mutate(var1 = case_when(is.na(var1) ~ "s2E",
                              TRUE ~ 's2A'))
    vcs <- setNames(foo$vcov, foo$var1)
    vcs <- vcs[order(names(vcs))]
    ll_REML <- logLik(model, REML = TRUE)[1]
    ll_ML <- logLik(model, REML = FALSE)[1]
    rss <- unname(getME(model, "devcomp")$cmp["pwrss"])
    list("type" = "lme4",
         "model" = model,
         "beta" = beta,
         "vcs" = vcs,
         "ll_REML" = ll_REML,
         "ll_ML" = ll_ML,
         "rss" = rss,
         "method" = ifelse(REML, "REML", "ML"))
  } else {
    d <- data.frame('y' = as.matrix(Y),
                    as.matrix(X))
    model <- lm(y ~ .,
                d)
    X <- model.matrix(model)
    beta <- model$coefficients
    r <- Y - X %*% beta
    rss <- {t(r) %*% r}[1, 1]
    if (!REML) {
      vcs <- c("s2E" = rss / nrow(X))
    } else {
      vcs <- c("s2E" = rss / (nrow(X) - ncol(X)))
    }
    lnum <- log(2 * pi * rss)
    d <- 0
    dd <- det(t(X) %*% X)
    ll_REML <- logLik(model, REML = TRUE)[1]
    ll_ML <- logLik(model, REML = FALSE)[1]
    beta <- beta[order(names(beta))]
    list("type" = "lm",
         "model" = model,
         "beta" = beta,
         "vcs" = vcs,
         "rss" = rss,
         "ll_REML" = ll_REML,
         "ll_ML" = ll_ML,
         "method" = ifelse(REML, "REML", "ML"))
  }
}

### Get critical values of 50:50 0 and Chi Square 1 ####
qmix <- function(q) {
  if (q < 0.5) {
    0
  } else {
    qchisq(2*(q - 0.5), 1)
  }
}
qmix <- Vectorize(qmix)

### Get time ####
get_time <- function() {
  format(Sys.time(), "%Y_%m_%d %I_%M_%S_%p")
}

## Create functions for simulation ####
### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Generate what variance components should be ####
gen_vc <- function(h2,
                   s2E) {
  s2A <- h2 / (1 - h2) * s2E
  list(h2 = h2,
       s2A = s2A,
       s2E = s2E)
}

### Generate data ####
generate_data <- function(K,
                          p,
                          ws,
                          vc,
                          beta) {
  fam_types <- sample(1:length(as),
                      K,
                      TRUE,
                      ws)
  a <- lapply(fam_types,
              function(x) as[[x]])
  A <- bdiag(a)
  AAt <- tcrossprod(A)
  n <- nrow(A)
  X <- Matrix(rbinom(n*p,1,0.5),
              ncol = p)
  P <- diag(n) - tcrossprod(qr.Q(qr(X)))
  rho <- eigen(t(A) %*% P %*% A,
               symmetric = TRUE,
               only.values = TRUE)$values

  phi <- eigen(AAt,
               symmetric = TRUE,
               only.values = TRUE)$values
  a <- rnorm(n, sd = sqrt(vc$s2A))
  e <- rnorm(n, sd = sqrt(vc$s2E))
  Y <- {X%*%beta + A %*% a + e}[, 1]
  colnames(X) <- xnames
  ids <- 'id'%+%{1:n}
  colnames(AAt) <- ids
  rownames(AAt) <- ids
  list(A = A,AAt = AAt,
       X = X,
       P = P,
       rho = rho,
       phi = phi,
       Y = Y,
       n = n,
       p = p)
}

### Get exact distributions ####
#- just one
get_exact_dist_one <- function(data,
                               statistics) {
  out <- list()
  ##-- common things
  u2 <- rnorm(data$n)^2
  logRSS0 <- log(sum(u2))
  ##-- LRT
  if('LRT' %in% statistics) {
    f_n <- function(lambda) {
      onepluslambdarho <- 1 + lambda * data$rho
      onepluslambdaphi <- 1 + lambda * data$phi
      data$n * log(sum(u2/onepluslambdarho)) + sum(log(onepluslambdaphi))
    }
    out['LRT'] <- data$n * logRSS0 -
      optim(par = 0.5,
            fn = f_n,
            lower = 0,
            method = "L-BFGS-B")$value
  }

  ##-- RLRT
  if ('RLRT' %in% statistics) {
    g_n <- function(lambda) {
      onepluslambdarho <- 1 + lambda * data$rho
      (data$n - data$p) * log(sum(u2/onepluslambdarho)) + sum(log(onepluslambdarho))
    }
    out["RLRT"] <- (data$n - data$p) * logRSS0 -
      optim(par = 0.5,
            fn = g_n,
            lower = 0,
            method = "L-BFGS-B")$value
  }
  ##-- Fstat
  if ("Fstat" %in% statistics) {
    r_n <- function(lambda) {

      onepluslambdarho <- 1 + lambda * data$rho
      (data$n - data$p) * log(sum(u2/onepluslambdarho)) + sum(log(onepluslambdarho))
    }
    lambda_hat <- optim(par = 0.5,
                        fn = r_n,
                        lower = 0,
                        method = "L-BFGS-B")$par
    lambda_hatrho <- lambda_hat * data$rho
    onepluslambda_hatrho <- 1 + lambda_hatrho
    out["Fstat"] <- sum(u2 * lambda_hatrho / onepluslambda_hatrho) / sum(u2 / onepluslambda_hatrho)
  }
  ##-- Q
  if ("Q" %in% statistics) {
    out["Q"] <- sum(data$rho * u2) / sum(u2)
  }
  ##-- Qtilde
  if ("Qtilde" %in% statistics) {
    out["Qtilde"] <- sum(data$rho * u2)
  }

  tibble(type = 'exact',
         statistic = names(out),
         value = unlist(out))
}

#- Get many
get_exact_dist <- function(data,
                           statistics,
                           samps_exact,
                           cores_exact) {
  bind_rows(mclapply(1:samps_exact,
                     function(x) {
                       o <- get_exact_dist_one(data = data,
                                               statistics = statistics)
                       o$sim <- x
                       o
                     },
                     mc.cores = cores_exact))
}

### Calculate test statistics ####
get_statistics <- function(data,
                           statistics) {
  stats <- list()
  ests <- list()
  ##-- LRT
  if('LRT' %in% statistics) {
    ###--- Fit alternative model

    ml_model_A <- fit_lmm(Y = data$Y,
                          X = data$X,
                          corr = data$AAt,
                          REML = FALSE)
    ests[['ml_A']] <- bind_rows(tibble(var = names(ml_model_A$vcs),
                                       value = ml_model_A$vcs),
                                tibble(var = 'beta_' %+% names(ml_model_A$beta),
                                       value = ml_model_A$beta)) %>%
      mutate(method = 'ml',
             hypoth = 'hA')

    ###--- Fit null model
    ml_model_0 <- fit_lmm(Y = data$Y,
                          X = data$X,
                          REML = FALSE)
    ests[['ml_0']] <- bind_rows(tibble(var = names(ml_model_0$vcs),
                                       value = ml_model_0$vcs),
                                tibble(var = 'beta_' %+% names(ml_model_0$beta),
                                       value = ml_model_0$beta)) %>%
      mutate(method = 'ml',
             hypoth = 'h0') %>%
      bind_rows(tibble(hypoth = rep(c('h0', 'hA'), each = 2),
                       var = rep(c('ll', 'rll'), 2),
                       method = 'ml',
                       value = c(ml_model_0$ll_ML,
                                 ml_model_0$ll_REML,
                                 ml_model_A$ll_ML,
                                 ml_model_A$ll_REML)))
    stats["LRT"] <- 2 * (ml_model_A$ll_ML - ml_model_0$ll_ML)
  }

  ##-- Things that depend on restricted likelihood (could be optimized)
  if(any(c('RLRT', 'Fstat', 'Q', 'Qtilde') %in% statistics)) {
    ###--- Fit alternative model
    rml_model_A <- fit_lmm(Y = data$Y,
                           X = data$X,
                           corr = data$AAt,
                           REML = TRUE)
    ests[['rml_A']] <- bind_rows(tibble(var = names(rml_model_A$vcs),
                                        value = rml_model_A$vcs),
                                 tibble(var = 'beta_' %+% names(rml_model_A$beta),
                                        value = rml_model_A$beta)) %>%
      mutate(method = 'rml',
             hypoth = 'hA')
    ###--- Fit null model
    rml_model_0 <- fit_lmm(Y = data$Y,
                           X = data$X,
                           REML = TRUE)
    ests[['rml_0']] <- bind_rows(tibble(var = names(rml_model_0$vcs),
                                      value = rml_model_0$vcs),
                                 tibble(var = 'beta_' %+% names(rml_model_0$beta),
                                        value = rml_model_0$beta)) %>%
      mutate(method = 'rml',
             hypoth = 'h0') %>%
      bind_rows(tibble(hypoth = rep(c('h0', 'hA'), each = 2),
                       var = rep(c('ll', 'rll'), 2),
                       method = 'rml',
                       value = c(rml_model_0$ll_ML,
                                 rml_model_0$ll_REML,
                                 rml_model_A$ll_ML,
                                 rml_model_A$ll_REML)))
    ###--- Calculate some things that are re-rused (could be optimized)
    r <- Matrix(rml_model_0$model$residuals, ncol = 1)
    Qs_num <- (t(r) %*% data$AAt %*% r)
    ###--- Q
    if ('Q' %in% statistics) {
      stats["Q"] <- {Qs_num / (t(r) %*% r)}[1, 1]
    }
    ###--- Qtilde
    if ("Qtilde" %in% statistics) {
      stats["Qtilde"] <- {Qs_num}[1, 1]
    }
    ###--- RLRT
    if ("RLRT" %in% statistics) {
      stats["RLRT"] <- 2 * (rml_model_A$ll_REML - rml_model_0$ll_REML)
    }
    ###--- Fstat
    if ("Fstat" %in% statistics) {
      stats["Fstat"] <- (rml_model_0$rss - rml_model_A$rss) / rml_model_A$rss
    }
  }
  ##-- output
  list('stats' = tibble(statistic = names(stats),
                        value = unlist(stats)),
       'ests' = bind_rows(ests))
}

one_simulation <- function(sim,
                           K,
                           p,
                           ws,
                           vc,
                           beta,
                           statistics,
                           samps_exact,
                           cores_exact,
                           crit_asym,
                           alphas) {
  ##-- Generate data
  data <- generate_data(K = K,
                        p = p,
                        ws = ws,
                        vc = vc,
                        beta = beta)
  ##-- Get exact distribution
  dist_exact <- get_exact_dist(data = data,
                               statistics = statistics,
                               samps_exact = samps_exact,
                               cores_exact = cores_exact)
  ##-- Get critical valuees
  crit <- dist_exact %>%
    group_by(type, statistic) %>%
    summarise(alpha = alphas,
              crit = quantile(value,
                              1 - alphas),
              .groups = 'keep') %>%
    ungroup %>%
    bind_rows(crit_asym)

  ##-- get statistics
  res <- get_statistics(data = data,
                        statistics = statistics)
  ##-- get reject
  res$reject <- res$stats %>%
    left_join(crit,
              by = 'statistic') %>%
    mutate(reject = value > crit) %>%
    select(statistic, type, alpha, reject)
  ##-- return
  res$ests <- res$ests %>%
    mutate(sim = sim)
  res$stats <- res$stats %>%
    mutate(sim = sim)
  res$reject <- res$reject %>%
    mutate(sim = sim)
  #-- record sample size:
  res$ssize<- c(data$n, sim)
  # %>% mutate(sim = sim)
  res
}

sims <- function(name,
                 nsim,
                 cores_sims,
                 K,
                 p,
                 ws,
                 s2E,
                 h2,
                 beta,
                 statistics = statistics,
                 samps_exact,
                 cores_exact,
                 alphas) {
  time <- get_time()
  ##-- Generate VCs
  vc <- gen_vc(h2 = h2,
               s2E = s2E)
  ##-- Get critical values
  crit_asym <- tibble(type = 'asym',
                      alpha = alphas,
                      crit = qmix(1 - alphas)) %>%
    full_join(tibble(statistic = c('RLRT', 'LRT')),
              by = character()) %>%
    filter(statistic %in% statistics)
  ##-- do simulation
  o <- pbmclapply(1:nsim,
                  one_simulation,
                  K = K,
                  p = p,
                  ws = ws,
                  vc = vc,
                  beta = beta,
                  statistics = statistics,
                  samps_exact = samps_exact,
                  cores_exact = cores_exact,
                  crit_asym = crit_asym,
                  alphas = alphas,
                  mc.cores = cores_sim,
                  mc.style = 'txt')
  stats <- bind_rows(lapply(o, function(x) x$stats))
  rej <- bind_rows(lapply(o, function(x) x$rej))
  ests <- bind_rows(lapply(o, function(x) x$ests))
  pow <- rej %>%
    group_by(statistic, type, alpha) %>%
    summarize(power = 100*mean(reject)) %>%
    ungroup
  #-- record sample size:
  ssize <- bind_rows(lapply(o, function(x) x$ssize))
  ##-- write output
  fwrite(stats,
         file = name %+% "_" %+% h2 %+% "_" %+% 'stats' %+% "_" %+% K %+% ".csv")
  fwrite(rej,
         file = name %+% "_" %+% h2 %+% "_" %+% 'rej' %+% "_" %+% K %+% ".csv")
  fwrite(ests,
         file = name %+% "_" %+% h2 %+% "_" %+% 'ests' %+% "_" %+% K %+% ".csv")
  fwrite(pow,
         file = name %+% "_" %+% h2 %+% "_" %+% 'pow' %+% "_" %+% K %+% ".csv")
  fwrite(ssize,
         file = name %+% "_" %+% h2 %+% "_" %+% 'ssize' %+% "_" %+% K %+% ".csv")
  ##-- save settings
  settings <- list('name' = name,
                   'nsim' = nsim,
                   'K' = K,
                   'p' = p,
                   'ws' = ws,
                   's2E' = s2E,
                   'h2' = h2,
                   'beta' = beta,
                   'statistics' = statistics,
                   'samps_exact' = samps_exact,
                   'as' = as)
  saveRDS(settings,
          name %+% "_" %+% K %+% ".settings")
}
sims(name = name,
     nsim = nsim,
     cores_sims = cores_sim,
     K = K,
     p = p,
     ws = ws,
     s2E = s2E,
     h2 = h2,
     beta = beta,
     statistics = statistics,
     samps_exact = samps_exact,
     cores_exact = cores_exact,
     alphas = alphas)
