library(brms)
library(posterior)
library(tidybayes)
library(tidyverse)


# Stan model code
stan_code <- "
data {
  int<lower=1> N;                    // rows
  int<lower=1> K;                    // # extra covariates
  vector[N]  pulse;                  // NOPULSE
  matrix[N,K] Z;                     // age, HTN, DM, EF, LAP, RAP, preCR
  vector[N]  y;                      // deltaCr24hr
}
parameters {
  real  alpha;
  real  beta_pulse;                  // slope before hinge
  real  delta;                       // slope change after cp
  vector[K] gamma;                   // coefficients for Z
  real<lower=30, upper=70> cp;       // *** change-point is now a parameter ***
  real<lower=0> sigma;
}
transformed parameters {
  vector[N] mu;
  for (n in 1:N)
    mu[n] = alpha
            + beta_pulse * pulse[n]
            + delta * fmax(0, pulse[n] - cp)
            + Z[n] * gamma;
}
model {
  alpha       ~ normal(0, 5);
  beta_pulse  ~ normal(0, 0.1);
  delta       ~ normal(0, 0.1);
  gamma       ~ normal(0, 0.5);
  cp          ~ uniform(30, 70);          // flat prior; edit if needed
  sigma       ~ student_t(3, 0, 0.1);
  y           ~ normal(mu, sigma);
}
generated quantities {
  vector[N] y_rep;
  for (n in 1:N)
    y_rep[n] = normal_rng(mu[n], sigma);
}
"

stan_file <- write_stan_file(stan_code, basename = "hinge_cov_estcp")
mod       <- cmdstan_model(stan_file)


make_standata_est <- function(dat) {
  dat <- dat %>% mutate(across(c(Hypertension, Diabetes), as.numeric))
  Z   <- as.matrix(dat[, c("age","Hypertension","Diabetes",
                           "EF","LAP","RAP","preCR")])
  list(
    N     = nrow(dat),
    K     = ncol(Z),
    pulse = dat$NOPULSE,
    Z     = Z,
    y     = dat$deltaCr24hr
  )
}

imputed_list <- complete(imp, "all")

fits <- map(imputed_list, \(dat)
            mod$sample(
              data            = make_standata_est(dat),
              chains          = 4,
              parallel_chains = 4,
              iter_warmup     = 1000,
              iter_sampling   = 1000,
              adapt_delta     = 0.95,
              refresh         = 0           # silent sampling
            )
)


draws <- do.call(                       # keeps the draws_* class + metadata
  bind_draws,
  c(
    unname(lapply(fits, \(f) f$draws())),   # list of draws_array
    list(along = "chain")                   # stack along chain dimension
  )
)

draws <- as_draws_df(draws)     


summary_tbl <- summarise_draws(draws) |>          # summarise everything …
  filter(variable %in% c("cp", "beta_pulse",      # … then keep the 4 rows
                         "delta", "sigma"))

print(summary_tbl)

cp_med <- summary_tbl %>%               # posterior-median cut-point
  filter(variable == "cp") %>% pull(median)


## 3·1  grid of NOPULSE values
pulse_vals <- seq(
  min(df$NOPULSE, na.rm = TRUE),
  max(df$NOPULSE, na.rm = TRUE),
  length.out = 150
)

## 3·2  cohort-average covariate vector (age, HTN, DM, EF, LAP, RAP, preCR)
Zbar <- df %>%
  select(age, Hypertension, Diabetes, EF, LAP, RAP, preCR) %>%
  mutate(across(c(Hypertension, Diabetes), as.numeric)) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
  as.numeric()

## 3·3  grab draws as a matrix for fast algebra
dm <- as_draws_matrix(draws, variables = c(
  "alpha", "beta_pulse", "delta",
  paste0("gamma[", 1:7, "]"), "cp"
))

alpha      <- dm[, "alpha"      ]
beta_pulse <- dm[, "beta_pulse" ]
delta      <- dm[, "delta"      ]
gamma_mat  <- dm[, grep("^gamma\\[", colnames(dm))]
cp_draw    <- dm[, "cp"         ]

## 3·4  posterior predictions at each pulse value (rows = pulses, cols = draws)
mu_mat <- sapply(
  pulse_vals,
  \(x) alpha +
    beta_pulse * x +
    delta * pmax(0, x - cp_draw) +
    gamma_mat %*% Zbar          # average over covariates
)

## 3·5  summarise to median & 95 % CrI

## mu_mat is (draws × pulses): 80 000 × 150
curve_df <- tibble(
  NOPULSE = pulse_vals,
  mu      = colMedians(mu_mat),                       # posterior median
  .lower  = colQuantiles(mu_mat, probs = 0.025),      # 2.5 %
  .upper  = colQuantiles(mu_mat, probs = 0.975)       # 97.5 %
)



## optional: where you want the blue dotted target
target_y <- 0.30    # ΔCr = 0·30 mg/dL

dat <- df %>%                                   # or use one imputed set
  select(NOPULSE, deltaCr24hr) %>%              # keep the two columns
  drop_na()      

ggplot() +
  ## raw data points
  geom_point(data = dat, aes(NOPULSE, deltaCr24hr),
             colour = "black", size = 1.6) +
  
  ## 95 % CrI ribbon
  geom_ribbon(data = curve_df,
              aes(x = NOPULSE, ymin = .lower, ymax = .upper),
              fill = "grey80", alpha = 0.45) +
  
  ## fitted posterior-median line
  geom_line(data = curve_df, aes(x = NOPULSE, y = mu),
            colour = "#8B0000", linewidth = 0.9) +
  
  ## target horizontal reference
  geom_hline(yintercept = target_y,
             colour = "#1f78b4", linetype = "dotted", linewidth = 0.6) +
  
  ## change-point vertical reference
  geom_vline(xintercept = cp_med,
             colour = "black", linetype = "dashed", linewidth = 0.6) +
  
  ## labels & title
  labs(
    title = "Figure 4. Bayesian Change-Point Linear Regression",
    x     = "NO.PULSE",
    y     = expression(Delta*'Cr (mg/dL, 24 h)')
  ) +
  
  ## theme tweaks to mimic the screenshot
  theme_bw(base_size = 13) +
  theme(
    plot.title      = element_text(face = "plain", hjust = 0,
                                   size = 14, margin = margin(b = 10)),
    panel.grid.major = element_line(colour = "grey88", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    axis.title.x    = element_text(margin = margin(t = 6)),
    axis.title.y    = element_text(margin = margin(r = 6))
  )


cp_draws <- draws$cp                       # numeric vector length = draws

cp_ci <- tibble(
  cp_median = median(cp_draws),
  cp_lower  = quantile(cp_draws, 0.025),
  cp_upper  = quantile(cp_draws, 0.975)
)

print(cp_ci, digits = 2)


r2_one <- function(dat, fit) {
  # ── 0 · tidy covariates & response ────────────────────────────────
  dat <- dat |>
    mutate(across(c(Hypertension, Diabetes), as.numeric))
  
  Z      <- as.matrix(dat[, c("age","Hypertension","Diabetes",
                              "EF","LAP","RAP","preCR")])   # N × 7
  pulse  <- dat$NOPULSE                                    # length N
  N      <- length(pulse)
  
  # ── 1 · pull posterior draws we need ──────────────────────────────
  dm <- as_draws_matrix(fit$draws(
    variables = c("alpha","beta_pulse","delta","cp","sigma",
                  paste0("gamma[", 1:7, "]"))
  ))
  
  alpha_vec <- as.numeric(dm[,"alpha"])        # length = Ndraw
  beta_vec  <- as.numeric(dm[,"beta_pulse"])
  delt_vec  <- as.numeric(dm[,"delta"])
  cp_vec    <- as.numeric(dm[,"cp"])
  sigma_vec <- as.numeric(dm[,"sigma"])        # ensure plain vector
  gamma_mat <- dm[, grep("^gamma\\[", colnames(dm)), drop = FALSE]   # draws × 7
  
  Ndraw <- length(alpha_vec)
  
  # ── 2 · build draw × row matrices ---------------------------------
  pulse_mat <- matrix(rep(pulse,  each = Ndraw), nrow = Ndraw)       # draws × N
  cp_mat    <- matrix(cp_vec,  nrow = Ndraw, ncol = N)               # draws × N
  
  mu_mat <-                                            # draws × N
    matrix(alpha_vec, nrow = Ndraw, ncol = N)             +        # α
    matrix(beta_vec , nrow = Ndraw, ncol = N) * pulse_mat  +       # β·x
    matrix(delt_vec , nrow = Ndraw, ncol = N) *
    pmax(0, pulse_mat - cp_mat)                              +   # δ·hinge
    gamma_mat %*% t(Z)                                           # Zγ
  
  # ── 3 · one y_rep per draw ----------------------------------------
  eps <- matrix(rnorm(Ndraw * N), nrow = Ndraw)
  eps <- sweep(eps, 1, sigma_vec, "*")              # row-wise scale by σ
  yrep_mat <- mu_mat + eps
  
  # ── 4 · Gelman R² per draw ----------------------------------------
  var_mu   <- matrixStats::rowVars(mu_mat)
  var_eps  <- matrixStats::rowVars(yrep_mat - mu_mat)
  var_mu / (var_mu + var_eps)          # returns vector length = Ndraw
}


R2_draws <- map2(imputed_list, fits, r2_one) |> unlist()

R2_tbl <- tibble(
  R2_mean   = mean(R2_draws),
  R2_median = median(R2_draws),
  R2_lower  = quantile(R2_draws, 0.025),
  R2_upper  = quantile(R2_draws, 0.975)
)

print(R2_tbl, digits = 3)