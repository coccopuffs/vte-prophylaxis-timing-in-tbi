library(parallel)      # to detect available cores
library(data.table)
library(missRanger)
library(WeightIt)
library(cobalt)
library(qs)
library(ggplot2)
library(cowplot)
library(gridExtra)

## -----------------------------------------------------------------------
## 0 · PARALLEL SETTINGS  ––––– choose N = max(1, total_cores - 1)
## -----------------------------------------------------------------------
n_threads <- max(1, parallel::detectCores() - 1)
message("Using ", n_threads, " parallel threads out of ",
        parallel::detectCores(), " available.")

setDTthreads(n_threads)                    # data.table
options(mc.cores = n_threads)              # generic parallel option
Sys.setenv("OMP_NUM_THREADS" = n_threads)  # OpenMP (ranger, gbm, etc.)

set.seed(42)

## -----------------------------------------------------------------------
## 1 · READ DATA
## -----------------------------------------------------------------------
df <- fread("df_weighting.csv")                      # data.table auto-threads

## -----------------------------------------------------------------------
## 2 · RANDOM-FOREST + PMM IMPUTATION  (missRanger, OpenMP)
## -----------------------------------------------------------------------
df_imp <- missRanger(
  df,
  num.trees   = 1000,
  pmm.k       = 5,
  maxiter     = 1,
  num.threads = n_threads,                 # <—— PARALLEL
  verbose     = 1,
  seed        = 42,
  respect.unordered.factors = "partition"
)

## -----------------------------------------------------------------------
## 3 · ADD MISSINGNESS INDICATORS
## -----------------------------------------------------------------------
na_flags <- df[ , lapply(.SD, function(x) as.integer(is.na(x)))]
setnames(na_flags, paste0(names(df), "__NA"))
df_ready <- cbind(df_imp, na_flags)

## -----------------------------------------------------------------------
## 4 · GBM WEIGHTING  (WeightIt / twang)
## -----------------------------------------------------------------------
grid <- list(
  n.trees           = c(1500, 2500, 3500),
  shrinkage         = c(0.01, 0.02),
  interaction.depth = 1:3,
  bag.fraction      = c(0.6, 0.7)
)

W <- weightit(
  VTEPROPHYLAXISHRS ~ .,
  data        = df_ready,
  method      = "gbm",
  density     = "kernel",
  estimand    = "ATE",
  criterion   = "p.max",
  tune        = grid,
  stop.method = "es.mean",
  stabilize   = TRUE,
  trim.at     = 0.99,
  na.action   = na.pass,
  n.cores     = n_threads,                # <—— PARALLEL
  parallel.method = "multicore",
  verbose     = TRUE
)

## -----------------------------------------------------------------------
## 5 · SAVE WEIGHTIT OBJECT (compact, high compression)
## -----------------------------------------------------------------------
qsave(W, "iptw_weightit_full.qs", preset = "high")

## -----------------------------------------------------------------------
## 6 · DIAGNOSTIC PLOTS (unchanged except for balance table call)
## -----------------------------------------------------------------------
bal_res <- bal.tab(W, un = TRUE, disp.v.ratio = TRUE)

love_plot <- plot(W, type = "balance") +
  theme_minimal(base_size = 12) +
  ggtitle("Covariate Balance (Love Plot)")

ps_plot <- plot(W, type = "density") +
  theme_minimal(base_size = 12) +
  ggtitle("Propensity Score Overlap")

weight_df <- data.frame(weight = W$weights)
weight_plot <- ggplot(weight_df, aes(weight)) +
  geom_histogram(bins = 50, alpha = 0.8) +
  coord_cartesian(xlim = c(0, quantile(weight_df$weight, 0.995))) +
  labs(title = "Weight Distribution",
       x = "Stabilised Kernel Weight",
       y = "Count") +
  theme_minimal(base_size = 12)

tbl_grob <- tableGrob(
  bal_res$Balance[1:min(25, nrow(bal_res$Balance)), ],
  rows = NULL, theme = ttheme_minimal(base_size = 10)
)

row1 <- plot_grid(love_plot, ps_plot, ncol = 2, rel_widths = c(1, 1))
row2 <- plot_grid(weight_plot, tbl_grob, ncol = 2, rel_widths = c(1, 1.1))
full_fig <- plot_grid(row1, row2, nrow = 2)

pdf("supplementary_iptw_diagnostics.pdf", width = 11, height = 8.5)
grid::grid.draw(full_fig)
dev.off()

message("✓ Saved model  → iptw_weightit_full.qs")
message("✓ Saved figure → supplementary_iptw_diagnostics.pdf")