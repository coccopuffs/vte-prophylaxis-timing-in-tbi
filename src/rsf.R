library(data.table)          
library(randomForestSRC)     
library(survival)          
library(missRanger)          # imputation
library(qs)                  # fast read
library(parallel)            # detectCores()

#TODO
# - make ventilator variable better with GCSQ



## eproducibility & Parallel
set.seed(2025)
n_threads   <- max(1, parallel::detectCores() - 1)   # leave 1 core head-room
options(mc.cores = n_threads)
Sys.setenv(OMP_NUM_THREADS = n_threads)               # honour in C/OpenMP code
data.table::setDTthreads(n_threads)

message("Using ", n_threads, " parallel threads.")

## Load data & weights =
df <- fread("df_clean.csv")                     #  ❶ main analytic cohort
W  <- qread("iptw_weightit_full.qs")            #  ❷ WeightIt object (GBM IPTW)

# sanity check: same number of rows and in the same order
stopifnot(nrow(df) == length(W$weights))

df[, iptw := W$weights[]]                       # materialise weights

## optional: inspect weight tails / ESS
weight_summary <- summary(df$iptw)
ess            <- sum(df$iptw)^2 / sum(df$iptw^2)
print(weight_summary)
message(sprintf("Effective sample size (ESS): %0.0f (%.1f%% of N)",
                ess, 100 * ess / nrow(df)))

# I need somehow this model to work when certain variables are positive like TBIMIDLINESHIFT,
# do I also include TOTALGCS
# So I make an indicator variable for TBIMIDLINESHIFT for certain NAs when 

## ══ 4 · Define modelling formula  ============================================
surv_formula <- Surv(FINALDISCHARGEHRS, WITHDRAWALLST) ~
  # Baseline Knowledge
  AGEyears +
  SEX +
  Hx_AnticoagulantTherapy +
  Hx_BleedingDisorder +
  Hx_DisseminatedCancer +
  Hx_CVA +
  Hx_FunctionallyDependent +
  TOTALGCS +
  ISS +
  # Blood Products Received at 4 Hour Mark
  PRBC_units +
  PLT_units +
  FFP_units +
  CRYO_units +
  # Events that are recorded at 24 hours s/p injury
  TBIHIGHESTTOTALGCS +
  TBIMIDLINESHIFT +
  # Events that happened before VTEPROPHYLAXISHRS
  Ventilator_Before_VTE +
  ICU_Before_VTE +
  ED_VTE +
  ANTIBIOTICTHERAPY + # Filtered for before VTEPROPHYLAXISHRS, could be at any time
  HMRRHGCTRLSURG + # same as above
  # Variable of Interest
  VTEPROPHYLAXISHRS
  
  
  
#gcsqintuated -> for vent variable

## ══ 5 · Fit weighted Random-Survival-Forest =================================
rsf_ctl <- rfsrc(
  formula      = surv_formula,
  data         = df,
  ntree        = 2000,          # ↑ trees → stabler survival curves
  mtry         = floor(sqrt(ncol(df) - 3)),  # heuristic; adapt if needed
  nodesize     = 5,             # small → finer splits for large N
  nsplit       = 10,
  splitrule    = "logrank",
  block.size   = 5e3,           # memory vs speed trade-off
  na.action    = "na.impute",   # RF internal imputation
  importance   = "permute",     # permutation-based VIMP
  case.wt      = df$iptw,
  n.cores      = n_threads,
  seed         = 2025
)

## ══ 6 · Model diagnostics ====================================================
# 6.1 Variable importance
vimp <- data.table(variable = names(rsf_ctl$importance),
                   vimp     = rsf_ctl$importance)[order(-vimp)]

# 6.2 OOB C-index
cindex <- rsf_ctl$err.rate[length(rsf_ctl$err.rate)]
message(sprintf("Out-of-bag C-index: %.3f", 1 - cindex))

# save diagnostics to disk for manuscript supplement
fwrite(vimp, "rsf_vimp.tsv", sep = "\t")
saveRDS(rsf_ctl, "shinyapp/rsf_iptw.RDS")

message("✓ RSF model saved  → shinyapp/rsf_iptw.RDS")
message("✓ Variable importance table → rsf_vimp.tsv")