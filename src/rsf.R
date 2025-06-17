library(randomForestSRC)
library(splines)  # for ns()
library(survival)
library(parallel) # for detectCores()
library(data.table)


# Creating weighted df
W_tuned <- readRDS("W_tuned.rds")
df <- fread("df_clean.csv")

df$wtps <- W_tuned$weights


# set.seed(2025)            
# test <- df[sample(nrow(df), 50), ]

# n - 1 cores
nc <- max(1, parallel::detectCores() - 1)   # leave 1 core free

# weighted RSF with imputation + multicore processing
rsf_wt_imp <- rfsrc(
  Surv(FINALDISCHARGEHRS, WITHDRAWALLST) ~
    VTEPROPHYLAXISHRS +
    TBIHIGHESTTOTALGCS +
    ISS +
    AGEyears +
    Hx_AnticoagulantTherapy +
    ANTIBIOTICTHERAPY_Before_VTE +
    HMRRHGCTRLSURG_Before_VTE +
    Ventilator_Assisted_Day1 +
    ICU_Early_Intervention +
    BLOOD4ML +
    PLATELETS4ML +
    SEX,
  data      = df,
  case.wt   = df$wtps,   #
  ntree     = 1000,
  nsplit    = 10,
  splitrule = "logrank",
  na.action = "na.impute",
  n.cores   = nc          
)


saveRDS(rsf_wt_imp, file="shinyapp/rsf_wt.RDS")