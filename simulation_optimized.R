##############################################################################

# Optimized Trauma Time-Varying Simulation (data.table + parallel, v3)

# – matrix-based hazard lookup, single-threaded workers, lean exports –

##############################################################################

suppressPackageStartupMessages({
library(data.table)
library(splines)
library(survival)
library(ggplot2)
library(patchwork)
library(parallel)
})

set.seed(2025)

#---------------------------------------------------------------------------

# GLOBAL PARAMETERS

#---------------------------------------------------------------------------
n_patients   <- 2000          # feel free to scale
max_hours    <- 168
total_hours  <- 730           # follow-up horizon (≈ 1 month)

# Dense 0–168 h, then every 24 h

hrs_detailed <- 0:168
hrs_sparse   <- seq(192, total_hours, by = 24)
hrs          <- c(hrs_detailed, hrs_sparse)

# Severity, GCS etc.

severity_levels <- c("mild", "moderate", "severe")
severity_probs  <- c(0.30, 0.45, 0.25)
gcs_ranges <- list(mild = 13:15, moderate = 9:12, severe = 3:8)

bleed_plateau_duration <- 30     # h
bleed_baseline         <- 0.001
peak_bleed_by_severity <- c(mild = 0.015, moderate = 0.025, severe = 0.040)

surgery_prob     <- c(mild = 0.15, moderate = 0.40, severe = 0.70)
surgery_timing   <- data.table(severity = severity_levels,
mean = c(48, 24, 6), sd = c(24, 12, 3))
repeat_surgery_prob  <- c(mild = 0.05, moderate = 0.15, severe = 0.30)
repeat_surgery_days  <- c(3, 7, 14, 21) * 24   # h

# VTE & AC

base_vte_haz   <- 0.0008
vte_acc_rate   <- 6e-6
vte_surg_boost <- 2
ac_vte_reduc   <- 0.70
ac_bleed_incr  <- 1.50
ac_thresholds  <- list(platelets = 50, hgb = 8)

# Mortality

mort_base        <- c(mild = 0.0001, moderate = 0.0002, severe = 0.0004)
mort_bleed_mult  <- 2
mort_dvt_mult    <- 3
mort_pe_mult     <- 5

#---------------------------------------------------------------------------

# 1. PRE-COMPUTE BASE HAZARDS (matrix – O(1) lookup)

#---------------------------------------------------------------------------
compute_base_hazards <- function() {
hrs_vec <- 0:total_hours
n_h     <- length(hrs_vec)
n_s     <- length(severity_levels)

bleed_mat <- matrix(bleed_baseline, n_s, n_h)
vte_vec   <- base_vte_haz + vte_acc_rate * hrs_vec

for (s in seq_len(n_s)) {
pk <- peak_bleed_by_severity[s]
plateau <- hrs_vec <= bleed_plateau_duration
trans   <- hrs_vec > bleed_plateau_duration &
hrs_vec <= bleed_plateau_duration + 30
bleed_mat[s, plateau] <- pk
bleed_mat[s, trans] <- pk - (pk - bleed_baseline) *
((hrs_vec[trans] - bleed_plateau_duration) / 30)
# beyond trans already baseline
}
list(bleed = bleed_mat, vte = vte_vec, hrs = hrs_vec)
}
base_haz <- compute_base_hazards()

fast_hazard <- function(sev_idx, h) {           # h ≥ 0, integer
h_idx <- min(h + 1L, ncol(base_haz$bleed))    # bounds check
list(
bleed = base_haz$bleed[sev_idx, h_idx],
vte   = base_haz$vte[h_idx]
)
}

#---------------------------------------------------------------------------

# 2. SINGLE-PATIENT SIMULATION (vectorised + pre-drawn uniforms)

#---------------------------------------------------------------------------
simulate_patient_optimized <- function(patient_id) {

severity <- sample(severity_levels, 1, prob = severity_probs)
sev_idx  <- match(severity, severity_levels)
initial_gcs <- sample(gcs_ranges[[severity]], 1)

#------------------- surgery schedule ------------------------------------
has_surgery <- runif(1) < surgery_prob[severity]
surgery_time <- if (has_surgery) {
pars <- surgery_timing[surgery_timing$severity == severity]
val  <- rnorm(1, pars$mean, pars$sd)
pmax(1, pmin(val, max_hours - 48))
} else NA_real_

#------------------- allocate result table --------------------------------
n_hrs <- length(hrs)
dt <- data.table(
patient_id = patient_id, hour = hrs,
severity, gcs = NA_real_, on_ac = FALSE,
surgery_scheduled = FALSE, had_surgery = FALSE,
time_since_surgery = NA_real_, platelets = NA_real_, hgb = NA_real_,
bleed_hz = NA_real_, dvt_hz = NA_real_, pe_hz = NA_real_,
new_bleed = 0L, new_dvt = 0L, new_pe = 0L, new_vte = 0L,
cum_bleed = 0L, cum_dvt = 0L, cum_pe = 0L,
had_dvt = FALSE, had_pe = FALSE,
alive = TRUE, died = 0L
)

#------------------- pre-draw random numbers ------------------------------
u_bleed <- runif(n_hrs); u_dvt <- runif(n_hrs); u_pe <- runif(n_hrs); u_mort <- runif(n_hrs)
plts_noise <- cumsum(rnorm(n_hrs, 0, 5))
hgb_noise  <- cumsum(rnorm(n_hrs, 0, 0.1))
plts_base  <- switch(severity, mild = 280, moderate = 250, severe = 220)
hgb_base   <- switch(severity, mild = 14, moderate = 13, severe = 12)
u_repeat   <- runif(length(repeat_surgery_days))

#------------------- mutable state ----------------------------------------
on_ac <- FALSE; cum_bleed <- cum_dvt <- cum_pe <- 0L
last_bleed_hour <- NA_real_; had_surgery <- FALSE
surgery_count <- 0L; next_surgery <- NA_real_
had_dvt <- had_pe <- FALSE; alive <- TRUE; death_hour <- NA_real_

#------------------- main hourly loop -------------------------------------
for (i in seq_len(n_hrs)) {
if (!alive) break
h <- hrs[i]

```
# ---- GCS ----
dt[i, gcs := min(initial_gcs + floor(h/24)*2, 15)]

# ---- surgery flags ----
if (!is.na(surgery_time) && h >= surgery_time && !had_surgery) had_surgery <- TRUE
dt[i, had_surgery := had_surgery]
dt[i, surgery_scheduled := !is.na(surgery_time) && h < surgery_time]
if (had_surgery) dt[i, time_since_surgery := h - surgery_time]

# ---- labs ----
dt[i, platelets := pmax(20, pmin(600, plts_base - cum_bleed*30 + plts_noise[i]))]
dt[i, hgb       := pmax(6,  pmin(18, hgb_base  - cum_bleed*1.5 + hgb_noise[i]))]

# ---- AC management (within 0-168 h) ----
ac_allowed <- TRUE
if (!is.na(surgery_time) && h >= surgery_time-12 && h <= surgery_time+24) {
  ac_allowed <- FALSE; on_ac <- FALSE
}
if (!on_ac && ac_allowed && h <= 168) {
  safe_to_start <- dt[i, platelets] > ac_thresholds$platelets &&
                   dt[i, hgb]       > ac_thresholds$hgb &&
                   h >= 1 &&
                   (is.na(last_bleed_hour) || h - last_bleed_hour > 24)
  if (safe_to_start || had_dvt || had_pe) on_ac <- TRUE
}
dt[i, on_ac := on_ac]

# ---- repeat-surgery trigger ----
if (had_surgery && surgery_count < length(repeat_surgery_days) && !is.na(surgery_time)) {
  if ((h - surgery_time) >= repeat_surgery_days[surgery_count+1] &&
      u_repeat[surgery_count+1] < repeat_surgery_prob[severity]) {
    surgery_count <- surgery_count + 1L
    next_surgery  <- h + sample(12:36,1)
  }
}
if (!is.na(next_surgery) && h >= next_surgery) {
  surgery_time  <- next_surgery
  next_surgery  <- NA_real_
}

# ---- hazards (fast lookup) ----
bh <- fast_hazard(sev_idx, as.integer(h))
bleed_hz <- bh$bleed
vte_hz   <- bh$vte

if (!is.na(surgery_time) && h >= surgery_time && h < surgery_time+24) {
  bleed_hz <- bleed_hz + 0.03 * exp(-(h - surgery_time)/12)
}
if (!is.na(surgery_time) && h > surgery_time+48) vte_hz <- vte_hz * vte_surg_boost

dvt_hz <- vte_hz * 0.75; pe_hz <- vte_hz * 0.25
if (had_dvt) pe_hz <- pe_hz * 2

if (on_ac) {
  bleed_hz <- bleed_hz * ac_bleed_incr
  dvt_hz   <- dvt_hz   * (1 - ac_vte_reduc)
  pe_hz    <- pe_hz    * (1 - ac_vte_reduc)
}

dt[i, `:=`(bleed_hz = bleed_hz, dvt_hz = dvt_hz, pe_hz = pe_hz)]

# ---- event simulation ----
if (u_bleed[i] < bleed_hz) {
  cum_bleed <- cum_bleed + 1L; last_bleed_hour <- h; on_ac <- FALSE
  dt[i, new_bleed := 1L]
}
if (!had_dvt && u_dvt[i] < dvt_hz) {
  cum_dvt <- cum_dvt + 1L; had_dvt <- TRUE; on_ac <- on_ac || (ac_allowed && h<=168)
  dt[i, new_dvt := 1L]
}
if (!had_pe && u_pe[i] < pe_hz) {
  cum_pe <- cum_pe + 1L; had_pe <- TRUE; on_ac <- on_ac || (ac_allowed && h<=168)
  dt[i, new_pe := 1L]
}
dt[i, new_vte := as.integer(new_dvt==1L | new_pe==1L)]

dt[i, `:=`(cum_bleed = cum_bleed, cum_dvt = cum_dvt, cum_pe = cum_pe,
           had_dvt = had_dvt, had_pe = had_pe)]

# ---- mortality ----
mort_hz <- mort_base[severity] *
           (mort_bleed_mult^cum_bleed) *
           (if (had_dvt) mort_dvt_mult else 1) *
           (if (had_pe)  mort_pe_mult  else 1)

if (u_mort[i] < mort_hz) { alive <- FALSE; death_hour <- h; dt[i, `:=`(alive=FALSE, died=1L)] }

```

}

if (!alive) dt[hour <= death_hour] else dt
}

#---------------------------------------------------------------------------

# 3. PARALLEL POPULATION SIMULATION

#---------------------------------------------------------------------------
data.table::setDTthreads(1)            # each worker single-threaded

cat("→ Simulating", n_patients, "patients …\n")
t0 <- Sys.time()

n_cores <- min(detectCores() - 1, 10)
cl <- makeCluster(n_cores)

clusterEvalQ(cl, { suppressPackageStartupMessages(library(data.table)) })

# Export only necessary objects to workers

clusterExport(cl, c(
"simulate_patient_optimized", "base_haz", "fast_hazard",
"severity_levels", "severity_probs", "gcs_ranges",
"peak_bleed_by_severity", "surgery_prob", "surgery_timing",
"repeat_surgery_prob", "repeat_surgery_days", "hrs",
"bleed_plateau_duration", "bleed_baseline", "mort_base",
"mort_bleed_mult", "mort_dvt_mult", "mort_pe_mult",
"base_vte_haz", "vte_acc_rate", "vte_surg_boost",
"ac_vte_reduc", "ac_bleed_incr", "ac_thresholds",
"max_hours", "total_hours"
))

patient_dt_list <- parLapply(cl, 1:n_patients, simulate_patient_optimized)
stopCluster(cl)

full_data <- rbindlist(patient_dt_list)
t1 <- Sys.time()
cat("✓ Done in", round(difftime(t1, t0, units="secs"),1), "sec ⇒",
nrow(full_data), "rows\n")

#---------------------------------------------------------------------------

# 4. (Optional) ANALYSIS / VISUALS — keep your existing analysis blocks

#---------------------------------------------------------------------------

# ↓ Example quick sanity: cumulative fatalities

cat("\nDeaths:", full_data[died==1, uniqueN(patient_id)], "of", n_patients, "\n")

# Convert to tibble for compatibility with existing visualization code

library(tidyverse)
full_data_tbl <- as_tibble(full_data)

##############################################################################

# ANALYSIS: Cox Models with thinned data

##############################################################################

cat("\n→ Fitting time-varying Cox models...\n")

# For Cox models, expand back to hourly data for first 168 hours only

cox_data <- full_data[hour <= 168]

# Prepare survival data

surv_data <- cox_data[, .(
patient_id, hour, died, new_vte, new_bleed, on_ac, severity,
had_surgery, gcs, platelets, hgb, alive
)]

surv_data[, tstart := hour]
surv_data[, tstop := shift(hour, type = "lead", fill = max(hour) + 1), by = patient_id]

# Death model with time-varying effects

death_cox <- coxph(
Surv(tstart, tstop, died) ~
ns(hour, 3) * new_vte +          # Time-varying VTE effect
ns(hour, 3) * new_bleed +        # Time-varying bleed effect
on_ac +                          # AC status
severity +                       # Injury severity
had_surgery +                    # Surgery occurred
gcs +                            # Current GCS
platelets + hgb,                 # Labs
data = surv_data
)

# Secondary models for VTE and bleeding (on alive patients only)

surv_data_alive <- surv_data[alive == TRUE]

vte_cox <- coxph(
Surv(tstart, tstop, new_vte) ~
ns(hour, 4) +                    # Time trend
on_ac +                          # AC effect
severity +                       # Injury severity
had_surgery +                    # Surgery occurred
platelets + hgb,                 # Labs
data = surv_data_alive
)

bleed_cox <- coxph(
Surv(tstart, tstop, new_bleed) ~
ns(hour, 4) +                    # Time trend
on_ac +                          # AC effect
severity +                       # Injury severity
had_surgery +                    # Surgery occurred
platelets + hgb,                 # Labs
data = surv_data_alive
)

# Print summaries

cat("\n=== DEATH MODEL RESULTS ===\n")
print(summary(death_cox))

cat("\n=== Key Hazard Ratios ===\n")
cat(sprintf("VTE HR for AC: %.2f (95%% CI: %.2f-%.2f)\n",
exp(coef(vte_cox)["on_acTRUE"]),
exp(confint(vte_cox)["on_acTRUE", 1]),
exp(confint(vte_cox)["on_acTRUE", 2])))
cat(sprintf("Bleed HR for AC: %.2f (95%% CI: %.2f-%.2f)\n",
exp(coef(bleed_cox)["on_acTRUE"]),
exp(confint(bleed_cox)["on_acTRUE", 1]),
exp(confint(bleed_cox)["on_acTRUE", 2])))

##############################################################################

# VISUALIZATION (using tidyverse for compatibility)

##############################################################################

cat("\n→ Creating visualizations...\n")

# Color scheme

col_vte   <- "#0066CC"
col_bleed <- "#CC0000"
col_ac    <- "#00AA44"

# Panel 1: Average hazards over time by severity

hazard_summary <- full_data_tbl %>%
filter(!on_ac, hour <= 168) %>%
group_by(hour, severity) %>%
summarise(
mean_bleed_hz = mean(bleed_hz),
mean_vte_hz = mean(dvt_hz + pe_hz),
.groups = "drop"
)

p1 <- ggplot(hazard_summary %>% filter(severity == "severe"), aes(x = hour)) +
geom_line(aes(y = mean_bleed_hz), color = col_bleed, size = 1.5) +
geom_line(aes(y = mean_vte_hz), color = col_vte, size = 1.5) +
geom_vline(xintercept = c(30, 48), linetype = "dashed", alpha = 0.3) +
annotate("text", x = 30, y = 0.03, label = "Bleed plateau ends", angle = 90, vjust = -0.5, size = 3) +
annotate("text", x = 48, y = 0.002, label = "~Equilibration", angle = 90, vjust = -0.5, size = 3) +
scale_y_log10() +
labs(
title = "A. Baseline Hazards Over Time (Severe Trauma)",
subtitle = "Realistic bleeding pattern with plateau and drop",
x = "Hours from admission",
y = "Hazard (log scale)"
) +
theme_minimal() +
theme(plot.title = element_text(face = "bold"))

# Panel 2: AC uptake by severity

ac_uptake <- full_data_tbl %>%
filter(hour <= 168) %>%
group_by(hour, severity) %>%
summarise(
pct_on_ac = mean(on_ac) * 100,
.groups = "drop"
)

p2 <- ggplot(ac_uptake, aes(x = hour, y = pct_on_ac, color = severity)) +
geom_line(size = 1.2) +
scale_color_manual(values = c("mild" = "#66AA00", "moderate" = "#FF8800", "severe" = "#CC0000")) +
labs(
title = "B. AC Uptake by Injury Severity",
subtitle = "More severe injuries delay AC initiation",
x = "Hours from admission",
y = "Percent on AC",
color = "Severity"
) +
theme_minimal() +
theme(
plot.title = element_text(face = "bold"),
legend.position = "bottom"
)

# Panel 3: Example patient trajectories

example_patients <- full_data_tbl %>%
filter(patient_id %in% c(
sample(unique(patient_id[severity == "mild"]), 1),
sample(unique(patient_id[severity == "moderate"]), 1),
sample(unique(patient_id[severity == "severe"]), 1),
sample(unique(patient_id[had_surgery]), 1)
), hour <= 168)

p3 <- ggplot(example_patients, aes(x = hour)) +
geom_line(aes(y = bleed_hz, color = "Bleeding"), size = 1) +
geom_line(aes(y = dvt_hz + pe_hz, color = "VTE"), size = 1) +
geom_rect(data = filter(example_patients, on_ac),
aes(xmin = hour, xmax = hour + 1, ymin = 0, ymax = Inf),
fill = col_ac, alpha = 0.2) +
geom_vline(data = filter(example_patients, had_surgery) %>%
group_by(patient_id) %>%
slice(1),
aes(xintercept = hour),
linetype = "dashed", color = "orange") +
scale_color_manual(values = c("Bleeding" = col_bleed, "VTE" = col_vte)) +
scale_y_log10() +
facet_wrap(~ paste0("Pt ", patient_id, " (", severity, ")"), scales = "free_y") +
labs(
title = "C. Example Patient Trajectories",
subtitle = "Green shade = on AC, Orange line = surgery",
x = "Hours from admission",
y = "Hazard (log scale)",
color = "Risk"
) +
theme_minimal() +
theme(
plot.title = element_text(face = "bold"),
legend.position = "bottom"
)

# Panel 4: Cumulative outcomes over full follow-up

cumulative_summary <- full_data_tbl %>%
group_by(hour) %>%
summarise(
mean_cum_bleed = mean(cum_bleed),
mean_cum_vte = mean(cum_dvt + cum_pe),
mean_cum_dvt = mean(cum_dvt),
mean_cum_pe = mean(cum_pe),
mortality = 1 - mean(alive),
pct_on_ac = mean(on_ac) * 100,
.groups = "drop"
)

p4 <- ggplot(cumulative_summary, aes(x = hour)) +
geom_line(aes(y = mean_cum_bleed, color = "Bleeding events"), size = 1.2) +
geom_line(aes(y = mean_cum_vte * 10, color = "VTE events (×10)"), size = 1.2) +
geom_line(aes(y = mortality * 100, color = "Mortality (×100)"), size = 1.2) +
scale_color_manual(values = c(
"Bleeding events" = col_bleed,
"VTE events (×10)" = col_vte,
"Mortality (×100)" = "black"
)) +
labs(
title = "D. Cumulative Outcomes",
subtitle = "Population-level event rates over 1 month follow-up",
x = "Hours from admission",
y = "Mean events per patient",
color = "Outcome"
) +
theme_minimal() +
theme(
plot.title = element_text(face = "bold"),
legend.position = "bottom"
)

# Combine plots

combined_plot <- (p1 | p2) / (p3) / (p4) +
plot_annotation(
title = "Optimized Trauma Time-Varying Analysis with Death Outcomes",
subtitle = sprintf("n = %d patients; 1-month follow-up; AC window 1-168h; Runtime: %.1f sec",
n_patients, difftime(t1, t0, units = "secs")),
theme = theme(
plot.title = element_text(size = 16, face = "bold"),
plot.subtitle = element_text(size = 12)
)
)

# Save plot

ggsave("trauma_optimized_analysis.png",
combined_plot,
width = 14, height = 14,
dpi = 300)

##############################################################################

# SUMMARY

##############################################################################

cat("\n=== FINAL SUMMARY ===\n")

# Find equilibration points by severity

equilibration_summary <- hazard_summary %>%
group_by(severity) %>%
mutate(diff = abs(mean_bleed_hz - mean_vte_hz)) %>%
arrange(severity, diff) %>%
slice(1) %>%
select(severity, hour, mean_bleed_hz, mean_vte_hz)

cat("\nBleeding-VTE equilibration by severity:\n")
print(equilibration_summary)

# AC timing summary

ac_timing <- full_data_tbl %>%
filter(on_ac, hour <= 168) %>%
group_by(patient_id) %>%
slice(1) %>%
ungroup() %>%
group_by(severity) %>%
summarise(
median_ac_start = median(hour),
q25_ac_start = quantile(hour, 0.25),
q75_ac_start = quantile(hour, 0.75),
.groups = "drop"
)

cat("\nAC initiation timing by severity:\n")
print(ac_timing)

# Death analysis

death_summary <- full_data_tbl %>%
filter(died == 1) %>%
group_by(patient_id) %>%
slice(1) %>%
ungroup() %>%
summarise(
n_deaths = n(),
median_death_hour = median(hour),
q25_death = quantile(hour, 0.25),
q75_death = quantile(hour, 0.75),
pct_had_vte = mean(had_dvt | had_pe) * 100,
pct_had_bleed = mean(cum_bleed > 1) * 100,
pct_on_ac = mean(on_ac) * 100
)

cat("\nDeath analysis:\n")
cat(sprintf("  Total deaths: %d (%.1f%%)\n", death_summary$n_deaths,
death_summary$n_deaths / n_patients * 100))
cat(sprintf("  Median time to death: %.0f hours (IQR: %.0f-%.0f)\n",
death_summary$median_death_hour, death_summary$q25_death, death_summary$q75_death))
cat(sprintf("  Deaths with VTE: %.1f%%\n", death_summary$pct_had_vte))
cat(sprintf("  Deaths with major bleeding: %.1f%%\n", death_summary$pct_had_bleed))
cat(sprintf("  Deaths on AC: %.1f%%\n", death_summary$pct_on_ac))

cat("\n✓ Analysis complete! Plot saved to: trauma_optimized_analysis.png\n")
cat(sprintf("Total runtime: %.1f seconds\n", difftime(t1, t0, units = "secs")))
