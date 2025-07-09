##############################################################################

# Improved Trauma-Only Time-Varying Simulation

# Realistic bleeding patterns, GCS-based severity, and dynamic surgery

##############################################################################

suppressPackageStartupMessages({
library(tidyverse)
library(splines)
library(survival)
library(ggplot2)
library(patchwork)
})

set.seed(2025)

##############################################################################

# PARAMETERS

##############################################################################

n_patients   <- 500  # Reduced for faster testing
max_hours    <- 168
hrs          <- 0:730


# GCS-based severity levels

severity_levels <- c("mild", "moderate", "severe")
severity_probs  <- c(0.30, 0.45, 0.25)

# Initial GCS ranges by severity

gcs_ranges <- list(
mild     = c(13, 15),
moderate = c(9, 12),
severe   = c(3, 8)
)

# Bleeding parameters by severity

bleed_plateau_duration <- 30  # Hours of high bleeding post-trauma
bleed_baseline <- 0.001       # Background bleeding after acute phase

# Death risk multipliers

mort_base <- c(mild = 0.0001, moderate = 0.0002, severe = 0.0004)
mort_bleed_mult <- 2.0
mort_dvt_mult <- 3.0

mort_pe_mult <- 5.0

# Peak bleeding rates during plateau phase

peak_bleed_by_severity <- c(
mild     = 0.015,
moderate = 0.025,
severe   = 0.040
)

# Surgery probabilities by severity

surgery_prob <- c(
mild     = 0.15,
moderate = 0.40,
severe   = 0.70
)

# Repeat surgery probabilities

repeat_surgery_prob <- c(
mild     = 0.05,
moderate = 0.15,
severe   = 0.30
)

# Repeat surgery timing (days after initial)

repeat_surgery_days <- c(3, 7, 14, 21)


# Surgery timing (hours from admission)

surgery_timing <- list(
mild     = c(mean = 48, sd = 24),
moderate = c(mean = 24, sd = 12),
severe   = c(mean = 6, sd = 3)
)

# VTE parameters

base_vte_haz   <- 0.0008
vte_acc_rate   <- 6e-6      # Base acceleration
vte_surg_boost <- 2.0       # Multiplier after surgery bleeding resolves

# AC effects

ac_vte_reduc   <- 0.70
ac_bleed_incr  <- 1.50

# Lab thresholds for AC initiation

ac_thresholds <- list(
platelets = 50,
hgb = 8
)

##############################################################################

# HELPER FUNCTIONS

##############################################################################

# Calculate GCS trajectory (improves over time)

get_current_gcs <- function(initial_gcs, hour) {

# GCS improves by 2 points per 24 hours, capped at 15

improvement <- floor(hour / 24) * 2
min(initial_gcs + improvement, 15)
}

# Realistic bleeding hazard pattern

calc_bleed_hazard <- function(hour, severity, surgery_time = NA, surgery_bleed_duration = 24) {
peak_rate <- peak_bleed_by_severity[severity]

# Trauma bleeding: high plateau then drop

if (hour <= bleed_plateau_duration) {
trauma_bleed <- peak_rate
} else {
# Gradual transition to baseline over next 24 hours
transition_period <- 30
time_since_plateau <- hour - bleed_plateau_duration
if (time_since_plateau <= transition_period) {
# Linear decrease from peak to baseline
trauma_bleed <- peak_rate - (peak_rate - bleed_baseline) * (time_since_plateau / transition_period)
} else {
trauma_bleed <- bleed_baseline
}
}

# Add surgical bleeding if applicable

surgical_bleed <- 0
if (![is.na](http://is.na/)(surgery_time) && hour >= surgery_time && hour < (surgery_time + surgery_bleed_duration)) {
# Surgical bleeding peaks immediately then decays
time_since_surgery <- hour - surgery_time
surgical_bleed <- 0.03 * exp(-time_since_surgery / 12)  # 12-hour half-life
}

trauma_bleed + surgical_bleed
}

# VTE hazard with post-surgical acceleration

calc_vte_hazard <- function(hour, surgery_time = NA) {

# Base linear increase

base_vte <- base_vte_haz + vte_acc_rate * hour

# Post-surgical boost after bleeding period (>48h post-op)

if (![is.na](http://is.na/)(surgery_time) && hour > (surgery_time + 48)) {
base_vte <- base_vte * vte_surg_boost
}

base_vte
}

# Generate time-varying labs

generate_labs <- function(hour, cum_bleed, severity) {

# Base values depend on severity

base_plts <- switch(severity, mild = 280, moderate = 250, severe = 220)
base_hgb <- switch(severity, mild = 14, moderate = 13, severe = 12)

# Platelets and Hgb drop with bleeding

## They drop initially after surgery

plts <- rnorm(1, base_plts - cum_bleed * 30, 20) %>% pmax(20) %>% pmin(600)
hgb <- rnorm(1, base_hgb - cum_bleed * 1.5, 0.5) %>% pmax(6) %>% pmin(18)

list(platelets = plts, hgb = hgb)
}

# Schedule surgery for a patient

schedule_surgery <- function(severity) {
if (runif(1) > surgery_prob[severity]) {
return(NA)  # No surgery
}

# Sample surgery time

params <- surgery_timing[[severity]]
surgery_time <- rnorm(1, params["mean"], params["sd"])
max(1, min(surgery_time, max_hours - 48))  # Ensure reasonable bounds
}

##############################################################################

# PATIENT SIMULATION

##############################################################################

simulate_trauma_patient <- function(patient_id) {

# Initialize patient

severity <- sample(severity_levels, 1, prob = severity_probs)
gcs_range <- gcs_ranges[[severity]]
initial_gcs <- sample(gcs_range[1]:gcs_range[2], 1)
surgery_time <- schedule_surgery(severity)

# Patient characteristics

patient <- list(
id = patient_id,
severity = severity,
initial_gcs = initial_gcs,
surgery_time = surgery_time,
has_surgery = ![is.na](http://is.na/)(surgery_time)
)

# Initialize state

state <- list(
alive = TRUE,
on_ac = FALSE,
cum_bleed = 0,
cum_dvt = 0,
cum_pe = 0,
last_bleed_hour = NA,
ac_start_hour = NA,
had_surgery = FALSE,
surgery_count = 0,
next_surgery = NA,
had_dvt = FALSE,
had_pe = FALSE,
death_hour = NA
)

# Simulate hour by hour

hourly_data <- list()

for (hour in hrs) {
if (!state$alive) break

```
# Update GCS
current_gcs <- get_current_gcs(initial_gcs, hour)

# Update surgery status
if (!is.na(surgery_time) && hour >= surgery_time && !state$had_surgery) {
  state$had_surgery <- TRUE
}

# Generate labs
labs <- generate_labs(hour, state$cum_bleed, severity)

# AC management decisions
ac_allowed <- TRUE

# Cannot give AC around surgery (12h before, 24h after)
if (!is.na(surgery_time)) {
  if (hour >= (surgery_time - 12) && hour <= (surgery_time + 24)) {
    ac_allowed <- FALSE
    if (state$on_ac) {
      state$on_ac <- FALSE  # Must stop for surgery
    }
  }
}

# Consider starting AC if allowed and not already on
if (!state$on_ac && ac_allowed && hour <= 168) {  # AC window constraint
  # Safety criteria
  safe_to_start <- labs$platelets > ac_thresholds$platelets &&
                  labs$hgb > ac_thresholds$hgb &&
                  hour >= 1  # Minimum 1 hour wait

  # Additional criteria based on last bleed
  if (!is.na(state$last_bleed_hour)) {
    safe_to_start <- safe_to_start && (hour - state$last_bleed_hour > 24)
  }

  # Start if safe OR if VTE occurred (indicated use)
  if (safe_to_start || (state$had_dvt || state$had_pe)) {
    state$on_ac <- TRUE
    state$ac_start_hour <- hour
  }
}

# Check for repeat surgery
if (state$had_surgery && state$surgery_count < length(repeat_surgery_days)) {
  days_since_surgery <- (hour - surgery_time) / 24
  if (days_since_surgery >= repeat_surgery_days[state$surgery_count + 1]) {
    if (runif(1) < repeat_surgery_prob[severity]) {
      state$surgery_count <- state$surgery_count + 1
      state$next_surgery <- hour + sample(12:36, 1)  # 12-36h from now
    }
  }
}

# Update surgery status for repeat surgery
if (!is.na(state$next_surgery) && hour >= state$next_surgery) {
  surgery_time <- state$next_surgery  # Update for bleeding calculation
  state$next_surgery <- NA
}

# Calculate current hazards
bleed_hz <- calc_bleed_hazard(hour, severity, surgery_time)
vte_hz <- calc_vte_hazard(hour, surgery_time)

# Split VTE into DVT and PE
dvt_hz <- vte_hz * 0.75  # 75% of VTE are DVT
pe_hz <- vte_hz * 0.25   # 25% are PE
if (state$had_dvt) pe_hz <- pe_hz * 2  # Higher PE risk after DVT

# Apply AC effects
if (state$on_ac) {
  bleed_hz <- bleed_hz * ac_bleed_incr
  dvt_hz   <- dvt_hz * (1 - ac_vte_reduc)
  pe_hz    <- pe_hz * (1 - ac_vte_reduc)
}

# Simulate events
new_bleed <- FALSE
new_dvt <- FALSE
new_pe <- FALSE

if (runif(1) < bleed_hz) {
  state$cum_bleed <- state$cum_bleed + 1
  state$last_bleed_hour <- hour
  new_bleed <- TRUE
  # Stop AC for significant bleeding
  if (state$on_ac) {
    state$on_ac <- FALSE
  }
}

if (!state$had_dvt && runif(1) < dvt_hz) {
  state$cum_dvt <- state$cum_dvt + 1
  state$had_dvt <- TRUE
  new_dvt <- TRUE
  # Start AC for DVT if allowed and safe
  if (!state$on_ac && ac_allowed && hour <= 168) {
    state$on_ac <- TRUE
    state$ac_start_hour <- hour
  }
}

if (!state$had_pe && runif(1) < pe_hz) {
  state$cum_pe <- state$cum_pe + 1
  state$had_pe <- TRUE
  new_pe <- TRUE
  # Must start AC for PE if allowed
  if (!state$on_ac && ac_allowed && hour <= 168) {
    state$on_ac <- TRUE
    state$ac_start_hour <- hour
  }
}

# Mortality with proper multipliers
mort_hz <- mort_base[severity]
if (state$cum_bleed > 0) mort_hz <- mort_hz * (mort_bleed_mult ^ state$cum_bleed)
if (state$had_dvt) mort_hz <- mort_hz * mort_dvt_mult
if (state$had_pe) mort_hz <- mort_hz * mort_pe_mult

if (runif(1) < mort_hz) {
  state$alive <- FALSE
  state$death_hour <- hour
}

# Record hourly data
hourly_data[[length(hourly_data) + 1]] <- tibble(
  patient_id = patient_id,
  hour = hour,
  severity = severity,
  gcs = current_gcs,
  on_ac = state$on_ac,
  surgery_scheduled = !is.na(surgery_time) && hour < surgery_time,
  had_surgery = state$had_surgery,
  time_since_surgery = if(state$had_surgery) hour - surgery_time else NA,
  platelets = labs$platelets,
  hgb = labs$hgb,
  bleed_hz = bleed_hz,
  dvt_hz = dvt_hz,
  pe_hz = pe_hz,
  new_bleed = as.numeric(new_bleed),
  new_dvt = as.numeric(new_dvt),
  new_pe = as.numeric(new_pe),
  new_vte = as.numeric(new_dvt || new_pe),
  cum_bleed = state$cum_bleed,
  cum_dvt = state$cum_dvt,
  cum_pe = state$cum_pe,
  had_dvt = state$had_dvt,
  had_pe = state$had_pe,
  alive = state$alive,
  died = as.numeric(!state$alive && hour == state$death_hour)
)

```

}

bind_rows(hourly_data)
}

##############################################################################

# MAIN SIMULATION

##############################################################################

cat("→ Simulating", n_patients, "trauma patients with improved dynamics...\n")

# Run simulation in parallel

library(parallel)
n_cores <- min(detectCores() - 1, 8)
cl <- makeCluster(n_cores)
clusterExport(cl, ls())
clusterEvalQ(cl, {
suppressPackageStartupMessages(library(tidyverse))
})

patient_data <- parLapply(cl, 1:n_patients, simulate_trauma_patient)
stopCluster(cl)

# Combine all patient data

full_data <- bind_rows(patient_data)

cat("→ Generated", nrow(full_data), "patient-hour observations\n")

# Surgery summary

surgery_summary <- full_data %>%
group_by(patient_id) %>%
summarise(
had_surgery = any(had_surgery),
severity = first(severity),
.groups = "drop"
)

cat(sprintf("→ Surgery rates by severity:\n"))
surgery_summary %>%
group_by(severity) %>%
summarise(
pct_surgery = mean(had_surgery) * 100,
.groups = "drop"
) %>%
print()

##############################################################################

# ANALYSIS: Cox Models

##############################################################################

cat("\n→ Fitting time-varying Cox models...\n")

# Prepare survival data

surv_data <- full_data %>%
group_by(patient_id) %>%
mutate(
tstart = hour,
tstop = lead(hour, default = max(hour) + 1)
) %>%
ungroup()

# Death model with time-varying effects of VTE and bleeding

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

surv_data_alive <- surv_data %>% filter(alive)

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

# VISUALIZATION

##############################################################################

cat("\n→ Creating visualizations...\n")

# Color scheme

col_vte   <- "#0066CC"
col_bleed <- "#CC0000"
col_ac    <- "#00AA44"

# Panel 1: Average hazards over time by severity

hazard_summary <- full_data %>%
filter(!on_ac) %>%  # Baseline hazards only
group_by(hour, severity) %>%
summarise(
mean_bleed_hz = mean(bleed_hz),
mean_vte_hz = mean(vte_hz),
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

ac_uptake <- full_data %>%
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

example_patients <- full_data %>%
filter(patient_id %in% c(
sample(unique(patient_id[severity == "mild"]), 1),
sample(unique(patient_id[severity == "moderate"]), 1),
sample(unique(patient_id[severity == "severe"]), 1),
sample(unique(patient_id[had_surgery]), 1)
))

p3 <- ggplot(example_patients, aes(x = hour)) +
geom_line(aes(y = bleed_hz, color = "Bleeding"), size = 1) +
geom_line(aes(y = vte_hz, color = "VTE"), size = 1) +
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

# Panel 4: Cumulative outcomes

cumulative_summary <- full_data %>%
group_by(hour) %>%
summarise(
mean_cum_bleed = mean(cum_bleed),
mean_cum_vte = mean(cum_vte),
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
title = "Improved Trauma Time-Varying Analysis with Death Outcomes",
subtitle = sprintf("n = %d patients; 1-month follow-up; AC window 1-168h", n_patients),
theme = theme(
plot.title = element_text(size = 16, face = "bold"),
plot.subtitle = element_text(size = 12)
)
)

# Save plot

ggsave("trauma_improved_analysis.png",
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

ac_timing <- full_data %>%
filter(on_ac) %>%
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

death_summary <- full_data %>%
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

cat("\n✓ Analysis complete! Plot saved to: trauma_improved_analysis.png\n")
