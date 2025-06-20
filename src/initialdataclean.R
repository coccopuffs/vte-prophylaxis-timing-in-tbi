library(tidyverse)
library(data.table)
library(fastDummies)

#TODO
# - look at TBICEREBRALMONITORHRS and the rest of the NA table to fix it

start <- fread("df.csv")

# MSM Dataset Creation ----------------------------------------------------

# Load necessary libraries
# ── 1 · Load Raw Data ──────────────────────────────────────────────────────────
# Using fread for speed with large datasets
start <- fread("df.csv")

# ── 2 · Initial Cleaning and Filtering ───────────────────────────────────────
# This block combines your initial data prep steps into one pipeline
# for clarity and efficiency.
df_processed <- start %>% 
  # Recode binary variables (2 -> 0) and flip RESPIRATORYASSISTANCE logic
  mutate(
    across(
      c(SEX, WITHDRAWALLST, TBIMIDLINESHIFT, RESPIRATORYASSISTANCE, 
        SUPPLEMENTALOXYGEN, INTERFACILITYTRANSFER, PREHOSPITALCARDIACARREST, 
        ANTIBIOTICTHERAPY),
      ~ if_else(.x == 2, 0, .x)
    ),
    RESPIRATORYASSISTANCE = 1 - RESPIRATORYASSISTANCE # Flips 0 to 1 and 1 to 0
  ) %>%
  ## ── A · recode the prophylaxis type ─────────────────────────────
  mutate(
    VTEPROPHYLAXISTYPE = dplyr::recode(
      VTEPROPHYLAXISTYPE,
      "LWMH"                   = 1L,   # typo “LWMH” fixed here
      "Unfractionated Heparin" = 2L,
      "Xa Inhibitor"           = 3L,
      .default                 = NA_integer_
    )
  ) %>% 
  # Apply cohort filters, now including Xa Inhibitors
  # Using numeric codes: 6=LMWH, 11=Unfractionated Heparin, 8=Xa Inhibitor
  filter(
    SEX != 3,
    TBIMIDLINESHIFT != 3,
    VTEPROPHYLAXISTYPE %in% c(1, 2, 3),
    VTEPROPHYLAXISHRS < 168 # Less than one week
  ) %>%
  # Drop variables that are not needed or have extensive missingness
  select(
    -inc_key, -V1, -EMSNOTIFYHRS, -EMSLEFTHRS, -EMSDepartureHrs, -InpatientHrs,
    -HOSPITALDISCHARGEHRS, -HIGHESTACTIVATION, -TRAUMASURGEONARRIVALHRS,
    -WITHDRAWALLSTHRS, -LOWESTSBP, -ICP_NA, -DEATHINED, -RACE_NA, -TCC_NA,
    -VPO_NA, -PREHOSPITALCARDIACARREST, -GCSQ_NA, -PMGCSQ_NA, -starts_with("HE_")  
  ) %>%
  # Convert all columns to numeric, suppressing warnings for NAs
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.x))))

# ── 3 · Blood Product Unit Conversion ────────────────────────────────────────
# Helper function remains the same
## ── Blood-product helpers ─────────────────────────────────────────────────
## ── Helper: mL ➜ units (base‐R internals) ──────────────────────────────────
ml2unit_fix <- function(x,
                        measure,          # 1 = “already units”, 2 = millilitres
                        default_vol,      # mL per unit if conv is NA
                        conv        = NA_real_,
                        max_units   = 60,
                        min_units   = 0.5,
                        round_dec   = 1) {
  
  # pick the appropriate divisor row-by-row
  ml_per_unit <- ifelse(is.na(conv), default_vol, conv)
  
  # convert to units or leave as-is
  units <- ifelse(
    measure == 1 & !is.na(x),                         x,
    ifelse(measure == 2 & !is.na(x), x / ml_per_unit, NA_real_)
  )
  
  # tiny positive amounts (<0.5 U) count as 1 whole unit
  units <- ifelse(units > 0 & units < min_units, 1, units)
  
  # clamp implausible or negative values
  units[units < 0 | units > max_units] <- NA_real_
  
  # optional rounding
  round(units, round_dec)
}

## ── Apply to each product, tidy clean-up, plausibility gate ───────────────
df_cleaned <- df_processed %>% 
  mutate(
    PRBC_units  = ml2unit_fix(BLOOD4HOURS,          BLOODMEASURE,          300, BLOODCONVERSION),
    PLT_units   = ml2unit_fix(PLATELETS4HOURS,      PLATELETSMEASURE,      250, PLATELETSCONVERSION),
    CRYO_units  = ml2unit_fix(CRYOPRECIPITATE4HOURS,CRYOPRECIPITATEMEASURE,  15, CRYOPRECIPITATECONVERSION),
    FFP_units   = ml2unit_fix(PLASMA4HOURS,         PLASMAMEASURE,         250, PLASMACONVERSION)
  ) %>% 
  
  # drop the raw mL columns
  select(
    -starts_with("BLOOD"), -starts_with("PLATELETS"),
    -starts_with("CRYOPRECIPITATE"), -starts_with("PLASMA"),
    -WHOLEBLOOD4HOURS
  ) %>% 
  
  # replace remaining NAs in *_units with 0
  mutate(across(ends_with("_units"), \(x) tidyr::replace_na(x, 0))) %>%
  # plausibility filter
  filter(
    !(
      (PRBC_units >= 10 & (PLT_units < 1 | FFP_units < 1 | CRYO_units < 1)) |
        PRBC_units  > 60 | PLT_units > 30 | FFP_units > 60 | CRYO_units > 50
    )
  )

setDT(df_cleaned)

## ── 1 · Pre-compute key variables (no copy) ────────────────────────────────
df_cleaned[, `:=`(
  id        = .I,                       # row number → id
  end_time4 = floor(VTEPROPHYLAXISHRS/4)   # store once; avoids recomputing
)]

## ── 2 · Expand to long (4-h grid) with CJ() ────────────────────────────────
# For each id build the sequence 0,4,8,…,end_time
grid <- df_cleaned[, .(time = seq(0, end_time4*4, by = 4)), by = id]

# Merge the patient-level columns onto the time grid
long <- merge(grid, df_cleaned, by = "id", all.x = TRUE)

## ── 3 · time-varying vars (first set) ──────────────────────────
long[, `:=`(
  VTE_prophylaxis  = +(time == end_time4*4),
  GCS_Motor_tv     = fifelse(time < 24, GCSMOTOR, TBIGCSMOTOR),
  On_Ventilator_tv = fifelse(time < 24,
                             +(GCSQ_INTUBATED==1 | RESPIRATORYASSISTANCE==1),
                             +(PMGCSQ_INTUBATED==1 | GCSQ_INTUBATED==1 |
                                 RESPIRATORYASSISTANCE==1)),
  Midline_Shift_tv = fifelse(time < 24, 0, TBIMIDLINESHIFT),
  Patient_Location_tv = fcase(
    time < EDDISCHARGEHRS,            "ED",
    EDDISCHARGEDISPOSITION==1,        "Floor",
    EDDISCHARGEDISPOSITION==2,        "Observation",
    EDDISCHARGEDISPOSITION==3,        "Step-down",
    EDDISCHARGEDISPOSITION==7,        "Operating_Room",
    EDDISCHARGEDISPOSITION==8,        "ICU",
    default = "Discharged_from_ED")
), by = id]

## ── 3b · add *_units_tv in a separate call ────────────────────
unitCols <- grep("_units$", names(long), value = TRUE)

long[, (paste0(unitCols, "_tv")) :=
       lapply(.SD, \(x) fifelse(time < 4, 0, x)),
     .SDcols = unitCols, by = id]

## ── 4 · Time-dependent interventions (compare once, then recycle) ──────────
long[, Angiography_tv :=
       fifelse( is.na(ANGIOGRAPHYHRS), 0L,        # never got it → 0
                as.integer(time >= ANGIOGRAPHYHRS) ),
     by = id]

long[, Surgery_for_Hemorrhage_tv :=
       fifelse( is.na(HMRRHGCTRLSURGHRS), 0L,
                as.integer(time >= HMRRHGCTRLSURGHRS) ),
     by = id]

long[, Antibiotic_Therapy_tv :=
       fifelse( is.na(ANTIBIOTICTHERAPYHRS), 0L,
                as.integer(time >= ANTIBIOTICTHERAPYHRS) ),
     by = id]

## ── 5 · Last-observation-carried-forward fill (fast rolling join) ─────────
# choose the columns that need LOCF
fill_vars <- c("GCS_Motor_tv", "On_Ventilator_tv",
               "Midline_Shift_tv", "Patient_Location_tv",
               grep("_tv$", names(long), value = TRUE))

setorderv(long, c("id","time"))          # ensure sorted
# encode as factor → integer
long[, Patient_Location_code := as.integer(factor(Patient_Location_tv))]

# numeric LOCF
long[, Patient_Location_code := nafill(Patient_Location_code, type = "locf"),
     by = id]

# decode back to original labels
labs <- levels(factor(long$Patient_Location_tv))
long[, Patient_Location_tv := labs[Patient_Location_code]]

# optional: drop helper column
long[, Patient_Location_code := NULL]

# drop helper column
long[, end_time4 := NULL]


# ── 5 · Final Cleanup and Export for Weighting ───────────────────────────────
df_weighting <- long %>%
  # One-hot encode the key categorical variables
  fastDummies::dummy_cols(
    select_columns = c("VTEPROPHYLAXISTYPE", "Patient_Location_tv", "ETHNICITY", 
                       "TBIPUPILLARYRESPONSE", 
                       "HMRRHGCTRLSURGTYPE", "ANGIOGRAPHY"),
    remove_selected_columns = TRUE,
    ignore_na = TRUE
  ) %>%
  # Remove raw variables that are now encoded in the time-varying format or are outcomes
  select(
    -VTEPROPHYLAXISHRS, -EDDISCHARGEHRS, -EDDISCHARGEDISPOSITION,
    -TBIGCSMOTOR, -TBIMIDLINESHIFT, -GCSQ_INTUBATED, -RESPIRATORYASSISTANCE,
    -PMGCSQ_INTUBATED, -starts_with("ANGIOGRAPHY", ignore.case = FALSE), 
    -starts_with("HMRRHGCTRLSURG", ignore.case = FALSE),
    -starts_with("ANTIBIOTIC", ignore.case = FALSE), -ends_with("_units"),
    -FINALDISCHARGEHRS, -WITHDRAWALLST # Remove outcomes
  )


na_table <- df_weighting %>% 
  summarise(across(everything(),
                   ~ sum(is.na(.x)))) %>% 
  tidyr::pivot_longer(everything(),
                      names_to  = "variable",
                      values_to = "n_NA") %>% 
  arrange(desc(n_NA))

print(na_table, n = Inf)   # show all rows


# Save the weighting dataset
fwrite(df_weighting, "df_weighting.csv")

# Clustering Dataset Creation--------------------------------------------------------------

df_km <- df_processed %>% 
  mutate(
    PRBC_units  = ml2unit_fix(BLOOD4HOURS,          BLOODMEASURE,          300, BLOODCONVERSION),
    PLT_units   = ml2unit_fix(PLATELETS4HOURS,      PLATELETSMEASURE,      250, PLATELETSCONVERSION),
    CRYO_units  = ml2unit_fix(CRYOPRECIPITATE4HOURS,CRYOPRECIPITATEMEASURE,  15, CRYOPRECIPITATECONVERSION),
    FFP_units   = ml2unit_fix(PLASMA4HOURS,         PLASMAMEASURE,         250, PLASMACONVERSION)
  ) %>% 
  select(
    -starts_with("BLOOD"), -starts_with("PLATELETS"),
    -starts_with("CRYOPRECIPITATE"), -starts_with("PLASMA"),
    -WHOLEBLOOD4HOURS
  ) %>% 
  mutate(across(ends_with("_units"), \(x) replace_na(x, 0)))

## ── 2 · OPTIONAL plausibility filter (same rules) ─────────────────────────
df_km <- df_km %>% 
  filter(
    !(
      (PRBC_units >= 10 & (PLT_units < 1 | FFP_units < 1 | CRYO_units < 1)) |
        PRBC_units  > 60 | PLT_units > 30 | FFP_units > 60 | CRYO_units > 50
    )
  )

## ── 3 · ONE-HOT ENCODE key categoricals for k-means  ──────────────────────
df_km <- df_km %>% 
  fastDummies::dummy_cols(
    select_columns = c("VTEPROPHYLAXISTYPE", "ETHNICITY",
                       "TBIPUPILLARYRESPONSE", "HMRRHGCTRLSURGTYPE"),
    remove_selected_columns = TRUE,
    ignore_na = TRUE
  )

## ── 4 · (Optional) replace remaining NA with column median 
df_km <- df_km %>% 
  mutate(across(where(is.numeric),
                \(x) replace_na(x, median(x, na.rm = TRUE))))


fwrite(as.data.table(df_km), "df_kmeans.csv")
