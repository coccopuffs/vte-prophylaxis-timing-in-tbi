library(tidyverse)
library(data.table)

start <- fread("df.csv")

df <- start %>% 
  ## recode 2 → 0 for listed binary variables
  mutate(across(
    c(SEX, WITHDRAWALLST, TBIMIDLINESHIFT,
      RESPIRATORYASSISTANCE, SUPPLEMENTALOXYGEN,
      INTERFACILITYTRANSFER, PREHOSPITALCARDIACARREST,
      DEATHINED, BLOODMEASURE, ANTIBIOTICTHERAPY),
    ~ ifelse(.x == 2, 0, .x)
  )) %>%
  ## recode respiratory assistance since its backwards
  mutate(RESPIRATORYASSISTANCE = if_else(RESPIRATORYASSISTANCE == 1, 0, 1)) %>%
  ## apply cohort filters -
  filter(
    SEX != 3,
    TBIMIDLINESHIFT != 3,
    VTEPROPHYLAXISTYPE %in% c("LWMH", "Unfractionated Heparin"),
    VTEPROPHYLAXISHRS < 168
  ) %>% 
  ## set reference level for prophylaxis type
  mutate(
    VTEPROPHYLAXISTYPE = relevel(factor(VTEPROPHYLAXISTYPE),
                                 ref = "Unfractionated Heparin")
  ) %>%
  ## remove variables that are extensively missing
  select(-EMSNOTIFYHRS, -EMSLEFTHRS, -EMSDepartureHrs, -InpatientHrs,
         -HOSPITALDISCHARGEHRS, -HIGHESTACTIVATION, -TRAUMASURGEONARRIVALHRS,
         -inc_key, -V1, -WITHDRAWALLSTHRS, -LOWESTSBP, -ICP_NA, -DEATHINED,
         -RACE_NA, -TCC_NA, -VPO_NA, -PREHOSPITALCARDIACARREST, -GCSQ_NA,
         -PMGCSQ_NA, -V1, -inc_key) %>%
  ## Recoding relevant variable with NAs --> 0
  ## - Need to add blood variables here?
  mutate(across(
    c(
      "TOTALVENTDAYS", "TOTALICULOS",
      "HE_UnplannedIntubation", "HE_Stroke.CVA", "HE_SevereSepsis",
      "HE_DVT", "HE_UnplannedAdmissiontoICU", "HE_UnplannedVisittoOR",
      "Hx_Pregnancy"
    ),
    ~ tidyr::replace_na(.x, 0)
  )) %>%
  ## Making everything numeric
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.x)))) %>%
  ## Remove Hospital Events before weighting
  select(-starts_with("HE")) %>%
  ## Only keeping time variables before VTEPROPHYLAXISHRS
  mutate(
    # keep value only if event happened before VTE prophylaxis; otherwise 0
    ANGIOGRAPHY =
      if_else(!is.na(ANGIOGRAPHYHRS) & ANGIOGRAPHYHRS < VTEPROPHYLAXISHRS,
              ANGIOGRAPHY, 0),
    
    HMRRHGCTRLSURG =
      if_else(!is.na(HMRRHGCTRLSURGHRS) & HMRRHGCTRLSURGHRS < VTEPROPHYLAXISHRS,
              HMRRHGCTRLSURGTYPE, 0),
    
    ANTIBIOTICTHERAPY =
      if_else(!is.na(ANTIBIOTICTHERAPYHRS) & ANTIBIOTICTHERAPYHRS < VTEPROPHYLAXISHRS,
              ANTIBIOTICTHERAPY, 0)
  ) %>% 
  
  # convert any remaining NA in the three columns to 0
  mutate(across(c(ANGIOGRAPHY, HMRRHGCTRLSURGTYPE, ANTIBIOTICTHERAPY),
                ~ replace_na(.x, 0))) %>% 
  # anticoagulation in the ER
  mutate(
    ED_VTE = if_else(
      !is.na(EDDISCHARGEHRS) & !is.na(VTEPROPHYLAXISHRS) &
        EDDISCHARGEHRS > VTEPROPHYLAXISHRS,
      1L,              
      0L              
    )
  ) %>%
  # drop the now-redundant *_HRS columns
  select(-ANTIBIOTICTHERAPYHRS, -ANGIOGRAPHYHRS, -HMRRHGCTRLSURGHRS)
  



# Converting all transfusion products to units given
## Lookup table 
vol_par <- tibble::tribble(
  ~prod,   ~unit_vol, ~hours_var,        ~measure_var,       ~conv_var,
  "PRBC",  300,       "BLOOD4HOURS",     "BLOODMEASURE",     "BLOODCONVERSION",
  "PLT",   250,       "PLATELETS4HOURS", "PLATELETSMEASURE", "PLATELETSCONVERSION",
  "CRYO",   15,       "CRYOPRECIPITATE4HOURS",  "CRYOPRECIPITATEMEASURE","CRYOPRECIPITATECONVERSION",
  "FFP",   250,       "PLASMA4HOURS",    "PLASMAMEASURE",    "PLASMACONVERSION"
  #,"WB",   500,       "WHOLEBLOOD4HOURS","WHOLEBLOODMEASURE","WHOLEBLOODCONVERSION"
)

# helper (unchanged)
ml2unit_fix <- function(x, measure, default_vol,
                        conv        = NA_real_,
                        max_units   = 60,
                        min_u_thres = 0.5) {
  u <- dplyr::case_when(
    measure == 1 & x >= min_u_thres               ~ x,
    measure == 1 & x <  min_u_thres               ~ x / default_vol,
    measure == 2 & !is.na(conv)                   ~ x / conv,
    measure == 2 &  is.na(conv) & x >= default_vol~ x / default_vol,
    measure == 2 &  is.na(conv) & x <  default_vol~ x,
    TRUE                                          ~ NA_real_
  )
  ifelse(u > max_units | u < 0, NA_real_, u)
}

# main pipeline
df <- df %>% 
  mutate(
    PRBC_units  = ml2unit_fix(BLOOD4HOURS,          BLOODMEASURE,          300, BLOODCONVERSION),
    PLT_units   = ml2unit_fix(PLATELETS4HOURS,      PLATELETSMEASURE,      250, PLATELETSCONVERSION),
    CRYO_units  = ml2unit_fix(CRYOPRECIPITATE4HOURS,CRYOPRECIPITATEMEASURE,15,  CRYOPRECIPITATECONVERSION),
    FFP_units   = ml2unit_fix(PLASMA4HOURS,         PLASMAMEASURE,         250, PLASMACONVERSION)
    # ,WB_units = ml2unit_fix(WHOLEBLOOD4HOURS, WHOLEBLOODMEASURE, 500, WHOLEBLOODCONVERSION)
  ) %>% 
  select(
    -BLOOD4HOURS,      -BLOODMEASURE,      -BLOODCONVERSION,
    -PLATELETS4HOURS,  -PLATELETSMEASURE,  -PLATELETSCONVERSION,
    -CRYOPRECIPITATE4HOURS, -CRYOPRECIPITATEMEASURE, -CRYOPRECIPITATECONVERSION,
    -PLASMA4HOURS,     -PLASMAMEASURE,     -PLASMACONVERSION,
    -WHOLEBLOOD4HOURS #, -WHOLEBLOODMEASURE, -WHOLEBLOODCONVERSION
  ) %>% 
  mutate(across(ends_with("_units"), ~ replace_na(.x, 0)))

# flag rows that still look impossible / unbalanced
bad_tbl <- df %>% 
  mutate(bad_row =
           (PRBC_units >= 10 & (PLT_units < 1 | FFP_units < 1 | CRYO_units < 1)) |
           PRBC_units  > 60 | PLT_units  > 30 | FFP_units > 60 | CRYO_units > 50 |
           PRBC_units  < 0  | PLT_units  < 0  | FFP_units < 0  | CRYO_units < 0) %>% 
  filter(bad_row) %>% 
  select(-bad_row)

# drop those rows from the working data
df <- df %>% 
  filter(
    !( (PRBC_units >= 10 & (PLT_units < 1 | FFP_units < 1 | CRYO_units < 1)) |
         PRBC_units  > 60 | PLT_units  > 30 | FFP_units > 60 | CRYO_units > 50 |
         PRBC_units  < 0  | PLT_units  < 0  | FFP_units < 0  | CRYO_units < 0 )
  )


# Making sure all timed interventions or events are before VTEPROPHYLAXISHRS
icp_cols <- grep("^ICP", names(df), value = TRUE)

df <- df %>% 
  ## Zero-out ICP values taken *after* VTE prophylaxis started
  mutate(
    across(
      all_of(icp_cols),
      ~ if_else(
        !is.na(TBICEREBRALMONITORHRS) &
          !is.na(VTEPROPHYLAXISHRS)      &
          TBICEREBRALMONITORHRS >= VTEPROPHYLAXISHRS,
        0,        # after prophylaxis  → 0
        .x        # otherwise keep
      )
    )
  ) %>% 
  mutate(
    ICP_Unknown = if_else(
      !is.na(TBICEREBRALMONITORHRS) &                       # monitor time recorded
        if_all(all_of(icp_cols), ~ is.na(.x) | .x == 0),      # every ICP value NA or 0
      1L,                                                   # flag as 1
      0L                                                    # else 0
    )
  ) %>% 
  select(-TBICEREBRALMONITORHRS) %>%
  mutate(ICP_Unknown = replace_na(ICP_Unknown, 0))

# Adding time dependent variables, making sure not to add variables after vteprophylaxishrs
df <- df %>% 
  # creating ICU Variable
  mutate(
    ICU_Before_VTE = ifelse(
      !is.na(EDDISCHARGEDISPOSITION) & EDDISCHARGEDISPOSITION == 8 &
        !is.na(EDDISCHARGEHRS) & !is.na(VTEPROPHYLAXISHRS) &
        EDDISCHARGEHRS < VTEPROPHYLAXISHRS,
      1,
      0
    )
  ) %>%
  # creating ventilator variable
  mutate(
    Ventilator_Before_VTE = ifelse(
      !is.na(RESPIRATORYASSISTANCE) & RESPIRATORYASSISTANCE == 1 & # Assuming 1 means assisted based on previous recode logic
        !is.na(TOTALVENTDAYS) & TOTALVENTDAYS >= 1, # If any ventilator days are recorded
      1,
      0
    )
  ) %>%
  # CT scan shows shift in brain 24 hours after admission
  mutate(
    TBIMIDLINESHIFT = if_else(
      !is.na(VTEPROPHYLAXISHRS) & VTEPROPHYLAXISHRS < 19,
      0,
      TBIMIDLINESHIFT
    )
  )

df <- df %>%
  select(-TOTALVENTDAYS, -TOTALICULOS, -EDDISCHARGEHRS)

write.csv(df, "df_clean.csv")

# Removing remaining variables after vteprophylaxishrs
df <- df %>%
  select(-FINALDISCHARGEHRS, -WITHDRAWALLST)

write.csv(df, "df_weighting.csv")
