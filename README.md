**Title:** Optimizing the Timing of Venous Thromboembolism (VTE) Prophylaxis After Traumatic Brain Injury

**Authors:** Mahesh Challapalli, BA; Andrew Warburton, MD; Daniel Katz, MD

**Keywords:** National Trauma Data Bank (NTDB); MIMIC-IV; traumatic brain injury; VTE prophylaxis; marginal structural model; phenotyping

<h3>Methodology (Concise Overview):</h3>

**1. Data Preparation**
- National Trauma Data Standard cohort, missing values addressed with multivariate chained‐equations imputation


**2. IPTW (Inverse Probability of Treatment Weighting)**
- Stabilized weights estimated via generalized boosting machine (GBM) models using all pre-exposure covariates.


**3. Primary Outcome Modeling: Marginal structural models fitted with three separate techniques**     
- Bayesian change-point regression (in progress)
- Bayesian penalized spline regression (in progress))
- Random survival forests


**4. Clinical Phenotype Discovery**
- Unsupervised clustering (kmeans) to derive clinically meaningful TBI phenotypes


**5. Phenotype-Specific Analyses**
- Steps 2–3 repeated within each cluster


**6. External Validation (in progress)**
- Primary rsf model tested on MIMIC-IV patients and graded on outcome
