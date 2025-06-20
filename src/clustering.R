library(data.table)
library(sparcl)                # KMeansSparseCluster
library(doParallel)           
library(foreach)
library(missRanger)

# TODO
#- for some reason imputation is not detecting values to be impute
# I think its because i did median assignments in the dataset creation
#- do i still need the implementation of getting optimal parameters on a sample
# and then running this pararmeters on the whole dataset


# 1. Setting up parallel computing ----------
tot_cores <- parallel::detectCores()          
setDTthreads(tot_cores)                       
Sys.setenv(OMP_NUM_THREADS = 8)               

## load CSV 
df <- fread("df_kmeans.csv")                 

# 2. imputation ----------

if (anyNA(df)) {                      
  message("Imputing missing values with missRanger …")
  
  df <- missRanger(
    df,
    num.trees   = 500,
    maxiter     = 1,
    pmm.k       = 5,
    num.threads = parallel::detectCores(),
    verbose     = 1,
    seed        = 42,
    respect.unordered.factors = "partition"
  )
  
} else {
  message("No NA detected – skipping imputation step.")
}

# 3. preparation for kmeans ----------

X_full <- df %>%                              
  # keep only numeric columns
  select(where(is.numeric)) %>%               
  # drop zero-variance columns before scaling
  select(where(~ var(.x, na.rm = TRUE) > 1e-12)) %>% 
  # scaling
  mutate(across(everything(), scale)) %>%     
  # convert to matrix for sparse K-means pkg
  as.matrix()                                 

cat("Final matrix dims (full):",
    nrow(X_full), "×", ncol(X_full), "\n\n")

#  4 · permutation-gap search on 20 k subsample for now -------
set.seed(42)
sub_idx <- sample.int(nrow(X_full), 20000)
X_sub   <- X_full[sub_idx, ]

lambda_grid <- seq(1.2, 7, 0.4)      # sparsity penalties
k_grid      <- 2:9                   # clusters to test
n_perms     <- 10                    # perms per λ

best_lambda <- NA; best_k <- NA; best_gap <- -Inf

for (k in k_grid) {
  perm <- KMeansSparseCluster.permute(
    X_sub,
    K        = k,
    wbounds  = lambda_grid,
    nperms   = n_perms)
  
  λ_idx <- which.max(perm$gap)
  gap   <- perm$gap[λ_idx]
  cat("k", k,
      "→ λ", lambda_grid[λ_idx],
      "gap", round(gap, 3), "\n")
  
  
  # Finding the best combination of k and lambda
  if (gap > best_gap) {
    best_gap    <- gap
    best_k      <- k
    best_lambda <- lambda_grid[λ_idx]
  }
}

cat("\nChosen λ =", best_lambda,
    "| k =", best_k,
    "| gap =", round(best_gap, 3), "\n\n")

# 5 · final sparse K-means on ALL patients -----------
set.seed(42)
fit <- KMeansSparseCluster(X_full,
                           K       = best_k,
                           wbounds = best_lambda,
                           nstart  = 20)

clusters <- fit[[1]]$Cs              # patient labels
weights  <- fit[[1]]$ws              # feature weights

# 6 · save  + supplementary figures
fwrite(data.table(id = df$id, cluster = clusters),
       "kmeans_sparse_clusters.csv")

fwrite(data.table(variable = colnames(X_full),
                  weight   = weights)[weight > 0][order(-weight)],
       "kmeans_sparse_weights.csv")