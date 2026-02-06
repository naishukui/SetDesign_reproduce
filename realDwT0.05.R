# # For binary outcome, MAF:0.05-0.11, calculates power for true and misspecified  models using SKAT
# change n and output location for different sample size settings (n=500, 1000, 2000)

library(glmnet)
library(gdsfmt)
library(SeqArray)
library(STAAR)
library(SKAT)
library(rootSolve)

source("/rsrch6/home/biostatistics/nkui/simu/get_param_mis.R")
source("/rsrch6/home/biostatistics/nkui/simu/lee.R")


GDS_PATH <- "/rsrch6/home/biostatistics/nkui/ukbiobank_gds/ukb_imp_chr1_v3_qc.gds"
n <- 500
alpha <- 0.05
runs <- 100 
k <- 50

# Define ALL possible values for the simulation grid
ES_FIXED_SETS <- list(
  set0 = c(0,0,0,0),
  set1 = c(0.1, 0.1, -0.1, -0.1),
  set2 = c(0.25, 0.25, -0.25, -0.25)
)
ES_FIXED_NAMES <- c(0, 0.1, 0.25)

beta0_values <- c(-1, 0, 1)
es1_values <- c(0.1, 0.25, 0.5, 0.75)
es2_values <- c(-0.1, -0.25, -0.5, -0.75) 

# Calculate the total number of scenarios (jobs)
NUM_FIXED <- length(ES_FIXED_NAMES)
NUM_B0 <- length(beta0_values)
NUM_es1 <- length(es1_values)
NUM_es2 <- length(es2_values)
TOTAL_SCENARIOS <- NUM_B0 * NUM_es1 * NUM_es2 * NUM_FIXED
# TOTAL_SCENARIOS = 3 * 4 * 4 * 3 = 48*3

# Fixed Causal Architecture Definition
indices_multi <- c(1, 2)
indices_additional <- seq(from = 11, to = 41, by = 10)
causal_indices <- c(indices_multi, indices_additional) # [1, 2, 11, 21, 31, 41]

# Open and filter the GDS file for the fixed genotype matrix G
gdsfile <- seqOpen(GDS_PATH)
sample_all <- seqGetData(gdsfile, "sample.id")
set.seed(42)
sample_sub <- sample(sample_all, n, replace = FALSE)

# Apply Sample Filter
#seqSetFilter(gdsfile, sample.id = sample_sub,variant.id = c(37940:41940))
seqSetFilter(gdsfile, sample.id = sample_sub,variant.id = c(81573:83573))
#seqSetFilterCond(gdsfile, maf = 0.15) 
seqSetFilterCond(gdsfile, maf = 0.05) 
variantid <- seqGetData(gdsfile, "variant.id")
position <- seqGetData(gdsfile, "position")
Geno <- SeqArray::seqGetData(gdsfile, "annotation/format/DS")
alt <- colMeans(Geno) /2
MAF <- pmin(alt,1-alt)
#maf_filter <- which(MAF > 0.17 & MAF <0.24)
maf_filter <- which(MAF > 0.05 & MAF <0.11)
seqSetFilter(gdsfile, sample.id = sample_sub,variant.id = variantid[maf_filter][1:k]) 

position <- seqGetData(gdsfile,"position")
# Extract the fixed genotype matrices
Geno_raw <- SeqArray::seqGetData(gdsfile, "annotation/format/DS")
G <- matrix(2-Geno_raw, nrow = n, ncol = k) # The N x M (1000 x 50) genotype matrix

# The Multi-Allelic Combined Genotype Matrix (G2)
first <- pmin(G[, 1] + G[, 2], 2)
ck <-cbind(G[,c(1,2)],first)
G2 <- cbind(first, G[, 3:k])

# Center the genotype matrices
Gc <- t(t(G) - colMeans(G))
Gc2 <- t(t(G2) - colMeans(G2))
colnames(Gc2) <- c("V12", paste0("V", 3:k))

# Close GDS file early to save resources, as all necessary data is extracted
seqClose(gdsfile) 

# ----------------------------------------------------------------------
# SIMULATION FUNCTION FOR A SINGLE SCENARIO
# ----------------------------------------------------------------------


# ----------------------------------------------------------------------
# SCRIPT EXECUTION BASED ON JOB INDEX
# ----------------------------------------------------------------------

# 1. Get the job index from command line arguments
args <- commandArgs(trailingOnly = TRUE)
job_index <- as.numeric(args[1])

idx_rem <- job_index - 1  # zero-based

i_fixed <- floor(idx_rem / (NUM_B0 * NUM_es1 * NUM_es2)) + 1
idx_rem <- idx_rem %% (NUM_B0 * NUM_es1 * NUM_es2)

i_B0 <- floor(idx_rem / (NUM_es1 * NUM_es2)) + 1
idx_rem <- idx_rem %% (NUM_es1 * NUM_es2)

i_es1 <- floor(idx_rem / NUM_es2) + 1
idx_rem <- idx_rem %% NUM_es2

i_es2 <- idx_rem + 1

# --- Extract current parameter set ---
current_fixed_name <- ES_FIXED_NAMES[i_fixed]
current_es_fixed <- ES_FIXED_SETS[[i_fixed]]
current_beta0 <- beta0_values[i_B0]
current_es1 <- es1_values[i_es1]
current_es2 <- es2_values[i_es2]

cat("Job index:", job_index, "\n",
    "ES_fixed:", current_fixed_name, "\n",
    "beta0:", current_beta0, "\n",
    "es1:", current_es1, "\n",
    "es2:", current_es2, "\n")

power_sim_D <- function(beta0, es1, es2, Gc, Gc2, G, G2, causal_indices, alpha, runs, n, k, es_fixed, es_fixed_name) {
  
  # --- SCENARIO SETUP ---
  scenario_name <- paste0("beta0_", beta0, "_es1_", es1, "_es2_", es2)
  cat(sprintf("\n--- Running Scenario: %s ---\n", scenario_name))
  
  # A. Define Scenario-Specific Effect Sizes (beta)
  causal_effect <- c(es1, es2, es_fixed) # Combine fixed and varying effects
  beta <- rep(0, k)
  gamma <- rep(0, k-1)
  
  beta[causal_indices] <- causal_effect
  
  # B. Linear Predictors (Fixed for the entire scenario)
  eta <- beta0 + Gc %*% beta
  pii <- exp(eta) / (exp(eta) + 1)
  pi0 = exp(beta0) / (1 + exp(beta0))
  mubeta <- pii - pi0
  v <- as.vector(pii * (1 - pii))
  V <- diag(v)
  
  A <- t(Gc) %*% V %*% Gc/n
  B <- t(Gc) %*% mubeta %*% t(mubeta) %*% Gc/n^2
  A2 <- A %*% A; A3 <- A2 %*% A
  
  c1_s <- c(sum(diag(A)) * n, sum(A * t(A)) * n^2, sum(A2 * t(A)) * n^3, sum(A2 * t(A2)) * n^4)
  c2_s <- c(sum(diag(B)) * n^2, 2 * sum(A * t(B)) * n^3, 3 * sum(A2 * t(B)) * n^4, 4 * sum(A3 * t(B)) * n^5)
  c_a_s <- c1_s + c2_s
  
  null_s <- lee(c1_s)
  q0_s <- qchisq(1 - alpha, df = null_s$l, ncp = null_s$delta)
  qc_s <- (q0_s - null_s$mux) * null_s$sigmaq / null_s$sigmax + null_s$muq
  alt_s <- lee(c_a_s)
  q1_s <- (qc_s - alt_s$muq) * alt_s$sigmax / alt_s$sigmaq + alt_s$mux
  powerT <- pchisq(q1_s, df = alt_s$l, ncp = alt_s$delta, lower.tail = FALSE)
  
  # Theoretical power calculation ends here.
  
  # --- EMPIRICAL POWER (Monte Carlo Simulation) with robust coef extraction ---
  skat_sep_pvals <- numeric(runs)
  skat_comb_pvals <- numeric(runs)
  powerT2_sim <- numeric(runs)
  gamma_matrix <- matrix(NA, nrow = runs, ncol = k) 
  
  for (i in 1:runs) {
    set.seed(i)
    y_i <- rbinom(n, 1, pii)
    
    cv_fit <- cv.glmnet(Gc2, y_i, family = "binomial", alpha = 0)
    coef_mat <- coef(cv_fit, s = "lambda.min")
    gamma_matrix[i,]<- coef_mat[1:k]     # The 49 average slope parameters
    gamma0 <- gamma_matrix[i,1]
    gamma <- gamma_matrix[i,2:k]
    eta2 <- gamma0 + Gc2 %*% gamma   
    pii2 <- exp(eta2) / (exp(eta2) + 1)   
    pi02 = exp(gamma0) / (1 + exp(gamma0))   
    mugamma <- pii2 - pi02   
    v2 <- as.vector(pii2 * (1 - pii2))   
    V2 <- diag(v2) 
    
    C <- t(Gc2) %*% V2 %*% Gc2 / n
    D <- t(Gc2) %*% mugamma %*% t(mugamma) %*% Gc2 / n^2
    C2 <- C %*% C; C3 <- C2 %*% C
    
    c1_c <- c(sum(diag(C)) * n, sum(C * t(C)) * n^2, sum(C2 * t(C)) * n^3, sum(C2 * t(C2)) * n^4)
    c2_c <- c(sum(diag(D)) * n^2, 2 * sum(C * t(D)) * n^3, 3 * sum(C2 * t(D)) * n^4, 4 * sum(C3 * t(D)) * n^5)
    c_a_c <- c1_c + c2_c
    
    null_c <- lee(c1_c)
    q0_c <- qchisq(1 - alpha, df = null_c$l, ncp = null_c$delta)
    qc_c <- (q0_c - null_c$mux) * null_c$sigmaq / null_c$sigmax + null_c$muq
    alt_c <- lee(c_a_c)
    q1_c <- (qc_c - alt_c$muq) * alt_c$sigmax / alt_c$sigmaq + alt_c$mux
    powerT2_sim[i] <- pchisq(q1_c, df = alt_c$l, ncp = alt_c$delta, lower.tail = FALSE)
    
    # --- Separated model p-values (unchanged) ---
    MAF1 <- colMeans(G) / 2
    w1 <- dbeta(MAF1, 1, 1)
    all1 <- as.data.frame(cbind(y_i, G))
    null_mod1 <- SKAT_Null_Model(y_i ~ 1, out_type = "D", data = all1)
    skat_sep_pvals[i] <- SKAT(G, obj = null_mod1, weights = w1)$p.value
    
    # --- Combined model p-values (unchanged) ---
    MAF2 <- colMeans(G2) / 2
    w2 <- dbeta(MAF2, 1, 1)
    all2 <- as.data.frame(cbind(y_i, G2))
    null_mod2 <- SKAT_Null_Model(y_i ~ 1, out_type = "D", data = all2)
    skat_comb_pvals[i] <- SKAT(G2, obj = null_mod2, weights = w2)$p.value
  } 
  
  # Empirical powers
  powerSeparated <- mean(skat_sep_pvals < alpha)
  powerCombined <- mean(skat_comb_pvals < alpha)
  powerT2 <- mean(powerT2_sim)
  
  # --- 3. FINAL RESULTS CONSOLIDATION ---
  
  results_df <- data.frame(
    JobIndex = 0, # Placeholder, will be updated below
    Beta0 = beta0,
    es1 = es1,
    es2 = es2,
    es_other = es_fixed_name,
    PowerT = powerT,
    PowerT2 = powerT2,
    PowerSeparated = powerSeparated,
    PowerCombined = powerCombined
  )
  
  return(results_df)
}

# --- Run simulation ---

final_result_df <- power_sim_D(
  current_beta0, current_es1, current_es2,
  Gc, Gc2, G, G2, causal_indices,
  alpha, runs, n, k, current_es_fixed, current_fixed_name
)

# 4. Update and Save the result
final_result_df$JobIndex <- job_index

output_filename <- paste0("realDwT0.05_500_", job_index, ".csv")
write.csv(final_result_df, file = paste0("/rsrch6/home/biostatistics/nkui/simu/output/",output_filename), row.names = FALSE)

cat(sprintf("\n--- Scenario %d Complete. Results saved to %s ---\n", job_index, output_filename))
