# For continuous outcome, MAF:0.05-0.11, calculates power for true and misspecified  models using SKAT
# change n and output location for different sample size settings (n=500, 1000, 2000)
library(gdsfmt)
library(SeqArray)
library(STAAR)
library(SKAT)

# --- GLOBAL SETUP ---

GDS_PATH <- "/rsrch6/home/biostatistics/nkui/ukbiobank_gds/ukb_imp_chr1_v3_qc.gds"
n <- 1000
alpha <- 0.05
runs <- 100 
k <- 50

lee<-function(c){
  muq=c[1]
  sigmaq=sqrt(2*c[2])
  s1=c[3]/c[2]^(3/2)
  s2=c[4]/c[2]^2
  if (s1^2>s2)
  {a=1/(s1-sqrt(s1^2-s2))
  delta= (s2*a^4-a^2)/2
  }else
  {a=1/sqrt(s2)
  delta=0}
  l=a^2-2*delta
  mux=l+delta
  sigmax=sqrt(2)*sqrt(l+2*delta)
  return(list(l=l,delta=delta,mux=mux,sigmax=sigmax,muq=muq,sigmaq=sigmaq))
}

# --- Scenario Definitions ---

# Fixed Effect Sets: Indices 3 to 6 of the causal_indices vector
ES_FIXED_SETS <- list(
  set0 = c(0, 0, 0, 0),
  set1 = c(0.1, 0.1, -0.1, -0.1),
  set2 = c(0.2, 0.2, -0.2, -0.2)
)
ES_FIXED_NAMES <- c(0, 0.1, 0.2)


# Define ALL possible values for the simulation grid
es1_values <- c(0.1, 0.2, 0.3, 0.4) # Causal SNP 1
es2_values <- c(-0.1, -0.2, -0.3, -0.4) # Causal SNP 2

# Calculate the total number of scenarios (jobs)
NUM_FIXED_SETS <- length(ES_FIXED_SETS) # = 2
NUM_ES1 <- length(es1_values)           # = 4
NUM_ES2 <- length(es2_values)           # = 4

SCENARIOS_PER_SET <-  NUM_ES1 * NUM_ES2 # = 16
TOTAL_SCENARIOS <- NUM_FIXED_SETS * SCENARIOS_PER_SET # = 32

# Fixed Causal Architecture Definition
indices_multi <- c(1, 2)
indices_additional <- seq(from = 11, to = 41, by = 10)
causal_indices <- c(indices_multi, indices_additional) # [1, 2, 11, 21, 31, 41]

# --- GDS Data Extraction (Remains the same) ---
gdsfile <- seqOpen(GDS_PATH)
sample_all <- seqGetData(gdsfile, "sample.id")
set.seed(42)
sample_sub <- sample(sample_all, n, replace = FALSE)

# Apply Sample Filter
#seqSetFilter(gdsfile, sample.id = sample_sub,variant.id = c(37940:41940))
seqSetFilter(gdsfile, sample.id = sample_sub,variant.id = c(81573:83573))
seqSetFilterCond(gdsfile, maf = 0.05) 
variantid <- seqGetData(gdsfile, "variant.id")
position <- seqGetData(gdsfile, "position")
Geno <- SeqArray::seqGetData(gdsfile, "annotation/format/DS")
alt <- colMeans(Geno) /2
MAF <- pmin(alt,1-alt)
#maf_filter <- which(MAF > 0.17 & MAF <0.24)
maf_filter <- which(MAF > 0.05 & MAF <0.11)
seqSetFilter(gdsfile, sample.id = sample_sub,variant.id = variantid[maf_filter][1:k]) 

# Extract the fixed genotype matrices (Flipping dosage from Ref count to Alt count)
Geno_raw <- SeqArray::seqGetData(gdsfile, "annotation/format/DS")
G <- matrix(2-Geno_raw, nrow = n, ncol = k) # The N x M (1000 x 50) genotype matrix

# The Multi-Allelic Combined Genotype Matrix (G2)
first <- pmin(G[, 1] + G[, 2], 2)
G2 <- cbind(first, G[, 3:k])

# Center the genotype matrices
Gc <- t(t(G) - colMeans(G))
Gc2 <- t(t(G2) - colMeans(G2))

# Close GDS file early
seqClose(gdsfile) 

# ----------------------------------------------------------------------
# SIMULATION FUNCTION FOR A SINGLE SCENARIO
# ----------------------------------------------------------------------

power_sim_C <- function( es1, es2, Gc, Gc2, G, G2, causal_indices, alpha, runs, n, k, es_fixed_set, es_fixed_name) {
  
  # --- SCENARIO SETUP ---
  scenario_name <- paste0( "_ES1_", es1, "_ES2_", es2, "_Fixed_", es_fixed_name)
  cat(sprintf("\n--- Running Scenario: %s ---\n", scenario_name))
  
  # The first two effects are varying (es1, es2), the last four are fixed (es_fixed_set)
  causal_effect <- c(es1, es2, es_fixed_set) 
  beta <- rep(0, k)
  
  beta[causal_indices] <- causal_effect
  
  #---------------------analytical power-------------------------------#
  Var_G1 <- var(Gc[, 1])
  Var_G2 <- var(Gc[, 2])
  mubeta <-  Gc %*% beta
  
  # C. Mis-specified Causal Effect (Combined Model)
  
  
  eps <- 1e-8
  XtX <- t(Gc2) %*% Gc2
  XtX_reg <- XtX + diag(eps, nrow(XtX))
  gamma <- as.numeric(solve(XtX_reg, t(Gc2) %*% mubeta))
  mugamma <- Gc2 %*% gamma
  
  
  
  # --- 1a. Model Separated (G/Gc) ---
  A <- t(Gc) %*% Gc / n
  B <- t(Gc) %*% mubeta %*% t(mubeta) %*% Gc / n^2
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
  
  # --- 1b. Model Combined (G2/Gc2) ---
  C <- t(Gc2) %*% Gc2 / n
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
  powerT2 <- pchisq(q1_c, df = alt_c$l, ncp = alt_c$delta, lower.tail = FALSE)
  
  #---------------------simulated power-------------------------------#
  
  skat_sep_pvals <- numeric(runs)
  skat_comb_pvals <- numeric(runs)
  
  for (i in 1:runs) {
    set.seed(i) 
    
    y_i <- rnorm(n,mubeta,1)
    # --- Separated Model (G) ---
    MAF1 <- colMeans(G) / 2
    w1 <- dbeta(MAF1, 1, 1) # Beta(1,1) weights
    all1 <- as.data.frame(cbind(y_i, G))
    null_mod1 <- SKAT_Null_Model(y_i ~ 1, out_type = "C", data = all1) 
    skat_sep_pvals[i] <- SKAT(G, obj = null_mod1, weights = w1)$p.value
    
    # --- Combined Model (G2) ---
    MAF2 <- colMeans(G2) / 2
    w2 <- dbeta(MAF2, 1, 1) # Beta(1,1) weights
    all2 <- as.data.frame(cbind(y_i, G2))
    null_mod2 <- SKAT_Null_Model(y_i ~ 1, out_type = "C", data = all2)
    skat_comb_pvals[i] <- SKAT(G2, obj = null_mod2, weights = w2)$p.value
  }
  
  powerSeparated <- mean(skat_sep_pvals < alpha)
  powerCombined <- mean(skat_comb_pvals < alpha)
  
  # --- 3. FINAL RESULTS CONSOLIDATION ---
  
  results_df <- data.frame(
    JobIndex = 0, 
    es1 = es1,
    es2 = es2,
    es_other = es_fixed_name,
    powerT = powerT,
    powerT2 = powerT2,
    PowerSeparated = powerSeparated,
    PowerCombined = powerCombined
  )
  
  return(results_df)
}

# ----------------------------------------------------------------------
# SCRIPT EXECUTION BASED ON JOB INDEX (Handles 96 Scenarios)
# ----------------------------------------------------------------------

# 1. Get the job index from command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Error: No job index provided. Use Rscript snpC.R ${LSB_JOBINDEX}", call. = FALSE)
}
job_index <- as.numeric(args[1])

if (job_index < 1 || job_index > TOTAL_SCENARIOS) {
  stop(paste("Error: Job index", job_index, "is out of range (1 to", TOTAL_SCENARIOS, ")"), call. = FALSE)
}

# 2. Map the single job index (1 to 32) directly to all three factors

# Calculate the 0-based index for mapping
idx_rem <- job_index - 1 

# A. Map to the Fixed Set (i_fixed)
# i_fixed determines which of the two fixed effect sets to use (1 or 2)
i_fixed <- floor(idx_rem / SCENARIOS_PER_SET) + 1
current_es_fixed_set <- ES_FIXED_SETS[[i_fixed]]
current_es_fixed_name <- ES_FIXED_NAMES[i_fixed]
# Reduce the index to within the fixed set (0 to 15)
idx_rem <- idx_rem %% SCENARIOS_PER_SET


# B. Map to ES1 (Causal SNP 1)
# i_es1 determines which of the four es1_values to use (1 to 4)
i_es1 <- floor(idx_rem / NUM_ES2) + 1
# Reduce the index to within the es1 level (0 to 3)
idx_rem <- idx_rem %% NUM_ES2
current_es1 <- es1_values[i_es1]


# C. Map to ES2 (Causal SNP 2)
# i_es2 determines which of the four es2_values to use (1 to 4)
i_es2 <- idx_rem + 1
current_es2 <- es2_values[i_es2]

# 4. Execute the single simulation
final_result_df <- power_sim_C(
  current_es1, current_es2,
  Gc, Gc2, G, G2, causal_indices,
  alpha, runs, n, k, current_es_fixed_set, current_es_fixed_name
)

# 5. Update and Save the result
final_result_df$JobIndex <- job_index

output_filename <- paste0("realCwT0.05_1k_", current_es_fixed_name, "_", job_index, ".csv")

output_path <- "/rsrch6/home/biostatistics/nkui/simu/output/"
write.csv(final_result_df, 
          file = paste0(output_path, output_filename), 
          row.names = FALSE)

cat(sprintf("\n--- Scenario %d Complete. Results saved to %s ---\n", job_index, output_filename))
