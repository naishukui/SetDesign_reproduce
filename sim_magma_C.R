# For continuous outcome, prepare  genotype and phenotype data for true and misspecified  models in order to run MAGMA
library(gdsfmt)
library(SeqArray)
library(STAAR)
library(SKAT)
library(rootSolve)
library(genio)


GDS_PATH <- "/rsrch6/home/biostatistics/nkui/ukbiobank_gds/ukb_imp_chr1_v3_qc.gds"
n <- 1000
alpha <- 0.05
runs <- 100 
k <- 50

# Define ALL possible values for the simulation grid
ES_FIXED_SETS <- list(
  set0 = c(0,0,0,0),
  set1 = c(0.1, 0.1, -0.1, -0.1),
  set2 = c(0.2, 0.2, -0.2, -0.2)
)
ES_FIXED_NAMES <- c(0, 0.1, 0.2)
es1_values <- c(0.1, 0.2, 0.3, 0.4) # Causal SNP 1
es2_values <- c(-0.1, -0.2, -0.3, -0.4) # Causal SNP 2

# Calculate the total number of scenarios (jobs)
NUM_FIXED_SETS <- length(ES_FIXED_SETS) # = 3
NUM_ES1 <- length(es1_values)           # = 4
NUM_ES2 <- length(es2_values)           # = 4

SCENARIOS_PER_SET <-  NUM_ES1 * NUM_ES2 # = 16
TOTAL_SCENARIOS <- NUM_FIXED_SETS * SCENARIOS_PER_SET # = 48

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
position[2] <- position[1]+1

CHR <- seqGetData(gdsfile,"chromosome")
REF <- seqGetData(gdsfile,"$ref")
ALT <- seqGetData(gdsfile,"$alt")
sampleid <- sample_sub
rsid <- seqGetData(gdsfile,"annotation/id")
rsid[2] <- "rs11111111"


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

seqClose(gdsfile) 


#---------------------------------geino-------------------------------#
genotypes  <- t(G)

rownames(genotypes) <- rsid
colnames(genotypes) <- sample_sub
genotypes2 <- t(G2)
rownames(genotypes2) <- rsid[-2]
colnames(genotypes2) <- sample_sub

bim_data <- data.frame(
  chr = CHR,
  id = rsid,
  pos = position, # Base-pair position
  alt = ALT, # Allele 1 (alternate)
  ref = REF, # Allele 2 (reference)
  posg = rep(0, k) # Genetic distance (0 if unknown)
  
)

bim_data2<- data.frame(
  chr = CHR[-2],
  id = rsid[-2],
  pos = position[-2], # Base-pair position
  alt = ALT[-2], # Allele 1 (alternate)
  ref = REF[-2], # Allele 2 (reference)
  posg = rep(0, k-1) # Genetic distance (0 if unknown)
  
)



# ----------------------------------------------------------------------
# SIMULATION FUNCTION FOR A SINGLE SCENARIO
# ----------------------------------------------------------------------

# 1. Get the job index from command line arguments
args <- commandArgs(trailingOnly = TRUE)
job_index <- as.numeric(args[1])
idx_rem <- job_index - 1 
i_fixed <- floor(idx_rem / SCENARIOS_PER_SET) + 1
current_fixed <- ES_FIXED_SETS[[i_fixed]]
current_fixed_name <- ES_FIXED_NAMES[i_fixed]
idx_rem <- idx_rem %% SCENARIOS_PER_SET
i_es1 <- floor(idx_rem / NUM_ES2) + 1
idx_rem <- idx_rem %% NUM_ES2
current_es1 <- es1_values[i_es1]
i_es2 <- idx_rem + 1
current_es2 <- es2_values[i_es2]



pheno_<- function( es1, es2, Gc, Gc2, G, G2, causal_indices, runs, n, k, es_fixed, es_fixed_name) {
  # --- SCENARIO SETUP ---
  scenario_name <- paste0("_fixed_", es_fixed_name,"_es1_", es1, "_es2_", es2)
  cat(sprintf("\n--- Running Scenario: %s ---\n", scenario_name))
  causal_effect <- c(es1, es2, es_fixed) 
  beta <- rep(0, k)
  beta[causal_indices] <- causal_effect
  #eta <- beta0 + Gc %*% beta
  #pii <- exp(eta) / (exp(eta) + 1)
  mubeta <-  Gc %*% beta
  
  output_path <- paste0("/rsrch6/home/biostatistics/nkui/magma/maf0.05C/","job",job_index)
  system(paste0("mkdir ",output_path))
  for (i in 1:runs) {
    set.seed(i) # Set seed for each run to ensure different y_i
    #y_i <- rbinom(n, 1, pii) +1 # Simulate binary y (Noise is here!)
    y_i <- rnorm(n,mubeta,1)
    
    fam_data <- data.frame(
      fam = sample_sub, # Family ID
      id = sample_sub, # Individual ID
      pat = rep(0,n), # Paternal ID (0 = unknown)
      mat = rep(0,n), # Maternal ID (0 = unknown)
      sex = rep(0,n), # Sex (1=M, 2=F, 0=Unknown)
      pheno = y_i# Phenotype (1=control, 2=case, -9=missing)
    )
    
    write_plink(
      file = paste0(output_path,"/true_", i),
      X = genotypes,
      bim = bim_data,
      fam = fam_data
    )
    
    write_plink(
      file = paste0(output_path,"/mis_", i),
      X = genotypes2,
      bim = bim_data2,
      fam = fam_data
    )
    
    
  }
  
}

pheno_(
  es1         = current_es1,
  es2         = current_es2,
  Gc          = Gc,
  Gc2         = Gc2,
  G           = G,
  G2          = G2,
  causal_indices = causal_indices,
  runs        = runs,
  n           = n,
  k           = k,
  es_fixed    = current_fixed,
  es_fixed_name = current_fixed_name
)


