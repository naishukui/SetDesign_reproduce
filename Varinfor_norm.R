# this code was based on the original code Varinfo_gds.R from STAAR pipeline tutorial
# and added normalization algorithm for REF and ALT.
##########################################################################
#           Input
##########################################################################

### DB split information
file_DBsplit <- "/path_to_the_file/FAVORdatabase_chrsplit.csv"
### Targeted GDS
dir_geno <- "/path_to_the_GDS_file/"
gds_file_name_1 <- "chr"
gds_file_name_2 <- ".gds"

### output
output_path <- "/rsrch6/home/biostatistics/nkui/ukbiobank/output/"

chr <- as.numeric(commandArgs(TRUE)[1])

###########################################################################
#           Main Function
###########################################################################

### make directory
system(paste0("mkdir ",output_path,"chr",chr))

### R package
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)

### chromosome number
## read info
DB_info <- read.csv(file_DBsplit,header=TRUE)
DB_info <- DB_info[DB_info$Chr==chr,]

## open GDS
gds.path <- paste0(dir_geno,gds_file_name_1,chr,gds_file_name_2)
genofile <- seqOpen(gds.path)

CHR <- as.numeric(seqGetData(genofile, "chromosome"))
positionold <- as.integer(seqGetData(genofile, "position"))
REFold <- as.character(seqGetData(genofile, "$ref"))
ALTold <- as.character(seqGetData(genofile, "$alt"))
position<-positionold
REF<-REFold
ALT<-ALTold

# compare REFold and ALTold
for (i in seq_along(REFold)) {

  xs <- strsplit(REFold[i],"")[[1]]
  ys <- strsplit(ALTold[i],"")[[1]]

  deleted <- 0
  # Check for common characters from the end
  while ( tail(xs, 1) == tail(ys, 1) && length(xs) > 1 && length(ys) > 1 ) {
    # Delete the last element from both vectors
    xs <- head(xs, -1)
    ys <- head(ys, -1)
  }
  # Check for common characters from the beginning
  while(xs[1]==ys[1] & length(ys)>1 & length(xs)>1) {
    deleted <- deleted + 1
    xs=xs[-1]
    ys=ys[-1]
  }

  position[i]<-positionold[i]+deleted
  REF[i]<-paste0(xs, collapse = "")
  ALT[i]<-paste0(ys, collapse = "")
}

VarInfo_genome <- paste0(CHR,"-",position,"-",REF,"-",ALT)

seqClose(genofile)

## Generate VarInfo
for(kk in 1:dim(DB_info)[1])
{
  print(kk)
  VarInfo <- VarInfo_genome[(position>=DB_info$Start_Pos[kk])&(position<=DB_info$End_Pos[kk])]
  VarInfo <- data.frame(VarInfo)
  write.csv(VarInfo,paste0(output_path,"chr",chr,"/VarInfo_old_chr",chr,"_",kk,".csv"),quote=FALSE,row.names = FALSE)
  system(paste0("rm ",output_path,"chr",chr,"/VarInfo_old_chr",chr,"_",kk,".csv"))
}

system(paste0("rm ",output_path,"chr",chr,"/VarInfo_old*"))
