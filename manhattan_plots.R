# this code generates manhattan plots under different STAAR settings

library(qqman)
library(dplyr)
# library(STAAR)
# library(STAARpipeline)
library(ggplot2)
library(calibrate)


manhattan2<-function (x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c("gray10",
                                                                               "gray60"), chrlabs = NULL, suggestiveline = -log10(1e-05),
                      genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE,
                      annotatePval = NULL, annotateTop = TRUE, ...)
{
  CHR = BP = P = index = NULL
  if (!(chr %in% names(x)))
    stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(x)))
    stop(paste("Column", bp, "not found!"))
  if (!(p %in% names(x)))
    stop(paste("Column", p, "not found!"))
  if (!(snp %in% names(x)))
    warning(paste("No SNP column found. OK unless you're trying to highlight."))
  if (!is.numeric(x[[chr]]))
    stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!is.numeric(x[[bp]]))
    stop(paste(bp, "column should be numeric."))
  if (!is.numeric(x[[p]]))
    stop(paste(p, "column should be numeric."))
  if (!is.null(x[[snp]]))
    d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]],
                   pos = NA, index = NA, SNP = x[[snp]], stringsAsFactors = FALSE)
  else d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]],
                      pos = NA, index = NA)
  d <- d[order(d$CHR, d$BP), ]
  if (logp) {
    d$logp <- -log10(d$P)
  }
  else {
    d$logp <- d$P
  }
  d$index = rep.int(seq_along(unique(d$CHR)), times = tapply(d$SNP,
                                                             d$CHR, length))
  nchr = length(unique(d$CHR))
  if (nchr == 1) {
    d$pos = d$BP
    xlabel = paste("Chromosome", unique(d$CHR), "position")
  }
  else {
    lastbase = 0
    ticks = NULL
    for (i in unique(d$index)) {
      if (i == 1) {
        d[d$index == i, ]$pos = d[d$index == i, ]$BP
      }
      else {
        lastbase = lastbase + max(d[d$index == (i - 1),
                                    "BP"])
        d[d$index == i, "BP"] = d[d$index == i, "BP"] -
          min(d[d$index == i, "BP"]) + 1
        d[d$index == i, "pos"] = d[d$index == i, "BP"] +
          lastbase
      }
    }
    ticks <- tapply(d$pos, d$index, quantile, probs = 0.5)
    xlabel = "Chromosome"
    labs <- unique(d$CHR)
  }
  xmax = ceiling(max(d$pos) * 1.03)
  xmin = floor(max(d$pos) * -0.03)
  def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i",
                   las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(0,
                                                                     ceiling(max(d$logp))), xlab = xlabel, ylab = expression(-log[10](italic(p))))
  dotargs <- list(...)
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in%
                                            names(dotargs)]))
  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
      if (length(chrlabs) == length(labs)) {
        labs <- chrlabs
      }
      else {
        warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
      }
    }
    else {
      warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
    }
  }
  if (nchr == 1) {
    axis(1, ...)
  }
  else {
    axis(1, at = ticks, labels = labs, ...)
  }
  col = rep_len(col, max(d$index))
  if (nchr == 1) {
    with(d, points(pos, logp, pch = 20, col = col[1], ...))
  }
  else {
    icol = 1
    for (i in unique(d$index)) {
      points(d[d$index == i, "pos"], d[d$index == i, "logp"],
             col = col[icol], pch = 20, ...)
      icol = icol + 1
    }
  }
  if (suggestiveline)
    abline(h = suggestiveline, col = "blue")
  if (genomewideline)
    abline(h = genomewideline, col = "red")
  if (!is.null(highlight)) {
    if (any(!(highlight %in% d$SNP)))
      warning("You're trying to highlight SNPs that don't exist in your results.")
    d.highlight = d[which(d$SNP %in% highlight), ]
    with(d.highlight, points(pos, logp, col = "green3", pch = 20,
                             ...))
  }
  if (!is.null(annotatePval)) {
    if (logp) {
      topHits = subset(d, P <= annotatePval)
    }
    else topHits = subset(d, P >= annotatePval)
    par(xpd = TRUE)
    if (annotateTop == FALSE) {
      if (logp) {
        with(subset(d, P <= annotatePval), textxy(pos,
                                                  -log10(P), offset = 0.625, labs = topHits$SNP,
                                                  cex = 0.9), ...)
      }
      else with(subset(d, P >= annotatePval), textxy(pos,
                                                     P, offset = 0.625, labs = topHits$SNP, cex = 0.9),
                ...)
    }
    else {
      topHits <- topHits[order(topHits$P), ]
      topSNPs <- NULL
      for (i in unique(topHits$CHR)) {
        chrSNPs <- topHits[topHits$CHR == i, ]
        topSNPs <- rbind(topSNPs, chrSNPs[1, ])
      }
      if (logp) {
        textxy(topSNPs$pos, -log10(topSNPs$P), offset = 0.625,
               labs = topSNPs$SNP, cex = 0.9, ...)
      }
      else textxy(topSNPs$pos, topSNPs$P, offset = 0.625,
                  labs = topSNPs$SNP, cex = 0.9, ...)
    }
  }
  par(xpd = FALSE)
}

#genes <- genes_info
genes<-get(load("/Users/nkui/Library/CloudStorage/OneDrive-InsideMDAnderson(2)/STAAR/output/combined/gene.Rdata"))
#final setting
setting1<-get(load("/Users/nkui/Library/CloudStorage/OneDrive-InsideMDAnderson(2)/STAAR/output/output6-10/coding_combined6.Rdata"))
setting2<-get(load("/Users/nkui/Library/CloudStorage/OneDrive-InsideMDAnderson(2)/STAAR/output/output6-10/coding_combined7.Rdata"))
setting3<-get(load("/Users/nkui/Library/CloudStorage/OneDrive-InsideMDAnderson(2)/STAAR/output/output6-10/coding_combined8.Rdata"))
setting4<-get(load("/Users/nkui/Library/CloudStorage/OneDrive-InsideMDAnderson(2)/STAAR/output/output6-10/coding_combined9.Rdata"))
setting5<-get(load("/Users/nkui/Library/CloudStorage/OneDrive-InsideMDAnderson(2)/STAAR/output/output6-10/coding_combined10.Rdata"))
#original STAAR pipeline
setting6<-get(load("/Users/nkui/Library/CloudStorage/OneDrive-InsideMDAnderson(2)/STAAR/output/output6-10/coding_combined11.Rdata"))


colnames(genes) <- c("gene", "chr", "start_position", "end_position")
custom_colors <- c("#74c2e7","#dea3e5")

#-------------------------------------------manhattan plot for setting 1&6-----------------------------------------------#
prepare_data <- function(setting) {

  # Merge the current setting with gene data using "chr" and "gene" columns
  combo <- merge(setting, genes, by = c("chr", "gene"))

  # Ensure chromosome is numeric
  combo$chr <- as.numeric(as.character(combo$chr))

  # Extract values (adjust column indices if needed)
  pskat         <- as.numeric(combo[, 8])    # p-value column (assumed at col 8)
  gene_names    <- as.character(combo[, 2])  # gene names (assumed at col 2)
  start_position<- as.numeric(combo[, 9])    # start position (assumed at col 9)

  # Build a data frame for Manhattan plot
  allGeneP <- data.frame(chr = combo$chr,
                         pskat = pskat,
                         start_position = start_position,
                         gene = gene_names)

  # Order by chromosome and add a counter to simulate base-pair position
  bychr <- allGeneP[order(allGeneP$chr), ]
  df <- bychr %>%
    group_by(chr) %>%
    mutate(counter = row_number()) %>%
    as.data.frame()

  return(as.data.frame(df))
}
df1 <- prepare_data(setting1)
df2 <- prepare_data(setting2)
df3 <- prepare_data(setting3)
df4 <- prepare_data(setting4)
df5 <- prepare_data(setting5)
df6 <- prepare_data(setting6)


par(mfrow = c(2, 1))
#Manhattan plot showing associations of gene-centric analysis for pancreatic cancer versus −log10(P) of SKAT
#（n=7435）#
par(mar=c(0, 5, 3, 2))
manhattan2(df1,
          chr = "chr",
          bp = "counter",
          snp = "gene",
          p = "pskat",
          col = custom_colors,  # Use custom colors
          #suggestiveline = FALSE,  # Optional: significance threshold
          genomewideline = -log10(5e-8),  # Optional: genome-wide significance threshold
          annotatePval = 0.001,
          ylim = c(0, 6),
          las=2,
          main = "Setting 1 vs. Setting 6")  # Optional: add title
par(mar=c(0.5, 5, 3, 2))

manhattan2(df2,
          chr = "chr",
          bp = "counter",
          snp = "gene",
          p = "pskat",
          col = custom_colors,
          #suggestiveline = -log10(1e-5),
          genomewideline = -log10(5e-8),
          annotatePval = 0.001,
          ylim = c(6, 0),
          xlab="",xaxt="n")

outlier_idx <- which(-log10(combo2$SKAT1_1) > 10)
outlier_x <- which(df2$gene=="ULK1")
outlier_y <- -log10(combo2$SKAT1_1)[outlier_idx]
outlier_snp <- "ULK1"
points(outlier_x, rep(5.5, length(outlier_x)), pch = 17, col = "red", cex = 1)
text(outlier_x, rep(5.5, length(outlier_x)), labels = paste0(outlier_snp, " (~", round(outlier_y, 2), ")"), pos = 4, col = "red", cex = 1)

#--------------------------------------------setting 1&6-------------------------------------------#
par(mfrow = c(4, 1))
par(mar=c(0, 5, 1, 3))
#B L T R
manhattan2(df1,
           chr = "chr",
           bp = "counter",
           snp = "gene",
           p = "pskat",
           col = custom_colors,  # Use custom colors
           #suggestiveline = FALSE,  # Optional: significance threshold
           genomewideline = -log10(5e-8),  # Optional: genome-wide significance threshold
           annotatePval = 0.001,
           ylim = c(0, 6),
           las=2,
           main = " ")  # Optional: add title
mtext("Setting 1", side = 4, line = 1.2, las = 0, srt = 270, cex = 1)


par(mar=c(0, 5, 2, 3))

manhattan2(df6,
           chr = "chr",
           bp = "counter",
           snp = "gene",
           p = "pskat",
           col = custom_colors,
           #suggestiveline = -log10(1e-5),
           genomewideline = -log10(5e-8),
           annotatePval = 0.001,
           ylim = c(6, 0),
           xlab="",xaxt="n")
mtext("Setting 6", side = 4, line = 1.2, las = 0, srt = 270, cex = 1)

outlier_idx <- which(-log10(combo2$SKAT1_1) > 10)
outlier_x <- which(df2$gene=="ULK1")
outlier_y <- -log10(combo2$SKAT1_1)[outlier_idx]
outlier_snp <- "ULK1"
points(outlier_x, rep(5.5, length(outlier_x)), pch = 17, col = "red", cex = 1)
text(outlier_x, rep(5.5, length(outlier_x)), labels = paste0(outlier_snp, " (~", round(outlier_y, 2), ")"), pos = 4, col = "red", cex = 1)



#-------------------------------------------2&3-------------------#
par(mar=c(0, 5, 1, 3))

manhattan2(df2,
           chr = "chr",
           bp = "counter",
           snp = "gene",
           p = "pskat",
           col = custom_colors,  # Use custom colors
           #suggestiveline = FALSE,  # Optional: significance threshold
           genomewideline = -log10(5e-8),  # Optional: genome-wide significance threshold
           annotatePval = 0.001,
           ylim = c(0, 6),
           las=2,
           main = " ")  # Optional: add title
mtext("Setting 2", side = 4, line = 1.2, las = 0, srt = 270, cex = 1)

par(mar=c(0.5, 5, 2, 3))

manhattan2(df3,
           chr = "chr",
           bp = "counter",
           snp = "gene",
           p = "pskat",
           col = custom_colors,
           #suggestiveline = -log10(1e-5),
           genomewideline = -log10(5e-8),
           annotatePval = 0.001,
           ylim = c(6, 0),
           xlab="",xaxt="n")
mtext("Setting 3", side = 4, line = 1.2, las = 0, srt = 270, cex = 1)




#--------------------------------------------setting 4&5-------------------------------------------#
par(mfrow = c(2, 1))

par(mar=c(0, 5, 3, 3))
manhattan2(df4,
           chr = "chr",
           bp = "counter",
           snp = "gene",
           p = "pskat",
           col = custom_colors,  # Use custom colors
           #suggestiveline = FALSE,  # Optional: significance threshold
           genomewideline = -log10(5e-8),  # Optional: genome-wide significance threshold
           annotatePval = 0.001,
           ylim = c(0, 6),
           las=2,
           main = " ")  # Optional: add title
mtext("Setting 4", side = 4, line = 1.2, las = 0, srt = 270, cex = 1)

par(mar=c(0.5, 5, 3, 3))

manhattan2(df5,
           chr = "chr",
           bp = "counter",
           snp = "gene",
           p = "pskat",
           col = custom_colors,
           #suggestiveline = -log10(1e-5),
           genomewideline = -log10(5e-8),
           annotatePval = 0.001,
           ylim = c(6, 0),
           xlab="",xaxt="n")
mtext("Setting 5", side = 4, line = 1.2, las = 0, srt = 270, cex = 1)


library(tidyverse)
library(ggfastman)
data("gwas_data")
gwas_data %>%
  mutate( gr= "Study 1") %>%
  mutate(pvalue=-log10(pvalue)) %>%
  # rbind a second study with pvalues with other sign.
  bind_rows(., mutate(., gr= "Study 2",
                      pvalue = -pvalue)) %>%
  # plot the points
  fast_manhattan(., build = "hg18", speed = "fast",log10p = F, dodge_x = T,pointsize = 2.1, pixels = c(1000,500)) +
  # add significance line
  geom_hline(data= . %>% group_by(gr) %>% slice(1), aes(yintercept = ifelse(pvalue>0, -log10(5e-08),log10(5e-08))),color ="deeppink") +
  facet_wrap(~gr, scales = "free_y",ncol = 1,strip.position = "right")+
  scale_y_continuous(expression(-log[10](italic(p))),breaks= seq(-90,80,10), labels = abs(seq(-90,80,10)), expand = c(0.01, 0))

