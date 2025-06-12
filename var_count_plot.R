# this code generates plots comparing number of variants under different STAAR settings
library(qqman)
library(dplyr)
library(ggplot2)


#genes <- genes_info
genes<-get(load("/path/to/file/gene.Rdata"))
setting1<-get(load("/path/to/file/setting1.Rdata"))
setting2<-get(load("/path/to/file/setting2.Rdata"))
setting3<-get(load("/path/to/file/setting3.Rdata"))
setting4<-get(load("/path/to/file/setting4.Rdata"))
setting5<-get(load("/path/to/file/setting5.Rdata"))
setting6<-get(load("/path/to/file/setting6.Rdata"))


colnames(genes) <- c("gene", "chr", "start_position", "end_position")


#--------------------------------------num of variants qq plot-----------------------------------------#
df <- setting1 %>% select(gene, SNV_s1 = SNV) %>%
  inner_join(setting2 %>% select(gene, SNV_s2 = SNV), by = "gene") %>%
  inner_join(setting3 %>% select(gene, SNV_s3 = SNV), by = "gene") %>%
  inner_join(setting4 %>% select(gene, SNV_s4 = SNV), by = "gene") %>%
  inner_join(setting5 %>% select(gene, SNV_s5 = SNV), by = "gene") %>%
  inner_join(setting6 %>% select(gene, SNV_s6 = SNV), by = "gene")

# Define a list of data pairs (each element is a list with x and y)
qq_data <- list(
  list(x = df$SNV_s2, y = df$SNV_s1),
  list(x = df$SNV_s6, y = df$SNV_s5)
)

# Define labels for x and y axes for each chart
qq_labels <- list(
  c("SNV", "SNV+MNV"),
  c("Single Composite", "Multiple")
)

# Define titles for each chart
titles <- c(
  "Mutation Type",
  "Multi-alleles"
)

# Set up a 2x2 plotting layout
par(mfrow = c(1, 2))
#par(mfrow = c(1, 1))

# Loop through each pair to create a QQ plot
for (i in seq_along(qq_data)) {
  x <- log(as.numeric(qq_data[[i]]$x))
  y <- log(as.numeric(qq_data[[i]]$y))
  x <- as.numeric(qq_data[[i]]$x)
  y <- as.numeric(qq_data[[i]]$y)

  # Create the QQ plot comparing the quantiles of x and y
  plot(x, y,
       main = titles[i],
       xlab = qq_labels[[i]][1],
       ylab = qq_labels[[i]][2],
       pch = 19,
       cex.lab = 1.5,     # Increase axis label size
       cex.axis = 1.2,    # Increase axis tick label size
       cex.main = 1.5  ,# Use solid circle markers
       col = "#74c2e7")   # Set point color

  # Add a 45-degree reference line (slope 1, intercept 0)
  abline(0, 1, col = "#dea3e5", lwd = 2, lty = 2)
}
