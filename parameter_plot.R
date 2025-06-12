# supplementary figure 1
# compare parameters in true model beta1 and beta2 for multi-allelic variants to parameters in misspecified model
# alpha1 for uncorrelated genotypes and alpha2 for correlated genotypes
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(rootSolve)

# Define the vectors
beta1 <- c(-0.9, -0.75, -0.5, -0.25, -0.1,
           -0.9, -0.75, -0.5, -0.25, -0.1,
           -0.9, -0.75, -0.5, -0.25, -0.1,
           -0.9, -0.75, -0.5, -0.25, -0.1,
           -0.9, -0.75, -0.5, -0.25, -0.1)
beta2 <- c(rep(0.9, 5), rep(0.75, 5), rep(0.5, 5), rep(0.25, 5), rep(0.1, 5))
alpha1<-get_param_mis(cor=F,p=0.01,beta1list=beta1,beta2list=beta2)$alpha1
alpha2<-get_param_mis(cor=T,p=0.05,rho=0.15,beta1list=beta1,beta2list=beta2)$alpha1


# Create a data frame
data <- data.frame(
  Index = seq_along(beta1),
  beta1 = beta1,
  beta2 = beta2,
  alpha1 = alpha1,
  alpha2 = alpha2,

  Group = rep(c(0.9, 0.75, 0.5, 0.25, 0.1), each = 5)
)

# Transform the data to a long format
data_long <- pivot_longer(
  data,
  cols = c(beta1, beta2, alpha1,alpha2),
  names_to = "Variable",
  values_to = "Value"
)

# Define Greek labels for the legend
variable_labels <- c(
  beta1 = expression(beta[49]),
  beta2 = expression(beta[50]),
  alpha1 = expression(alpha[49]),
  alpha2 = expression(alpha[49] * "'")
  )

ggplot(data_long, aes(x = Index, y = Value, color = Variable)) +
  geom_line() +
  geom_point() +
  # Facet in 5 subgroups, with Greek label for beta2
  facet_wrap(
    ~ Group,
    scales = "free_x",
    labeller = label_bquote(beta[50] == .(Group))
  ) +
  # Use default discrete color scale, remove legend title, but keep Greek labels
  scale_color_discrete(
    name = NULL,          # Remove the legend title
    labels = variable_labels
  ) +
  # Main title with Greek subscripts
  labs(
    title = expression("Parameters in True Models and Misspecified Models"),
    x = NULL,       # remove the x-axis label
    y = "Effect Size"
  ) +
  # theme_classic() to remove grid lines; then tweak strips and x-axis
  theme_classic(base_size = 13) +
  theme(
    # Remove the box around panel titles
    strip.background = element_blank(),
    # Remove x-axis text and ticks
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

