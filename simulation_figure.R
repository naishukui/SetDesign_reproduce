# generate figure for simulated power and analytical power under four scenarios

library(SKAT)
library(bindata)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(KATSP)
library(ggpubr)
library(grid)
library(SetDesign)

beta1list1<-c(-0.9,-0.75,-0.5,-0.25,-0.1,
             -0.9,-0.75,-0.5,-0.25,-0.1,
             -0.9,-0.75,-0.5,-0.25,-0.1,
             -0.9,-0.75,-0.5,-0.25,-0.1,
             -0.9,-0.75,-0.5,-0.25,-0.1)
beta2list1<-c( rep(0.9,5),rep(0.75,5),rep(0.5,5),rep(0.25,5),rep(0.1,5))

beta1list11<-c(-0.5,-0.4,-0.3,-0.2,-0.1,
               -0.5,-0.4,-0.3,-0.2,-0.1,
               -0.5,-0.4,-0.3,-0.2,-0.1,
               -0.5,-0.4,-0.3,-0.2,-0.1,
               -0.5,-0.4,-0.3,-0.2,-0.1)
beta2list11<-c( rep(0.5,5),rep(0.4,5),rep(0.3,5),rep(0.2,5),rep(0.1,5))

theme_set(theme_cowplot())
mytheme = list(
  theme_classic()+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          legend.position="none", )
)

##---------------------------------Continuous outcome, 50 uncorrelated SNPs, MAF=0.01, 0.1<=beta<=0.9
powerDerived<-powerC_derive(k=50,n=2000,alpha=0.05,p=0.01,list1=beta1list1,list2=beta2list1)
result<-get(load("/path_to_the_file/snpC.Rdata"))
powerE2 <- sapply(result, function(x) x$powerE2)
powerCombined <- sapply(result, function(x) x$powerCombined)
powerSeparated <- sapply(result, function(x) x$powerSeparated)
df<-as.data.frame(cbind(as.vector(beta1list1),as.vector(beta2list1),powerSeparated,powerCombined,powerDerived,powerE2))
colnames(df)<-c("beta1","beta2","power1","power2","power3","power4")

power1<-cbind(df[,1:3],rep("Sim,True",25))
power2<-cbind(df[,c(1,2,4)],rep("Sim,Mis",25))
power3<-cbind(df[,c(1,2,5)],rep("Analytic,True",25))
power4<-cbind(df[,c(1,2,6)],rep("Analytic,Mis",25))

colnames(power1)[3:4]<-c("power","group")
colnames(power2)[3:4]<-c("power","group")
colnames(power3)[3:4]<-c("power","group")
colnames(power4)[3:4]<-c("power","group")

dfnew<-rbind(power1,power2,power3,power4)
beta2<-as.vector(c(0.1,0.25,0.5,0.75))
dfsub<-list()
line_types <- c("solid", "solid", "dashed", "dashed")

for (i in 1:4){
  dfsub[[i]] <- dfnew[dfnew$beta2==beta2[i], ]
  assign(paste0("p", i), ggplot(dfsub[[i]], aes(x = beta1, y = power, group=group)) +
           geom_point(aes(color=group,shape=group)) +
           geom_line(aes(color=group,linetype=group)) +
           scale_linetype_manual(values = line_types) +  # Set linetype

           labs(title = bquote(beta[50] == .(beta2[i])),
                x = bquote("\u03B2"[49]),
                y = "Power",
                color = "",
                shape = "",
                linetype = "") +
           mytheme)
  mytheme
}

label <- textGrob("(a)", rot = 0, gp = gpar(fontsize = 12))
fig1 <- ggarrange(label,p1, p2, p3, p4, ncol = 5, nrow = 1,widths=c(0.1,0.2,0.2,0.2,0.2))
print(fig1)

##---------------------------------Continuous outcome, 50 correlated SNPs,MAF=0.05, rho=0.15, 0.1<=beta<=0.5

powerDerived<-powerC_deriveCor(k=50,n=2000,alpha=0.05,p=0.05, rho=0.15,list1=beta1list11,list2=beta2list11)
result<-get(load("/path_to_the_file/snpCcor.Rdata"))
powerE2 <- sapply(result, function(x) x$powerE2)
powerCombined <- sapply(result, function(x) x$powerCombined)
powerSeparated <- sapply(result, function(x) x$powerSeparated)
df<-as.data.frame(cbind(as.vector(beta1list11),as.vector(beta2list11),powerSeparated,powerCombined,powerDerived,powerE2))
colnames(df)<-c("beta1","beta2","power1","power2","power3","power4")

power1<-cbind(df[,1:3],rep("Sim,True",25))
power2<-cbind(df[,c(1,2,4)],rep("Sim,Mis",25))
power3<-cbind(df[,c(1,2,5)],rep("Analytic,True",25))
power4<-cbind(df[,c(1,2,6)],rep("Analytic,Mis",25))

colnames(power1)[3:4]<-c("power","group")
colnames(power2)[3:4]<-c("power","group")
colnames(power3)[3:4]<-c("power","group")
colnames(power4)[3:4]<-c("power","group")

dfnew<-rbind(power1,power2,power3,power4)
beta2<-as.vector(c(0.1,0.2,0.3,0.4))
dfsub<-list()
line_types <- c("solid", "solid", "dashed", "dashed")

for (i in 1:4){
  dfsub[[i]] <- dfnew[dfnew$beta2==beta2[i], ]
  assign(paste0("p", i), ggplot(dfsub[[i]], aes(x = beta1, y = power, group=group)) +
           geom_point(aes(color=group,shape=group)) +
           geom_line(aes(color=group,linetype=group)) +
           scale_linetype_manual(values = line_types) +  # Set linetype
           labs(title = bquote(beta[50] == .(beta2[i])),
                x = bquote("\u03B2"[49]),
                y = "Power",
                color = "",
                shape = "",
                linetype = "") +
           mytheme)
  mytheme


}
label <- textGrob("(b)", rot = 0, gp = gpar(fontsize = 12))
fig2 <- ggarrange(label,p1, p2, p3, p4, ncol = 5, nrow = 1,widths=c(0.1,0.2,0.2,0.2,0.2))
fig2

##---------------------------------Binary outcome, 50 uncorrelated SNPs, MAF=0.01, 0.1<=beta<=0.9

powerDerived<-powerD_derive(kk=50,n=2000,alpha=0.05,p=0.01,list1=beta1list1,list2=beta2list1)
result<-get(load("/path_to_the_output_file/snpD.Rdata"))
powerE2 <- sapply(result, function(x) x$powerE2)
powerCombined <- sapply(result, function(x) x$powerCombined)
powerSeparated <- sapply(result, function(x) x$powerSeparated)
df<-as.data.frame(cbind(as.vector(beta1list1),as.vector(beta2list1),powerSeparated,powerCombined,powerDerived,powerE2))
colnames(df)<-c("beta1","beta2","power1","power2","power3","power4")

power1<-cbind(df[,1:3],rep("Sim,True",25))
power2<-cbind(df[,c(1,2,4)],rep("Sim,Mis",25))
power3<-cbind(df[,c(1,2,5)],rep("Analytic,True",25))
power4<-cbind(df[,c(1,2,6)],rep("Analytic,Mis",25))

colnames(power1)[3:4]<-c("power","group")
colnames(power2)[3:4]<-c("power","group")
colnames(power3)[3:4]<-c("power","group")
colnames(power4)[3:4]<-c("power","group")

dfnew<-rbind(power1,power2,power3,power4)
beta2<-as.vector(c(0.1,0.25,0.5,0.75))
dfsub<-list()
line_types <- c("solid", "solid", "dashed", "dashed")

for (i in 1:4){
  dfsub[[i]] <- dfnew[dfnew$beta2==beta2[i], ]
  assign(paste0("p", i), ggplot(dfsub[[i]], aes(x = beta1, y = power, group=group)) +
           geom_point(aes(color=group,shape=group)) +
           geom_line(aes(color=group,linetype=group)) +
           scale_linetype_manual(values = line_types) +  # Set linetype
           labs(title = bquote(beta[50] == .(beta2[i])),
                x = bquote("\u03B2"[49]),
                y = "Power",
                color = "",
                shape = "",
                linetype = "") +
           mytheme)
  mytheme
}

label <- textGrob("(c)", rot = 0, gp = gpar(fontsize = 12))
fig3 <- ggarrange(label,p1, p2, p3, p4, ncol = 5, nrow = 1,widths=c(0.1,0.2,0.2,0.2,0.2))
fig3

##---------------------------------Binary outcome, 50 correlated SNPs,MAF=0.05, rho=0.15, 0.1<=beta<=0.9

powerDerived<-powerD_deriveCor(kk=50,n=2000,alpha=0.05,p=0.05, rho=0.15,list1=beta1list1,list2=beta2list1)
result<-get(load("/path_to_the_output_file/snpDcor.Rdata"))
powerE2 <- sapply(result, function(x) x$powerE2)
powerCombined <- sapply(result, function(x) x$powerCombined)
powerSeparated <- sapply(result, function(x) x$powerSeparated)


df<-as.data.frame(cbind(as.vector(beta1list1),as.vector(beta2list1),powerSeparated,powerCombined,powerDerived,powerE2))

colnames(df)<-c("beta1","beta2","power1","power2","power3","power4")

power1<-cbind(df[,1:3],rep("Sim,True",25))
power2<-cbind(df[,c(1,2,4)],rep("Sim,Mis",25))
power3<-cbind(df[,c(1,2,5)],rep("Analytic,True",25))
power4<-cbind(df[,c(1,2,6)],rep("Analytic,Mis",25))


colnames(power1)[3:4]<-c("power","group")
colnames(power2)[3:4]<-c("power","group")
colnames(power3)[3:4]<-c("power","group")
colnames(power4)[3:4]<-c("power","group")


dfnew<-rbind(power1,power2,power3,power4)

beta2<-as.vector(c(0.1,0.25,0.5,0.75))
dfsub<-list()
line_types <- c("solid", "solid", "dashed", "dashed")

for (i in 1:4){
  dfsub[[i]] <- dfnew[dfnew$beta2==beta2[i], ]
  assign(paste0("p", i), ggplot(dfsub[[i]], aes(x = beta1, y = power, group=group)) +
           geom_point(aes(color=group,shape=group)) +
           geom_line(aes(color=group,linetype=group)) +
           scale_linetype_manual(values = line_types) +  # Set linetype
           labs(title = bquote(beta[50] == .(beta2[i])),
                x = bquote("\u03B2"[49]),
                y = "Power",
                color = "",
                shape = "",
                linetype = "") +
           mytheme)
  mytheme


}

label <- textGrob("(d)", rot = 0, gp = gpar(fontsize = 12))
fig4 <- ggarrange(label,p1, p2, p3, p4, ncol = 5, nrow = 1,widths=c(0.1,0.2,0.2,0.2,0.2))
fig4


#-------------------------------------------------combine figure 1-4

dummy_plot <- ggplot(dfnew, aes(x = beta1, y = power,
                                color = group,
                                linetype = group,
                                shape = group)) +
  geom_point() +
  geom_line() +
  scale_linetype_manual(values = c("solid", "solid", "dashed", "dashed")) +
  theme(legend.position = "bottom") +
  labs(color = NULL, linetype = NULL, shape = NULL)  # Remove legend titles

# Extract the legend (now title-less)
legend <- cowplot::get_legend(dummy_plot)
combined_fig <- ggarrange(
  fig1, fig2, fig3, fig4,
  nrow = 4, ncol = 1,
  heights = c(1, 1, 1, 1)  # Adjust heights if needed
)

# Add the legend to the combined figure
final_plot <- ggarrange(
  combined_fig, legend,
  nrow = 2, ncol = 1,
  heights = c(10, 1)  # Allocate 10:1 ratio for plot vs legend
)

print(final_plot)


