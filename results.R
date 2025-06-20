# generate graphs for simulated power and theoretical power under four scenarios

library(SKAT)
library(bindata)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(KATSP)
library(ggpubr)

beta1list <- list(
  c(-0.9, -0.75, -0.5, -0.25, -0.1, -0.9, -0.75, -0.5, -0.25, -0.1,
    -0.9, -0.75, -0.5, -0.25, -0.1, -0.9, -0.75, -0.5, -0.25, -0.1,
    -0.9, -0.75, -0.5, -0.25, -0.1),
  c(-1.5, -1.25, -1, -0.75, -0.5, -1.5, -1.25, -1, -0.75, -0.5,
    -1.5, -1.25, -1, -0.75, -0.5, -1.5, -1.25, -1, -0.75, -0.5,
    -1.5, -1.25, -1, -0.75, -0.5),
  c(-3, -2.5, -2, -1.5, -1, -3, -2.5, -2, -1.5, -1,
    -3, -2.5, -2, -1.5, -1, -3, -2.5, -2, -1.5, -1,
    -3, -2.5, -2, -1.5, -1)
)

beta2list <- list(
  c(rep(0.9, 5), rep(0.75, 5), rep(0.5, 5), rep(0.25, 5), rep(0.1, 5)),
  c(rep(1.5, 5), rep(1.25, 5), rep(1, 5), rep(0.75, 5), rep(0.5, 5)),
  c(rep(3, 5), rep(2.5, 5), rep(2, 5), rep(1.5, 5), rep(1, 5))
)

beta1list1<-c(-0.9,-0.75,-0.5,-0.25,-0.1,
             -0.9,-0.75,-0.5,-0.25,-0.1,
             -0.9,-0.75,-0.5,-0.25,-0.1,
             -0.9,-0.75,-0.5,-0.25,-0.1,
             -0.9,-0.75,-0.5,-0.25,-0.1)
beta2list1<-c( rep(0.9,5),rep(0.75,5),rep(0.5,5),rep(0.25,5),rep(0.1,5))

beta1list2<-c(-1.5,-1.25,-1,-0.75,-0.5,
             -1.5,-1.25,-1,-0.75,-0.5,
             -1.5,-1.25,-1,-0.75,-0.5,
             -1.5,-1.25,-1,-0.75,-0.5,
             -1.5,-1.25,-1,-0.75,-0.5)
beta2list2<-c( rep(1.5,5),rep(1.25,5),rep(1,5),rep(0.75,5),rep(0.5,5))

beta1list3<-c(-3,-2.5,-2,-1.5,-1,
             -3,-2.5,-2,-1.5,-1,
             -3,-2.5,-2,-1.5,-1,
             -3,-2.5,-2,-1.5,-1,
             -3,-2.5,-2,-1.5,-1)
beta2list3<-c( rep(3,5),rep(2.5,5),rep(2,5),rep(1.5,5),rep(1,5))

theme_set(theme_cowplot())
mytheme = list(
  theme_classic()+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          legend.position="none", )
)

## Continuous outcome, 50 uncorrelated SNPs, MAF=0.01, 0.1<=beta<=0.9

powerDerived<-powerC_derive(k=50,n=2000,alpha=0.05,p=0.01,list1=beta1list1,list2=beta2list1)
result<-get(load("/path_to_the_file/snpC.Rdata"))
powerE2 <- sapply(result, function(x) x$powerE2)
powerCombined <- sapply(result, function(x) x$powerCombined)
powerSeparated <- sapply(result, function(x) x$powerSeparated)

df<-as.data.frame(cbind(as.vector(beta1list1),as.vector(beta2list1),powerSeparated,powerCombined,powerDerived,powerE2))
colnames(df)<-c("beta1","beta2","power1","power2","power3","power4")

power1<-cbind(df[,1:3],rep("S_true",25))
power2<-cbind(df[,c(1,2,4)],rep("S_misspecified",25))
power3<-cbind(df[,c(1,2,5)],rep("C_true",25))
power4<-cbind(df[,c(1,2,6)],rep("C_misspecified",25))

colnames(power1)[3:4]<-c("power","group")
colnames(power2)[3:4]<-c("power","group")
colnames(power3)[3:4]<-c("power","group")
colnames(power4)[3:4]<-c("power","group")

dfnew<-rbind(power1,power2,power3,power4)
beta2<-as.vector(c(0.1,0.25,0.5,0.75,0.9))
dfsub<-list()
for (i in 1:5){
  dfsub[[i]] <- dfnew[dfnew$beta2==beta2[i], ]
  assign(paste0("p", i), ggplot(dfsub[[i]], aes(x = beta1, y = power, group=group)) +
          geom_point(aes(color=group,shape=group)) +
          geom_line(aes(color=group,linetype=group)) +
            labs(title =bquote(beta[2] == .(beta2[i])),
                x=bquote("\u03B2"[1]),
                y = "power")) +
    mytheme
}

ggarrange(p1, p2, p3, p4,p5, ncol=3, nrow=2, common.legend = TRUE, legend="bottom")

## Continuous outcome, 50 correlated SNPs,MAF=0.05, rho=0.15, 0.1<=beta<=0.9

result<-get(load("/path_to_the_file/snpCcor.Rdata"))
powerE2 <- sapply(result, function(x) x$powerE2)
powerCombined <- sapply(result, function(x) x$powerCombined)
powerSeparated <- sapply(result, function(x) x$powerSeparated)
powerDerived<- powerC_deriveCor(k=50,n=2000,alpha=0.05,p=0.05,rho=0.15,list1=beta1list1,list2=beta2list1)

df<-as.data.frame(cbind(as.vector(beta1list1),as.vector(beta2list1),powerSeparated,powerCombined,powerDerived,powerE2))
colnames(df)<-c("beta1","beta2","power1","power2","power3","power4")

power1<-cbind(df[,1:3],rep("S_true",25))
power2<-cbind(df[,c(1,2,4)],rep("S_misspecified",25))
power3<-cbind(df[,c(1,2,5)],rep("C_true",25))
power4<-cbind(df[,c(1,2,6)],rep("C_misspecified",25))

colnames(power1)[3:4]<-c("power","group")
colnames(power2)[3:4]<-c("power","group")
colnames(power3)[3:4]<-c("power","group")
colnames(power4)[3:4]<-c("power","group")

dfnew<-rbind(power1,power2,power3,power4)
beta2<-as.vector(c(0.1,0.25,0.5,0.75,0.9))


dfsub<-list()


for (i in 1:5){
  dfsub[[i]] <- dfnew[dfnew$beta2==beta2[i], ]
  assign(paste0("p", i), ggplot(dfsub[[i]], aes(x = beta1, y = power, group=group)) +
          geom_point(aes(color=group,shape=group)) +
          geom_line(aes(color=group,linetype=group)) +
            labs(title =bquote(beta[2] == .(beta2[i])),
                x=bquote("\u03B2"[1]),
                y = "power")) +
    mytheme


}
ggarrange(p1, p2, p3, p4,p5, ncol=3, nrow=2, common.legend = TRUE, legend="bottom")

## Binary outcome, 50 uncorrelated SNPs, MAF=0.01, 0.1<=beta<=0.9

powerDerived<-powerD_derive(kk=50,n=2000,alpha=0.05,p=0.01,list1=beta1list1,list2=beta2list1)
result<-get(load("/path_to_the_output_file/snpD.Rdata"))
powerE2 <- sapply(result, function(x) x$powerE2)
powerCombined <- sapply(result, function(x) x$powerCombined)
powerSeparated <- sapply(result, function(x) x$powerSeparated)
df<-as.data.frame(cbind(as.vector(beta1list1),as.vector(beta2list1),powerSeparated,powerCombined,powerDerived,powerE2))
colnames(df)<-c("beta1","beta2","power1","power2","power3","power4")

power1<-cbind(df[,1:3],rep("S_true",25))
power2<-cbind(df[,c(1,2,4)],rep("S_misspecified",25))
power3<-cbind(df[,c(1,2,5)],rep("C_true",25))
power4<-cbind(df[,c(1,2,6)],rep("C_misspecified",25))

colnames(power1)[3:4]<-c("power","group")
colnames(power2)[3:4]<-c("power","group")
colnames(power3)[3:4]<-c("power","group")
colnames(power4)[3:4]<-c("power","group")

dfnew<-rbind(power1,power2,power3,power4)
beta2<-as.vector(c(0.1,0.25,0.5,0.75,0.9))
dfsub<-list()
for (i in 1:5){
  dfsub[[i]] <- dfnew[dfnew$beta2==beta2[i], ]
  assign(paste0("p", i), ggplot(dfsub[[i]], aes(x = beta1, y = power, group=group)) +
          geom_point(aes(color=group,shape=group)) +
          geom_line(aes(color=group,linetype=group)) +
            labs(title =bquote(beta[2] == .(beta2[i])),
                x=bquote("\u03B2"[1]),
                y = "power")) +
    mytheme
}
ggarrange(p1, p2, p3, p4,p5, ncol=3, nrow=2, common.legend = TRUE, legend="bottom")

## Binary outcome, 50 correlated SNPs,MAF=0.05, rho=0.15, 0.1<=beta<=0.9

powerDerived<-powerD_deriveCor(kk=50,n=2000,alpha=0.05,p=0.05, rho=0.15,list1=beta1list1,list2=beta2list1)
result<-get(load("/path_to_the_output_file/snpDcor.Rdata"))
powerE <- sapply(result, function(x) x$powerE)
powerCombined <- sapply(result, function(x) x$powerCombined)
powerSeparated <- sapply(result, function(x) x$powerSeparated)

df<-as.data.frame(cbind(as.vector(beta1list1),as.vector(beta2list1),powerSeparated,powerCombined,powerDerived,powerE2))
colnames(df)<-c("beta1","beta2","power1","power2","power3","power4")

power1<-cbind(df[,1:3],rep("S_true",25))
power2<-cbind(df[,c(1,2,4)],rep("S_misspecified",25))
power3<-cbind(df[,c(1,2,5)],rep("C_true",25))
power4<-cbind(df[,c(1,2,6)],rep("C_misspecified",25))

colnames(power1)[3:4]<-c("power","group")
colnames(power2)[3:4]<-c("power","group")
colnames(power3)[3:4]<-c("power","group")
colnames(power4)[3:4]<-c("power","group")

dfnew<-rbind(power1,power2,power3,power4)
beta2<-as.vector(c(0.1,0.25,0.5,0.75,0.9))
dfsub<-list()

#p1-p5
for (i in 1:5){
  dfsub[[i]] <- dfnew[dfnew$beta2==beta2[i], ]
  assign(paste0("p", i), ggplot(dfsub[[i]], aes(x = beta1, y = power, group=group)) +
          geom_point(aes(color=group,shape=group)) +
          geom_line(aes(color=group,linetype=group)) +
            labs(title =bquote(beta[2] == .(beta2[i])),
                x=bquote("\u03B2"[1]),
                y = "power")) +
    mytheme
}
ggarrange(p1, p2, p3, p4,p5, ncol=3, nrow=2, common.legend = TRUE, legend="bottom")




