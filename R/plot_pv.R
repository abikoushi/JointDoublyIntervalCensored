library(reshape2)
library(ggplot2)
library(dplyr)
library(readr)
path = "./jl/simdata/outsim_trunc_pv_b.csv"
out_pv_B <- read_csv(path)
dfB1 = dplyr::filter(out_pv_B, dist=="LogNormal")
dfB2 = dplyr::filter(out_pv_B, dist=="MixtureModel")
WSs=c(2,5,10)
tdists = c("LogNormal", "Mixture")

#pdf("pvalue_b.pdf", width=10, height=10)
p = ggplot(dfB1, aes(x=pv))+
  stat_ecdf()+
  geom_abline(slope = 1, intercept = 0, colour="grey", linetype=2)+
  facet_grid(tau+sigma ~ WS, labeller = label_both)+
  scale_x_continuous(n.breaks = 11)+
  theme_classic(16)+ labs(x="nominal", y="actual", title = "log-normal")+
  theme(strip.text.y = element_text(angle=0), axis.text.x = element_text(angle=90))
#print(p) 
  
p = ggplot(dfB2, aes(x=pv))+
  stat_ecdf()+
  geom_abline(slope = 1, intercept = 0, colour="grey", linetype=2)+
  facet_grid(tau+sigma ~ WS, labeller = label_both)+
  scale_x_continuous(n.breaks = 11)+
  theme_classic(16)+ labs(x="nominal", y="actual", title = "Mixture")+
  theme(strip.text.y = element_text(angle=0), axis.text.x = element_text(angle=90))
#print(p)
#dev.off()
####
#intensity
path = "./jl/simdata/outsim_trunc_pv_int.csv"
out_pv_int <- read_csv(path)
out_pv_int = rename(out_pv_int, quantile=pos)
#head(out_pv_int)
#3*3*3*4*100
df_int1 = dplyr::filter(out_pv_int, dist=="LogNormal")
df_int2 = dplyr::filter(out_pv_int, dist=="MixtureModel")

df_int_1_s = split(df_int1, df_int1$WS)
df_int_2_s = split(df_int2, df_int2$WS)
i=1
#pdf("pvalue_intensity.pdf", width=10, height=10)
for(i in 1:3){
p_i <- ggplot(df_int_1_s[[i]], aes(x=pv))+
  stat_ecdf()+
  geom_abline(slope = 1, intercept = 0, colour="grey", linetype=2)+
  facet_grid(tau+sigma ~ quantile, labeller = label_both)+
  scale_x_continuous(n.breaks = 11)+
  theme_classic(16)+ labs(x="nominal", y="actual", title = paste("log-normal dist. WS:",WSs[i]))+
  theme(strip.text.y = element_text(angle=0), axis.text.x = element_text(angle=90))  
print(p_i)
ggsave(plot = p_i, filename = paste0("pv_infection_lnorm_",i,".pdf"), width=10, height = 10)
}

for(i in 1:3){
  p_i <- ggplot(df_int_2_s[[i]], aes(x=pv))+
    stat_ecdf()+
    geom_abline(slope = 1, intercept = 0, colour="grey", linetype=2)+
    facet_grid(tau+sigma ~ quantile, labeller = label_both)+
    scale_x_continuous(n.breaks = 11)+
    theme_classic(16)+ labs(x="nominal", y="actual", title = paste("mixture dist. WS:",WSs[i]))+
    theme(strip.text.y = element_text(angle=0), axis.text.x = element_text(angle=90))  
  print(p_i)
  ggsave(plot = p_i, filename = paste0("pv_infection_mixture_",i,".pdf"), width=10, height = 10)
}
#dev.off()
###
#incubation
path = "./jl/simdata/outsim_trunc_pv_ccdf.csv"
out_pv_incubation <- read_csv(path)
out_pv_incubation = rename(out_pv_incubation, quantile=pos)
#dim(out_pv_int)
#3*3*3*4*100
df_inc1 = dplyr::filter(out_pv_incubation, dist=="LogNormal")
df_inc2 = dplyr::filter(out_pv_incubation, dist=="MixtureModel")
# df_inc1 = out_pv_incubation[1:10800,]
# df_inc2 = out_pv_incubation[10801:21600,]

df_inc_1_s = split(df_int1, df_int1$WS)
df_inc_2_s = split(df_int2, df_int2$WS)
#pdf("pvalue_ccdf.pdf", width=10, height=10)
for(i in 1:3){
  p_i <- ggplot(df_inc_1_s[[i]], aes(x=pv))+
    stat_ecdf()+
    geom_abline(slope = 1, intercept = 0, colour="grey", linetype=2)+
    facet_grid(tau+sigma ~ quantile, labeller = label_both)+
    scale_x_continuous(n.breaks = 11)+
    theme_classic(16)+ labs(x="nominal", y="actual", title = paste("log-normal dist. WS:",WSs[i]))+
    theme(strip.text.y = element_text(angle=0), axis.text.x = element_text(angle=90))  
  print(p_i)  
  ggsave(plot = p_i, filename = paste0("pv_incubation_lnorm_",i,".pdf"), width=10, height = 10)
}

for(i in 1:3){
  p_i <- ggplot(df_inc_2_s[[i]], aes(x=pv))+
    stat_ecdf()+
    geom_abline(slope = 1, intercept = 0, colour="grey", linetype=2)+
    facet_grid(tau+sigma ~ quantile, labeller = label_both)+
    scale_x_continuous(n.breaks = 11)+
    theme_classic(16)+ labs(x="nominal", y="actual", title = paste("mixture dist. WS:",WSs[i]))+
    theme(strip.text.y = element_text(angle=0), axis.text.x = element_text(angle=90))  
  print(p_i)  
  ggsave(plot = p_i, filename = paste0("pv_incubation_mixture_",i,".pdf"), width=10, height = 10)
}
#dev.off()
#ggsave(plot = p_i, filename = paste0("pv_incubation_mixture_",i,".pdf"))
