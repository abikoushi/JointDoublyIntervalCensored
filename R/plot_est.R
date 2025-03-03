library(readr)
library(dplyr)
library(ggplot2)
library(ggh4x)
#LogNormal <- list(meanlog=log(10)-0.5, sdlog=1, lower.tail=FALSE)
tccdf1 = function(x){
  plnorm(x, log(10)-0.5, 1, lower.tail = FALSE)
}

tccdf2 = function(x){
  r = c(2/3,1/3)
  r[1]*pweibull(x, 2, 5/gamma(1+1/2), lower.tail = FALSE) + 
    r[2]*pweibull(x,0.5, 10, lower.tail = FALSE)
}

Lambda <- function(x,sigma){
  x^sigma
}

path = "./jl/outsim_trunc_gibbs.csv"
res_gibbs = read_csv(path)
res_gibbs_s <- split(res_gibbs, res_gibbs$dist)


pdf("ccdf.pdf", width = 10, height=10)
p = ggplot(res_gibbs_s[[1]],aes(x=value, y=ccdf))+
  geom_step(aes(group=id), linewidth=0.1, alpha=0.1)+
  stat_function(fun=tccdf1, linetype=2, colour="royalblue", linewidth=1)+
  facet_grid(sigma+tau ~ WS, scales="free_x", labeller = label_both)+
  theme_classic(16) + labs(title = "log-normal")+
  theme(strip.text.y = element_text(angle = 0))
print(p)
#ggsave("ccdf_ln.pdf")

p = ggplot(res_gibbs_s[[2]],aes(x=value,y=ccdf))+
  geom_step(aes(group=id), linewidth=0.1, alpha=0.1)+
  stat_function(fun=tccdf2, linetype=2, colour="royalblue", linewidth=1)+
  facet_grid(sigma+tau ~ WS, scales="free_x", labeller = label_both)+
  theme_classic(16) + labs(title = "mixture")+
  theme(strip.text.y = element_text(angle = 0))
print(p)
dev.off()
#ggsave("ccdf_mixt.pdf")

sigmas = c(0.8, 1, 1.2)
res_gibbs_s_1 <- split(res_gibbs_s[[1]], res_gibbs_s[[1]]$sigma)
res_gibbs_s_2 <- split(res_gibbs_s[[2]], res_gibbs_s[[2]]$sigma)

####
#intensity
pdf("intensity.pdf")
for(i in 1:3){
  p = ggplot(res_gibbs_s_1[[i]], aes(x=value, y=intensity))+
  geom_step(aes(group=id), linewidth=0.1, alpha=0.1)+
  stat_function(fun=Lambda, args = list(sigma=sigmas[i]), linetype=2, colour="royalblue", linewidth=1)+
  facet_grid2(tau~WS, scales="free", labeller = label_both, independent = "x")+
  theme_classic(16) + labs(title = paste0("sigma = ",sigmas[i],"; log-normal"),x="time", y="cumulative intensity")+
  theme(strip.text.y = element_text(angle = 0))
  print(p)
}
#ggsave(paste0("cumint_ln_sigma",sigmas[i],".pdf"))

for(i in 1:3){
  p <- ggplot(res_gibbs_s_2[[i]], aes(x=value, y=intensity))+
    geom_step(aes(group=id), linewidth=0.1, alpha=0.1)+
    stat_function(fun=Lambda, args = list(sigma=sigmas[i]), linetype=2, colour="royalblue", linewidth=1)+
    facet_grid2(tau~WS, scales="free", labeller = label_both, independent = "x")+
    theme_classic(16) + labs(title = paste0("sigma = ",sigmas[i],"; mixture"),x="time", y="cumulative intensity")+
    theme(strip.text.y = element_text(angle = 0))  
  print(p)
}
dev.off()
#ggsave("cumint_mixt.pdf")

####
#bias & se
resEM = read_csv("./jl/outsim_trunc_em.csv")
resEM_s <- split(resEM, resEM$dist)

resVB = read_csv("./jl/outsim_trunc_vb.csv")
resVB_s <- split(resEM, resEM$dist)

df1_gibbs <- group_by(res_gibbs_s[[1]], dist, sigma, tau, WS) %>% 
  summarise(RMSE = sqrt(mean((ccdf-tccdf1(value))^2)),
            bias = mean(ccdf-tccdf1(value)),
            SE=sd(ccdf)) %>% 
  ungroup() %>% 
  mutate(method="EAP")

df1_em <- group_by(resEM_s[[1]], dist, sigma, tau, WS) %>% 
  summarise(RMSE = sqrt(mean((ccdf-tccdf1(value))^2)),
            bias = mean(ccdf-tccdf1(value)),
            SE=sd(ccdf)) %>% 
  ungroup() %>% 
  mutate(method="MLE")

df1_vb <- group_by(resVB_s[[1]], dist, sigma, tau, WS) %>% 
  summarise(RMSE = sqrt(mean((ccdf-tccdf1(value))^2)),
            bias = mean(ccdf-tccdf1(value)),
            SE=sd(ccdf)) %>% 
  ungroup() %>% 
  mutate(method="VEAP")

df1 = bind_rows(df1_gibbs, df1_em, df1_vb)

df2_gibbs <- group_by(res_gibbs_s[[2]], dist, sigma, tau, WS) %>% 
  summarise(RMSE = sqrt(mean((ccdf-tccdf2(value))^2)), 
            bias = mean(ccdf-tccdf2(value)),
            SE=sd(ccdf)) %>% 
  ungroup() %>% 
  mutate(method="EAP")

df2_em <- group_by(resEM_s[[2]], dist, sigma, tau, WS) %>% 
  summarise(RMSE = sqrt(mean((ccdf-tccdf2(value))^2)),
            bias = mean(ccdf-tccdf2(value)),
            SE=sd(ccdf)) %>% 
  ungroup() %>% 
  mutate(method="MLE")

df2_vb <- group_by(resVB_s[[2]], dist, sigma, tau, WS) %>% 
  summarise(RMSE = sqrt(mean((ccdf-tccdf1(value))^2)),
            bias = mean(ccdf-tccdf1(value)),
            SE=sd(ccdf)) %>% 
  ungroup() %>% 
  mutate(method="VEAP")

df2 = bind_rows(df2_gibbs, df2_em, df2_vb)

df_int_gibbs <- group_by(res_gibbs, dist, sigma, tau, WS) %>% 
  summarise(RMSE = sqrt(mean((intensity-Lambda(value,sigma))^2)),
            bias = mean(intensity-Lambda(value,sigma)),
            SE = sd(intensity)) %>% 
  ungroup() %>% 
  mutate(method="EAP", dist=gsub("Model","",dist))

df_int_em <- group_by(resEM, dist, sigma, tau, WS) %>% 
  summarise(RMSE = sqrt(mean((intensity-Lambda(value,sigma))^2)),
            bias = mean(intensity-Lambda(value,sigma)),
            SE=sd(intensity)) %>% 
  ungroup() %>% 
  mutate(method="MLE", dist=gsub("Model","",dist))

df_int_vb <- group_by(resVB, dist, sigma, tau, WS) %>% 
  summarise(RMSE = sqrt(mean((intensity-Lambda(value,sigma))^2)),
            bias = mean(intensity-Lambda(value,sigma)),
            SE=sd(intensity)) %>% 
  ungroup() %>% 
  mutate(method="VEAP", dist=gsub("Model","",dist))

df_int = bind_rows(df_int_gibbs, df_int_em, df_int_vb)

pdf("bias_and_se.pdf", width=10, height=12)
p = ggplot(df2, aes(x=factor(WS), y=bias, colour=method, shape = method))+
  geom_pointrange(aes(ymin=bias-SE, ymax=bias+SE), position = position_dodge(width = 0.5))+
  geom_hline(yintercept = 0, linetype=3)+
  facet_grid(tau~sigma, labeller=label_both) +
  scale_colour_grey()+scale_shape_manual(values = c(15:17))+
  theme_classic(16) + 
  labs(title = "CCDF: mixture dist.", x="WS", y="bias and se")
#ggsave("ccdf_bias_se_mixt.pdf")
print(p)

p = ggplot(df1, aes(x=factor(WS), y=bias, colour=method, shape = method))+
  geom_pointrange(aes(ymin=bias-SE, ymax=bias+SE), position = position_dodge(width = 0.5))+
  geom_hline(yintercept = 0, linetype=3)+
  facet_grid(tau~sigma, labeller=label_both) +
  scale_colour_grey()+scale_shape_manual(values = c(15:17))+
  theme_classic(16) + 
  labs(title = "CCDF: log-normal dist.", x="WS", y="bias and se")
#ggsave("ccdf_bias_se_lnorm.pdf")
print(p)

p = ggplot(df_int, aes(x=factor(WS), y=bias, colour=method, shape = method))+
  geom_pointrange(aes(ymin=bias-SE, ymax=bias+SE), position = position_dodge(width = 0.5))+
  geom_hline(yintercept = 0, linetype=3)+
  facet_grid2(tau+sigma~dist, labeller=label_both, scales="free_y") +
  scale_colour_grey()+scale_shape_manual(values = c(15:17))+
  theme_classic(16) + theme(strip.text.y = element_text(angle=0))+
  labs(title = "Cumulative intensity", x="WS", y="bias and se")
print(p)
#ggsave("intensity_bias_se.pdf")
dev.off()
