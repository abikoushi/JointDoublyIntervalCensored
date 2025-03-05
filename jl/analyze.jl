#=
Section 3 "Simulation Study"
Fig 4 - 5 
=#

using LinearAlgebra
using DataFrames
using CSV
using Plots
using Dates
using Distributions
using LogExpFunctions
using RCall
using StatsBase
using Random
using TexTables
include("./JointTools.jl")
include("./JointTools_Gibbs.jl")
dat = CSV.read("./data/covid19/dat_t.csv", DataFrame)



Random.seed!(1234)
res_gibbs = jointest.gibbssampler(dat.EL, dat.ER, dat.SL, dat.SR, maximum(dat.SR), 4000)

p = plot(2001:4000, res_gibbs.logprob[2001:4000], tickdirection=:out, legend=false, color="grey")
xlabel!(p, "step")
ylabel!(p, "log-probability")
savefig(p, "traceplot.pdf")


Shat = vec(mean(res_gibbs.ccdf[2001:end,:], dims=1))
CI_s = jointtool.takeCI(res_gibbs.ccdf[2001:end,:], [0.025, 0.975])
plot(res_gibbs.value, Shat)


df = DataFrame(value=res_gibbs.value, ccdf=Shat, lower = CI_s[:,1], upper = CI_s[:,2])

R"
library(coarseDataTools)
library(dplyr)

dic_type <- c(2,1,1,0)
dat_t2 <- mutate($dat_t, type=dic_type[ctype+1L]) %>% 
  mutate(EL = ifelse(type==0,0,EL)) %>% 
  as.data.frame()

covid_opt_w <- optim(c(0.1,2),loglikhd,dat=dat_t2,dist = 'W')
covid_opt_w <- optim(covid_opt_w$par,loglikhd,dat=dat_t2,dist = 'W', method = 'BFGS')
covid_opt_g <- optim(c(1.5,0.5),loglikhd,dat=dat_t2,dist = 'G')
covid_opt_g <- optim(covid_opt_g$par,loglikhd,dat=dat_t2,dist = 'G', method = 'BFGS')
covid_opt_l <- optim(c(1,0),loglikhd,dat=dat_t2,dist = 'L')
covid_opt_l <- optim(covid_opt_l$par,loglikhd,dat=dat_t2,dist = 'L', method = 'BFGS')
"

R"
library(ggplot2)
library(pammtools)
ggplot($df, aes(x=value))+
  geom_stepribbon(aes(ymin=lower, ymax=upper), alpha=0.1)+
  geom_step(aes(y=ccdf))+
  stat_function(fun=pgamma, args = list(shape=exp(covid_opt_g$par[1]), scale=exp(covid_opt_g$par[2]), lower.tail=FALSE),
                aes(linetype='gamma', colour='gamma'))+
  stat_function(fun=pgamma, args = list(shape=exp(covid_opt_g$par[1]), scale=exp(covid_opt_g$par[2]), lower.tail=FALSE),
                geom = 'point', n=20, aes(shape='gamma', colour='gamma'), size=2.5)+
  stat_function(fun=plnorm, args = list(meanlog=covid_opt_l$par[1], sdlog=exp(covid_opt_l$par[2]), lower.tail=FALSE),
                aes(linetype='log normal', colour='log normal'))+
  stat_function(fun=plnorm, args = list(meanlog=covid_opt_l$par[1], sdlog=exp(covid_opt_l$par[2]), lower.tail=FALSE),
                geom = 'point', n=20, aes(shape='log normal', colour='log normal'), size=2.5)+
  stat_function(fun=pweibull, args = list(shape=exp(covid_opt_w$par[1]), scale=exp(covid_opt_w$par[2]), lower.tail=FALSE),
                aes(linetype='weibull', colour='weibull'))+
  stat_function(fun=pweibull, args = list(shape=exp(covid_opt_w$par[1]), scale=exp(covid_opt_w$par[2]), lower.tail=FALSE),
                geom = 'point', n=20, aes(shape='weibull', colour='weibull'), size=2.5)+
  labs(x='incubation period', y='ccdf', shape='model', linetype='model', colour='model')+
  scale_colour_manual(values = c('orange2', 'royalblue', 'forestgreen'))+
  theme_classic(16)
  ggsave('fit_covid.pdf')
"


Lhat = vec(mean(res_gibbs.intensity[2001:end,:], dims=1))
CI_lambda = jointtool.takeCI(res_gibbs.intensity[2001:end,:], [0.025, 0.975])

bhat = vec(mean(res_gibbs.bt[2001:end,:], dims=1))
CI_b = jointtool.takeCI(res_gibbs.bt[2001:end,:], [0.025, 0.975])

df = DataFrame(date = Dates.Day.(res_gibbs.value) .+ dat_t.recmin[1], EAP = Lhat,lower = CI_lambda[:,1], upper = CI_lambda[:,2], bt=bhat, lower_b=CI_b[:,1], upper_b = CI_b[:,2])

R"
library(ggplot2)
library(pammtools)
ggplot($df, aes(x=date))+
  geom_stepribbon(aes(ymin=lower, ymax=upper), alpha=0.1)+
  geom_step(aes(y=EAP))+
  geom_stepribbon(aes(ymin=lower_b, ymax=upper_b, fill='truncated'), alpha=0.1)+
  geom_step(aes(y=bt, colour='truncated'), linetype=2)+
  labs(y='cumulative intensity', colour='', fill='')+
  theme_classic(16)+scale_x_date(date_labels = '%b/%d')
  ggsave('fit_covid_intensity.pdf')
"



####
#shifted estimator
EL = dat_t.EL - dat_t.EL
ER = dat_t.ER - dat_t.EL
SL = dat_t.SL - dat_t.EL
SR = dat_t.SR - dat_t.EL

res0 = jointest.gibbssampler(EL, ER, SL, SR, Inf, 4000)

Shat = vec(mean(res0.ccdf[2001:end,:], dims=1))
CI_s = jointtool.takeCI(res0.ccdf[2001:end,:], [0.025, 0.975])
plot(res0.value, Shat)

df0 = DataFrame(value=res0.value, ccdf=Shat, lower = CI_s[:,1], upper = CI_s[:,2])

R"
library(ggplot2)
library(pammtools)
ggplot($df0, aes(x=value))+
  geom_stepribbon(aes(ymin=lower, ymax=upper), alpha=0.1)+
  geom_step(data=$df0, aes(y=ccdf))+
  stat_function(fun=pgamma, args = list(shape=exp(covid_opt_g$par[1]), scale=exp(covid_opt_g$par[2]), lower.tail=FALSE),
  aes(linetype='gamma', colour='gamma'))+
stat_function(fun=pgamma, args = list(shape=exp(covid_opt_g$par[1]), scale=exp(covid_opt_g$par[2]), lower.tail=FALSE),
  geom = 'point', n=20, aes(shape='gamma', colour='gamma'), size=2.5)+
stat_function(fun=plnorm, args = list(meanlog=covid_opt_l$par[1], sdlog=exp(covid_opt_l$par[2]), lower.tail=FALSE),
  aes(linetype='log normal', colour='log normal'))+
stat_function(fun=plnorm, args = list(meanlog=covid_opt_l$par[1], sdlog=exp(covid_opt_l$par[2]), lower.tail=FALSE),
  geom = 'point', n=20, aes(shape='log normal', colour='log normal'), size=2.5)+
stat_function(fun=pweibull, args = list(shape=exp(covid_opt_w$par[1]), scale=exp(covid_opt_w$par[2]), lower.tail=FALSE),
  aes(linetype='weibull', colour='weibull'))+
stat_function(fun=pweibull, args = list(shape=exp(covid_opt_w$par[1]), scale=exp(covid_opt_w$par[2]), lower.tail=FALSE),
  geom = 'point', n=20, aes(shape='weibull', colour='weibull'), size=2.5)+
labs(x='incubation period', y='ccdf', shape='model', linetype='model', colour='model')+
scale_colour_manual(values = c('orange2', 'royalblue', 'forestgreen'))+
theme_bw(16)
ggsave('fit_covid0.pdf')
"
