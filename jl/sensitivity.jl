#=
Section 4 "Case Study"
Fig 
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
@time res_gibbs0 = jointest.gibbssampler(dat.EL, dat.ER, dat.SL, dat.SR, maximum(dat.SR)+0, 4000)
@time res_gibbs1 = jointest.gibbssampler(dat.EL, dat.ER, dat.SL, dat.SR, maximum(dat.SR)+1, 4000)
@time res_gibbs2 = jointest.gibbssampler(dat.EL, dat.ER, dat.SL, dat.SR, maximum(dat.SR)+2, 4000)
@time res_gibbs3 = jointest.gibbssampler(dat.EL, dat.ER, dat.SL, dat.SR, maximum(dat.SR)+3, 4000)
# p = plot(2001:4000, res_gibbs1.logprob[2001:4000], tickdirection=:out, legend=false, color="grey")
# xlabel!(p, "step")
# ylabel!(p, "log-probability")

Shat0 = vec(mean(res_gibbs0.ccdf[2001:end,:], dims=1))
Shat1 = vec(mean(res_gibbs1.ccdf[2001:end,:], dims=1))
Shat2 = vec(mean(res_gibbs2.ccdf[2001:end,:], dims=1))
Shat3 = vec(mean(res_gibbs3.ccdf[2001:end,:], dims=1))
p = plot(res_gibbs0.value, Shat0, label="0")
plot!(p, res_gibbs1.value, Shat1, label="1")
plot!(p, res_gibbs2.value, Shat1, label="2")
plot!(p, res_gibbs3.value, Shat1, label="3")
xlabel!(p, "day")
ylabel!(p, "CCDF")
savefig(p, "CCDF.pdf")

Lhat0= vec(mean(res_gibbs0.intensity[2001:end,:], dims=1))
Lhat1= vec(mean(res_gibbs1.intensity[2001:end,:], dims=1))
Lhat2= vec(mean(res_gibbs2.intensity[2001:end,:], dims=1))
Lhat3= vec(mean(res_gibbs3.intensity[2001:end,:], dims=1))
p = plot(res_gibbs0.value, Lhat0, label="0")
plot!(p, res_gibbs1.value, Lhat1, label="1")
plot!(p, res_gibbs2.value, Lhat2, label="2")
plot!(p, res_gibbs3.value, Lhat3, label="3")
xlabel!(p, "day")
ylabel!(p, "intensity")
savefig(p, "intensity.pdf")
