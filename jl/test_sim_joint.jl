using Random
using LinearAlgebra
using Distributions
using Statistics
using Plots
using SpecialFunctions
using DataFrames
using RCall
using Profile
using XLSX
using JLD
using SparseArrays
using HCubature
using CSV
using StatsBase

include("./JointTools_Gibbs.jl")
include("./JointTools.jl")

Lambda(t, shape) = t^shape
lambda(t, shape) = shape * (t^(shape - 1.0))

function integrand(X, d, tau, WS, sigma)
    Distributions.ccdf(d, tau - X[1]*WS - X[2]) * lambda(X[2], sigma)
end

tdist1 = LogNormal(log(10) - 0.5, 1)
tdist2 = MixtureModel([Weibull(2, 5/gamma(1+1/2)), Weibull(0.5, 10)], [2/3,1/3])
tdists = [tdist1, tdist2]

tdist = LogNormal(log(10)-0.5, 1)
tau=100.0
sigma=1.0
WE=2.0
WS=2.0
rng =Random.default_rng()
d = jointtool.makeDIC(rng, tdist, tau, WE, WS, sigma)
f = d[4] .<= tau
LE = floor.(d[1][f])
RE = ceil.(d[2][f])
LS = floor.(d[3][f])
RS = ceil.(d[4][f])
res = jointest.gibbssampler(LE, RE, LS, RS, tau, 4000)

plot(res.logprob, legend=false)
Shat = vec(mean(res.ccdf[2001:end,:], dims=1))
Lhat = vec(mean(res.intensity[2001:end,:], dims=1))

prob = [0.5, 0.8, 0.9, 0.95]
qs = quantile(tdist, prob)
pv = jointtool.simulated_pvalue1(res.ccdf[2001:end,:], res.value, qs, 1 .- prob)

p = plot(res.value, Shat, legend=false)
plot!(p, x -> ccdf(tdist,x))

for i in 1:4
      ind = findlast(res.value .<= qs[i])
CI = [quantile(res.ccdf[2001:end, ind], pv[i]/2), 
       quantile(res.ccdf[2001:end, ind], 1 - pv[i]/2)]
plot!(p, [qs[i] qs[i]]', [CI[1] CI[2]]', color=:red)
end
plot(p)
pv
qs
alpha=0.05

lower = [quantile(res.ccdf[2001:end, i], alpha/2) for i in axes(res.ccdf, 2)]
upper = [quantile(res.ccdf[2001:end, i], 1 - alpha/2) for i in axes(res.ccdf, 2)]

lower_err = Shat .- lower
upper_err = upper .- Shat
p = plot(res.value, Shat, ribbon = (lower_err, upper_err), color=:grey, fillalpha=0.3, legend=false)
plot!(p, x -> ccdf(tdist, x))

png("plot.png")
