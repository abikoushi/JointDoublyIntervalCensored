#=
Section 3 "Simulation Study"
Fig 3 - 6 
=#
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
# using QuadGK
using HCubature
using CSV

include("./JointTools_Gibbs.jl")
include("./JointTools.jl")

Lambda(t, shape) = t^shape
lambda(t, shape) = shape * (t^(shape-1.0))

function integrand(X, d, tau, WS, sigma)
    Distributions.ccdf(d, tau - X[1]*WS - X[2]) * lambda(X[2], sigma)
end

tdist1 = LogNormal(log(10) - 0.5, 1)
tdist2 = MixtureModel([Weibull(2, 5/gamma(1+1/2)), Weibull(0.5, 10)], [2/3,1/3])
tdists = [tdist1, tdist2]

#quantile(tdist1, [0.5, 0.8, 0.9, 0.95])
#round.(quantile(tdist1, [0.5, 0.8, 0.9, 0.95]), digits=1)
#quantile(tdist2, [0.5, 0.8, 0.9, 0.95])
#round.(quantile(tdist2, [0.5, 0.8, 0.9, 0.95]), digits=1)

function simfunc(iter, sigma, tau, WE, WS, tdist)
    trueB, _ = hcubature(x -> integrand(x, tdist, tau, WS, sigma), [0,0], [1,tau])
    prob = [0.5, 0.8, 0.9, 0.95] #4
    qs = quantile(tdist, prob)
    EN = tau ^ sigma #E[N]
    tvalue = 1.0 .- prob
    qs_intensity = prob*tau #4
    tvalue_intensity = qs_intensity .^ sigma
    pv_B_df = Vector{Any}(undef, iter)
    pv_inc_prob_df = Vector{Any}(undef, iter)
    pv_intensity_df = Vector{Any}(undef, iter)
    prop = Vector{Any}(undef, iter)
    prop_vb = Vector{Any}(undef, iter)
    prop_em = Vector{Any}(undef, iter)
    for j in 1:iter
        d = jointtool.makeDIC(Random.default_rng(j), tdist, tau, WE, WS, sigma)
        f = d[4] .<= tau
        LE = floor.(d[1][f])
        RE = ceil.(d[2][f])
        LS = floor.(d[3][f])
        RS = ceil.(d[4][f])

        res_em = jointtool.EM(LE, RE, LS, RS, tau)
        df = DataFrame(value = res_em.value, ccdf = jointtool.prob2ccdf(res_em[2]/sum(res_em[2])), intensity = cumsum(res_em[3]), id = j, dist = string(nameof(typeof(tdist))),  sigma=sigma, tau=tau, WE=WE, WS=WS)
        prop_em[j] = df

        res_vb = jointtool.vb(LE, RE, LS, RS, tau)
        df = DataFrame(value = res_em.value, ccdf = jointtool.prob2ccdf(res_em[2]/sum(res_em[2])), intensity = cumsum(res_em[3]), id = j, dist = string(nameof(typeof(tdist))),  sigma=sigma, tau=tau, WE=WE, WS=WS)
        prop_vb[j] = df

        res = jointest.gibbssampler(LE, RE, LS, RS, tau, 4000)
        Shat = vec(mean(res.ccdf[2001:end,:], dims=1))
        Lhat = vec(mean(res.intensity[2001:end,:], dims=1))
        df = DataFrame(value = res.value, ccdf = Shat, intensity = Lhat, id = j, dist = string(nameof(typeof(tdist))),  sigma=sigma, tau=tau, WE=WE, WS=WS)
        prop[j] = df

        pv_inc_prob_j = jointtool.simulated_pvalue(res.ccdf[2001:end,:], res.value, qs, tvalue)
        pv_intensity_j = jointtool.simulated_pvalue(res.intensity[2001:end,:], res.value, qs_intensity, tvalue_intensity)
        pv_B_j = jointtool.simulated_pvalue(res.b[2001:end], trueB)
        pv_B_df[j] = DataFrame(pv=pv_B_j, id=j, dist = string(nameof(typeof(tdist))), sigma=sigma, tau=tau, WE=WE, WS=WS)
        pv_inc_prob_df[j] = DataFrame(pv=pv_inc_prob_j, pos=prob, id=j, dist = string(nameof(typeof(tdist))), sigma=sigma, tau=tau, WE=WE, WS=WS)
        pv_intensity_df[j] = DataFrame(pv=pv_intensity_j, pos=prob, id=j,dist = string(nameof(typeof(tdist))), sigma=sigma, tau=tau, WE=WE, WS=WS)
   end
   prop = vcat(prop...)
   prop_em = vcat(prop_em...)
   prop_vb = vcat(prop_vb...)
   pv_B_df = vcat(pv_B_df...)
   pv_inc_prob_df = vcat(pv_inc_prob_df...)
   pv_intensity_df = vcat( pv_intensity_df...)
   return prop, prop_em, prop_vb, pv_inc_prob_df, pv_intensity_df, pv_B_df
end


function kicksim()
iter = 100
sigmas = [1, 1.2, 0.8]
taus = [50, 100, 200]
WE = 2
WSs = [2, 5, 10]

counter = 1
len = length(tdists)*length(sigmas)*length(taus)*length(WSs)
out_gibbs = Vector{Any}(undef, len)
out_em = Vector{Any}(undef, len)
out_vb = Vector{Any}(undef, len)

out_pv_b = Vector{Any}(undef, len)
out_pv_ccdf = Vector{Any}(undef, len)
out_pv_int = Vector{Any}(undef, len)

for l in eachindex(tdists)
for i in eachindex(sigmas)
    for j in eachindex(taus)
        for k in eachindex(WSs)
            prop, prop_em, prop_vb, pv_inc_prob, pv_intensity, pv_B = simfunc(iter, sigmas[i], taus[j], WE, WSs[k], tdists[l])
            out_gibbs[counter] = prop
            out_em[counter] = prop_em
            out_vb[counter] = prop_vb
            out_pv_ccdf[counter] = pv_inc_prob
            out_pv_int[counter] = pv_intensity
            out_pv_b[counter] = pv_B
            println(l,i,j,k)
            counter += 1
        end
    end
end 
end

out_gibbs = vcat(out_gibbs...)
out_em = vcat(out_em...)
out_vb = vcat(out_vb...)
out_pv_ccdf = vcat(out_pv_ccdf...)
out_pv_b = vcat(out_pv_b...)
out_pv_int = vcat(out_pv_int...)

CSV.write("outsim_trunc_gibbs.csv", out_gibbs)
CSV.write("outsim_trunc_em.csv", out_em)
CSV.write("outsim_trunc_vb.csv", out_vb)
CSV.write("outsim_trunc_pv_b.csv", out_pv_b)
CSV.write("outsim_trunc_pv_ccdf.csv", out_pv_ccdf)
CSV.write("outsim_trunc_pv_int.csv", out_pv_int)
end

@time kicksim()

#153947.154873 sec

###

# tdist = LogNormal(log(10)-0.5, 1)
# tau=100.0
# sigma=1.0
# WE=2.0
# WS=2.0
# d = jointtool.makeDIC(Random.default_rng(), tdist, tau, WE, WS, sigma)
# f = d[4] .<= tau
# LE = floor.(d[1][f])
# RE = ceil.(d[2][f])
# LS = floor.(d[3][f])
# RS = ceil.(d[4][f])
# res = jointest.gibbssampler(LE, RE, LS, RS, tau, 4000)

# plot(res.logprob, legend=false)
# Shat = vec(mean(res.ccdf[2001:end,:], dims=1))
# Lhat = vec(mean(res.intensity[2001:end,:], dims=1))

# prob = [0.5, 0.8, 0.9, 0.95] #4
# qs = quantile(tdist, prob)
# pv = jointtool.simulated_pvalue(res.ccdf[2001:end,:], res.value, qs, 1 .- prob)

# i=4
# alpha=pv[i]
# #alpha=0.7
# ind = findlast(res.value .< qs[i])
# CI = [quantile(res.ccdf[2001:end, ind], alpha/2), 
#       quantile(res.ccdf[2001:end, ind], 1 - alpha/2)]
# med = median(res.ccdf[2001:end, ind])
# plot(res.value, Shat, legend=false)
# plot!(x -> ccdf(tdist,x))
# plot!([qs[i] qs[i]]', [CI[1] CI[2]]')
