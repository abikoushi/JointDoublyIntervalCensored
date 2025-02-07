#=
Section 3 "Simulation Study"
Fig 3 - 5 
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
using QuadGK
using CSV
include("./JointTools.jl")




Lambda(t, shape) = t^shape
lambda(t, shape) = shape * (t^(shape-1.0))

function integrand(x, d, tau, sigma)
    ccdf(d, tau-x) * lambda(x, sigma)
end

tdist1 = LogNormal(log(10)-0.5, 1)
tdist2 = MixtureModel([Weibull(2, 5/gamma(1+1/2)), Weibull(0.5, 10)],[2/3,1/3])
tdists = [tdist1, tdist2]

# mean.(tdists)
# std.(tdists)
# var.(tdists)

#pl0 = plot(res[1], cumsum(res[2]/sum(res[2])), legend=false, ticks_direction=:out)
#plot!(pl0, x->cdf(tdist,x))
#for i in 1:5000
#    plot!(pl0, res[1], Rs_b[i,:], alpha=0.1, color=:grey)
#end


tdist = LogNormal(log(10)-0.5, 1)
sigma = 1.2
tau = 200
WE = 2
WS = 5



function simfunc(iter, sigma, tau, WE, WS, tdist)
    trueB = quadgk(x -> integrand(x, tdist, tau, sigma), 0.0, tau)
    prob = [0.5, 0.8, 0.9, 0.95] #4
    qs = quantile(tdist, prob)
    EN = tau ^ sigma #E[N]
    tvalue = prob*EN
    qs_intensity = [0.5, 0.8, 0.9, 0.95]*tau #4
    tvalue_intensity = qs_intensity .^ sigma
    pv_inc_freq = zeros(4, iter)
    pv_inc_prob = zeros(4, iter)
    pv_intensity = zeros(4, iter)
    pv_B = zeros(iter)
    prop = Vector{Any}(undef, iter)
    prop_em = Vector{Any}(undef, iter)
    for j in 1:iter
        d = jointest.makeDIC(MersenneTwister(j), tdist, tau, WE, WS, sigma)
        f = d[4] .<= tau
        LE = floor.(d[1][f])
        RE = ceil.(d[2][f])
        LS = floor.(d[3][f])
        RS = ceil.(d[4][f])
        res = jointest.vb(LE, RE, LS, RS, tau)
        df = DataFrame(value = res.value, ccdf = jointest.prob2ccdf(res[2]/sum(res[2])), intensity = cumsum(res[3]), id = j, dist = string(nameof(typeof(tdist))),  sigma=sigma, tau=tau, WE=WE, WS=WS)
        prop[j] = df

        res_em = jointest.EM(LE, RE, LS, RS, tau)
        df = DataFrame(value = res_em.value, ccdf = jointest.prob2ccdf(res_em[2]/sum(res_em[2])), intensity = cumsum(res_em[3]), id = j, dist = string(nameof(typeof(tdist))),  sigma=sigma, tau=tau, WE=WE, WS=WS)
        prop_em[j] = df

        Rs = jointest.randfreqn(res.alpha, 5000)
        pv_inc_freq[:,j] = jointest.simulated_pvalue(Rs, res.value, qs, tvalue)

        Rs ./= Rs[:,length(res.alpha)]
        pv_inc_prob[:,j] = jointest.simulated_pvalue(Rs, res.value, qs, prob)

        Rs = jointest.randfreqn(res.beta, 5000)
        pv_intensity[:,j] = jointest.simulated_pvalue(Rs, res.value, qs_intensity, tvalue_intensity)

        bs = jointest.sampleB(LE, RE, LS, RS, tau, res[2], res[3], 5000)
        pv_B[j] = jointest.simulated_pvalue(bs, trueB[1])
   end
   prop = vcat(prop...)
   prop_em = vcat(prop_em...)
   return prop, prop_em, pv_inc_freq, pv_inc_prob, pv_intensity, pv_B
end

test = simfunc(2, 1, 50, 2, 2, tdists[1])

# kick simfunc
iter = 200
sigmas = [1,1.2,0.8]
taus = [50,100,200]
WE = 2
WSs = [2, 5, 10]
#prop, prop_em, pv_inc_freq, pv_inc_prob, pv_intensity, pv_B = simfunc(5, sigmas[1], taus[1], WEs[1], WS, tdists[1])

out_pv_incubation = zeros(length(tdists), length(sigmas), length(taus), length(WSs), 4, iter)
out_pv_intensity  = zeros(length(tdists), length(sigmas), length(taus), length(WSs), 4, iter) 
out_pv_B = zeros(length(tdists), length(sigmas), length(taus), length(WSs), iter)

counter = 1
len = length(tdists)*length(sigmas)*length(taus)*length(WSs)
out_vb = Vector{Any}(undef, len)
out_em = Vector{Any}(undef, len)
for l in eachindex(tdists)
for i in eachindex(sigmas)
    for j in eachindex(taus)
        for k in eachindex(WSs)
            prop, prop_em, pv_inc_prob, pv_inc_freq, pv_intensity, pv_B = simfunc(iter, sigmas[i], taus[j], WE, WSs[k], tdists[l])
            out_vb[counter] = prop
            out_em[counter] = prop_em
            out_pv_incubation[l, i, j, k, :, :] = pv_incubation
            out_pv_intensity[l, i, j, k, :, :] = pv_intensity
            out_pv_B[l, i, j, k, :] = pv_B
            println(l,i,j,k)
            counter += 1
        end
    end
end 
end
    


CSV.write("outsim_trunc_vb.csv", out_vb)
CSV.write("outsim_trunc_em.csv", out_em)

@save "outsim_pv_incubation.jld" out_pv_incubation
@save "outsim_pv_intensity.jld" out_pv_intensity
@save "outsim_pv_B.jld" out_pv_B

