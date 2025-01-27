
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
include("./JointTools.jl")


Lambda(t, shape) = t^shape
lambda(t, shape) = shape * (t^(shape-1.0))


function integrand(x, d, tau, sigma)
    ccdf(d, tau-x) * lambda(x, sigma)
end

tdist1 = LogNormal(log(10)-0.5, 1)
τ = 20.0
trueB = quadgk(x -> integrand(x, tdist1, τ, 1.2), 0.0, τ)
b[1]

tdist1 = LogNormal(log(10)-0.5, 1)
tdist2 = MixtureModel([Weibull(2,5/gamma(1+1/2)), Weibull(0.5, 10)],[2/3,1/3])

tdists = [tdist1, tdist2]

mean.(tdists)
std.(tdists)

sigmas = [1,1.2,0.8]
taus = [50,100,200]
WEs = [2, 5, 10]

function simfunc(iter, sigma, tau, WE, WS, tdist)
    trueB = quadgk(x -> integrand(x, tdist, tau, sigma), 0.0, tau)
    prob = [0.5, 0.8, 0.9, 0.95]
    qs = quantile(tdist, prob)
    EN = tau ^ sigma #E[N]
    tvalue = cdf.(tdist, qs)*ENs
    qs_intensity = [0.5, 0.8, 0.9, 0.95]*tau
    tvalue_intensity = qs_intensity .^ sigma
    pv_incubation = zeros(4, iter)
    pv_intensity = zeros(4, iter)
    pv_B = zeros(iter)
    prop = []
    for j in 1:iter
        d = jointest.makeDIC(MersenneTwister(), tdist, tau, WE, WS, sigma)
        f = d[4] .<= tau
        LE = floor.(d[1][f])
        RE = ceil.(d[2][f])
        LS = floor.(d[3][f])
        RS = ceil.(d[4][f])
        res = jointest.vb(LE, RE, LS, RS, tau)
        push!(prop, (res[1], jointest.prob2ccdf(res[2]/sum(res[2])), cumsum(res[3])))
        Rs = mapreduce(permutedims, vcat,[jointest.randfreq(res.alpha) for _ in 1:5000])
        pv_incubation[:,iter] = jointest.simulated_pvalue(Rs, res.value, qs, tvalue)

        Rs = mapreduce(permutedims, vcat,[jointest.randfreq(res.beta) for _ in 1:5000])
        pv_incubation[:,iter] = jointest.simulated_pvalue(Rs, res.value, qs_intensity, tvalue_intensity)
   
        bs = jointest.sampleB(LE, RE, LS, RS, tau, res[2], res[3], 5000)
        pv_B[iter] = jointest.simulated_pvalue(Rs, trueB)
   end
   return prop, pv_incubation, pv_intensity, pv_B
end


#i,j,k = 1,2,1
function kicksim()
    tdist = MixtureModel([Weibull(2,5/gamma(1+1/2)), Weibull(0.5, 10)],[2/3,1/3])
    sigmas = [1,1.2,0.8]
    taus = [50,100,200]
    #taus .^ sigmas' #E[N]
    WEs = [2, 5, 10]
    WS = 2
    out_pv_incubation = zeros(length(sigmas),length(taus),length(WEs),4)
    out_pv_intensity  = zeros(length(sigmas),length(taus),length(WEs),4) 
    out_pv_B = zeros(length(sigmas),length(taus),length(WEs),4)
    out_prop = []
    for i in eachindex(sigmas)
        for j in eachindex(taus)
            for k in eachindex(WEs)
                prop, pv_incubation, pv_intensity, pv_B = simfunc(200, sigmas[i], taus[j], WEs[k], WS, tdist)
                append!(outprop,prop)
                out_pv_incubation[i,j,k,:] = pv_incubation
                out_pv_intensity[i,j,k,:] = pv_intensity
                out_pv_B[i,j,k] = pv_B
                println(i,j,k)
            end
        end
    end
    @save "outsim_trunc.jld" out_prop,out_pv_incubation,out_pv_intensity,out_pv_B
end

kicksim()
