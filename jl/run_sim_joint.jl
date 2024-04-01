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
include("./JointTools.jl")

#=

=#


tdist = MixtureModel([Weibull(2,5/gamma(1+1/2)), Weibull(0.5, 10)],[2/3,1/3])



function simfunc(iter, sigma, tau, WE, WS, tdist)
    prob = [0.5, 0.8, 0.9, 0.95]
    qs = quantile(tdist, prob)
    coverage_ccdf = zeros(4,4)
    coverage_intensity = zeros(4,4)
    coverage_B = zeros(4)
    prop = []
    for j in 1:iter
        d = jointest.makeDIC(MersenneTwister(), tdist, tau, WE, WS, sigma)
        f = d[4] .<= tau
        trueB = length(d[1]) - sum(f)
        LE = floor.(d[1][f])
        RE = ceil.(d[2][f])
        LS = floor.(d[3][f])
        RS = ceil.(d[4][f])
        res = jointest.vb(LE, RE, LS, RS, tau)
        CIs = jointest.sampleccdfci(res[2], 5000, [0.005,0.025,0.05,0.1,0.9,0.95,0.975,0.995])
        push!(prop, (res[1], jointest.prob2ccdf(res[2]/sum(res[2])), cumsum(res[3])))
        cover = zeros(4,4)
        for i in eachindex(qs)
            f = res[1] .< qs[i]'
            pos = findlast(f)
            if !isnothing(pos)
            cover[:,i] = CIs[pos,1:4] .<=  1 .-prob[i] .<= CIs[pos,8:-1:5]
            end
        end
        coverage_ccdf += cover
        cover = zeros(4,4)
        CIs = jointest.samplefreqci(res[3], 5000, [0.005,0.025,0.05,0.1,0.9,0.95,0.975,0.995])
        ps = [0.5, 0.8, 0.9, 0.95]*tau
        for i in eachindex(ps)
            f = res[1] .< ps[i]'
            pos = findlast(f)
            if !isnothing(pos)
            cover[:,i] = CIs[pos,1:4] .<=  ps[i]^sigma .<= CIs[pos,8:-1:5]
            end
        end
        coverage_intensity += cover
        bs = jointest.sampleB(LE, RE, LS, RS, tau, res[2], res[3], 5000)
        CIs = quantile(bs, [0.005,0.025,0.05,0.1,0.9,0.95,0.975,0.995])
        cover = CIs[1:4] .<=  trueB .<= CIs[8:-1:5]
        coverage_B += cover
   end
   return prop, coverage_ccdf, coverage_intensity, coverage_B
end


#i,j,k = 1,2,1
function kicksim()
    tdist = MixtureModel([Weibull(2,5/gamma(1+1/2)), Weibull(0.5, 10)],[2/3,1/3])
    sigmas = [1,1.2,0.8]
    taus = [50,100,200]
    #taus .^ sigmas' #E[N]
    WEs = [2, 5, 10]
    WS = 2
    outcoverage_ccdf = zeros(length(sigmas),length(taus),length(WEs),4,4)
    outcoverage_intensity = zeros(length(sigmas),length(taus),length(WEs),4,4)
    outcoverage_B = zeros(length(sigmas),length(taus),length(WEs),4)
    outprop = []
    for i in eachindex(sigmas)
        for j in eachindex(taus)
            for k in eachindex(WEs)
                prop, coverage_ccdf, coverage_intensity, coverage_B = simfunc(200, sigmas[i], taus[j], WEs[k], WS, tdist)
                append!(outprop,prop)
                outcoverage_ccdf[i,j,k,:,:] = coverage_ccdf
                outcoverage_intensity[i,j,k,:,:] = coverage_intensity
                outcoverage_B[i,j,k,:] = coverage_B
                println(i,j,k)
            end
        end
    end
    @save "outsim_trunc.jld" outprop, outcoverage_ccdf, outcoverage_intensity, outcoverage_B
end

kicksim()
