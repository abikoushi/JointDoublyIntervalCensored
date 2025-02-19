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
using StatsBase

include("./JointTools_Gibbs.jl")
include("./JointTools.jl")

Lambda(t, shape) = t^shape
lambda(t, shape) = shape * (t^(shape-1.0))

function integrand(x, d, tau, sigma)
    Distributions.ccdf(d, tau-x) * lambda(x, sigma)
end

tdist1 = LogNormal(log(10)-0.5, 1)
tdist2 = MixtureModel([Weibull(2, 5/gamma(1+1/2)), Weibull(0.5, 10)], [2/3,1/3])
tdists = [tdist1, tdist2]

#plot(res.value, vec(mean(res.intensity[2001:end,:], dims=1)), legend=false)
#plot!(res.value, vec(mean(res.bt[2001:end,:], dims=1)), legend=false)
#plot!(x->Lambda(x,1.1), linestyle=:dash, color="grey")

#p = plot(res.value, Shat, legend=false)
#plot!(p, x -> Distributions.ccdf(tdist1, x))
#for i in axes(res.ccdf[2000:end,:],1)
#    plot!(p, res.value, vec(res.ccdf[i,:]), color = "grey", alpha = 0.01)
#end

#p = plot(res.value, vec(mean(res.intensity[2001:end,:], dims=1)), legend=false)
#for i in axes(res.intensity[2001:end,:],1)
#    plot!(p, res.value, vec(res.intensity[i,:]), color = "grey", alpha = 0.1)
#end
#plot(p)

function simfunc(iter, sigma, tau, WE, WS, tdist)
    trueB, _ = quadgk(x -> integrand(x, tdist, tau, sigma), 0.0, tau)
    prob = [0.5, 0.8, 0.9, 0.95] #4
    qs = quantile(tdist, prob)
    EN = tau ^ sigma #E[N]
    tvalue = 1.0 .- prob
    qs_intensity = [0.5, 0.8, 0.9, 0.95]*tau #4
    tvalue_intensity = qs_intensity .^ sigma
    #pv_inc_freq = Matrix{Union{Missing, Float64}}(undef, 4, iter)
    pv_inc_prob = Matrix{Union{Missing, Float64}}(undef, 4, iter)
    pv_intensity = Matrix{Union{Missing, Float64}}(undef, 4, iter)
    pv_B = zeros(iter)
    prop = Vector{Any}(undef, iter)
    prop_em = Vector{Any}(undef, iter)
    for j in 1:iter
        d = jointtool.makeDIC(MersenneTwister(j), tdist, tau, WE, WS, sigma)
        f = d[4] .<= tau
        LE = floor.(d[1][f])
        RE = ceil.(d[2][f])
        LS = floor.(d[3][f])
        RS = ceil.(d[4][f])

        res_em = jointtool.EM(LE, RE, LS, RS, tau)
        df = DataFrame(value = res_em.value, ccdf = jointtool.prob2ccdf(res_em[2]/sum(res_em[2])), intensity = cumsum(res_em[3]), id = j, dist = string(nameof(typeof(tdist))),  sigma=sigma, tau=tau, WE=WE, WS=WS)
        prop_em[j] = df

        res = jointest.gibbssampler(LE, RE, LS, RS, tau, 4000)
        Shat = vec(mean(res.ccdf[2001:end,:], dims=1))
        Lhat = vec(mean(res.intensity[2001:end,:], dims=1))
        df = DataFrame(value = res.value, ccdf = Shat, intensity = Lhat, id = j, dist = string(nameof(typeof(tdist))),  sigma=sigma, tau=tau, WE=WE, WS=WS)
        prop[j] = df

        pv_inc_prob[:,j] = jointtool.simulated_pvalue(res.ccdf[2001:end,:], res.value, qs, tvalue)
        pv_intensity[:,j] = jointtool.simulated_pvalue(res.intensity[2001:end,:], res.value, qs_intensity, tvalue_intensity)
        pv_B[j] = jointtool.simulated_pvalue(res.b[2001:end], trueB)
   end
   prop = vcat(prop...)
   prop_em = vcat(prop_em...)
   return prop, prop_em, pv_inc_prob, pv_intensity, pv_B
end

#iter = 2
#sigmas = [1,1.2,0.8]
#taus = [50,100,200]
#WE = 2
#WSs = [2, 5, 10]
#@time test = simfunc(iter, sigmas[2], taus[3], WE, WSs[3], tdists[1])
#test[3]
#Lambda(200,1.2)
#sum(ismissing.(test[3][1,:]))
#plot(x -> StatsBase.ecdf( collect(skipmissing(test[3][2,:])) )(x) , xlim=[0,1])
#size(test[3])

#histogram(test[3][4,:])

function kicksim()
iter = 100
sigmas = [1,1.2,0.8]
taus = [50,100,200]
WE = 2
WSs = [2, 5, 10]
#prop, prop_em, pv_inc_freq, pv_inc_prob, pv_intensity, pv_B = simfunc(5, sigmas[1], taus[1], WEs[1], WS, tdists[1])

out_pv_incu_prob = Array{Union{Missing, Float64}}(undef, length(tdists), length(sigmas), length(taus), length(WSs), 4, iter)
#out_pv_incu_freq = Array{Union{Missing, Float64}}(undef, length(tdists), length(sigmas), length(taus), length(WSs), 4, iter)
out_pv_intensity  = Array{Union{Missing, Float64}}(undef, length(tdists), length(sigmas), length(taus), length(WSs), 4, iter) 
out_pv_B = Array{Union{Missing, Float64}}(undef, length(tdists), length(sigmas), length(taus), length(WSs), iter)

counter = 1
len = length(tdists)*length(sigmas)*length(taus)*length(WSs)
out_gibbs = Vector{Any}(undef, len)
out_em = Vector{Any}(undef, len)
for l in eachindex(tdists)
for i in eachindex(sigmas)
    for j in eachindex(taus)
        for k in eachindex(WSs)
            prop, prop_em, pv_inc_prob, pv_intensity, pv_B = simfunc(iter, sigmas[i], taus[j], WE, WSs[k], tdists[l])
            out_gibbs[counter] = prop
            out_em[counter] = prop_em
            out_pv_incu_prob[l, i, j, k, :, :] = pv_inc_prob
            out_pv_intensity[l, i, j, k, :, :] = pv_intensity
            out_pv_B[l, i, j, k, :] = pv_B
            println(l,i,j,k)
            counter += 1
        end
    end
end 
end

out_gibbs = vcat(out_gibbs...)
out_em = vcat(out_em...)

CSV.write("outsim_trunc_gibbs.csv", out_gibbs)
CSV.write("outsim_trunc_em.csv", out_em)
#@save "outsim_pv_incubation_freq.jld" out_pv_incu_freq
@save "outsim_pv_incubation_prob.jld" out_pv_incu_prob
@save "outsim_pv_intensity.jld" out_pv_intensity
@save "outsim_pv_B.jld" out_pv_B
end

@time kicksim()
