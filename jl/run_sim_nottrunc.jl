using Random
using LinearAlgebra
using Distributions
using Statistics
using Plots
using SpecialFunctions
using DataFrames
using Profile
using XLSX
using JLD
using RCall
using Base.Threads

nthreads()
#'Turnbull' function requires fixed truncate point
include("./JointTools.jl")

tdist = MixtureModel([Weibull(2,5/gamma(1+1/2)), Weibull(0.5, 10)],[2/3,1/3])
# mean(tdist) #10.0
# plot(x->pdf(tdist,x),0,50)
# plot(x->ccdf(tdist,x),0,50)
# tau=100
# WE=5
# sigma=1
function simfunc(iter, sigma, tau, WE, WS, tdist)
    prob = [0.5, 0.8, 0.9, 0.95]
    prop = []
    naive = []
    R"
    library(survival)
    "
    @threads for j in 1:iter
        d = jointest.makeDIC(MersenneTwister(), tdist, tau, WE, WS, sigma)
        LE = floor.(d[1])
        RE = ceil.(d[2])
        LS = floor.(d[3])
        RS = ceil.(d[4])
        res = jointest.vb(LE, RE, LS, RS, Inf)
        #plot(res.value, jointest.prob2ccdf(res.alpha/sum(res.alpha)))
        #plot!(x -> ccdf(tdist,x))
        push!(prop, (res[1], jointest.prob2ccdf(res[2]/sum(res[2])), cumsum(res[3])))
        df0 = DataFrame(L = [(LS[i]-RE[i] >0.) ? (LS[i]-RE[i]) : 0. for i in eachindex(LS)], R = [(RS[i]-LE[i] >0.) ? (RS[i]-LE[i]) : 0. for i in eachindex(LS)])

        R"
        sf = survfit(Surv(L, R, type='interval2')~1, data=$df0)
        res0 = with(sf, data.frame(time=time, surv=surv))
        "
        @rget res0

        #plot(naive.time, naive.surv)
        #plot!(x -> ccdf(tdist,x))
        
        push!(naive, res0)
        # cover = zeros(4,6)
        # cis = jointest.sampleccdfci(res[2], 5000, [0.005,0.025,0.05,0.1,0.9,0.95,0.975,0.995])
        # for i in eachindex(qs)
        #     f = res[1] .> qs[i]'
        #     pos = findfirst(f)
        #     cover[:,i] = cis[pos,1:4] .<=  1 .-prob[i] .<= cis[pos,8:-1:5]
        # end
        # coverage_ccdf += cover
   end
   return prop, naive
end

sigma = 1.0
taus = [50,100,200]
# taus .^ sigmas' #E[N]
WEs = [2, 5, 10]
WS = 2

#prop, naive = simfunc(5, sigmas[1], taus[1], WEs[1], WS, tdist)
# outcoverage_ccdf = zeros(length(sigmas),length(taus),length(WEs),4,6)
# outcoverage_intensity = zeros(length(sigmas),length(taus),length(WEs),4,6)
# outcoverage_B = zero(length(sigmas),length(taus),length(WEs),4)
outprop = []
outnaive = []

for j in eachindex(taus)
    for k in eachindex(WEs)
        prop, naive = simfunc(200, sigma, taus[j], WEs[k], WS, tdist)
        append!(outprop,prop)
        append!(outnaive,naive)
        println((j,k))
     end
end

@save "outsim_nottrunc.jld" outprop, outnaive

