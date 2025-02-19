module jointest

using LinearAlgebra
using Distributions
using Random
using LogExpFunctions
using SpecialFunctions
using Statistics
using SparseArrays


function prob_sample(x)
  out = [x[i] > 0.0 ? rand(Random.default_rng(), Gamma(x[i],1)) : 0.0 for i in eachindex(x)]
  out ./= sum(out)
  return out
end

function freq_sample(x)
  hs = [x[i] > 0.0 ? rand(Random.default_rng(), Gamma(x[i],1)) : 0.0 for i in eachindex(x)]
  return cumsum(hs)
end

prob2ccdf(x) = reverse(cumsum(reverse(x)))
colmarginal(x) = cumsum(vec(sum(x,dims=2)))

function asample(le_rank,re_rank,ls_rank,rs_rank,h,lam)
  m = size(lam,1)
  a = spzeros(m,m)
  for i in (le_rank):(re_rank)
    start = max(ls_rank,i+1)
    for j in start:rs_rank
        a[i,j] += h[j-i] * lam[i]
    end
  end
  rind, cind, val = findnz(a)
  val ./= sum(val)
  occ = rand(Random.default_rng(), Categorical(val))
  out = zero(a)
  out[rind[occ], cind[occ]] = 1.0
  return out
end

function Asample!(a, le_rank,re_rank,ls_rank,rs_rank,h,lam)
  m = size(lam,1)
  A = zero(a)
  for i in eachindex(rs_rank)
    A .+= asample(le_rank[i], re_rank[i], ls_rank[i], rs_rank[i], h, lam)
  end
  copy!(a, A)
end


function bsample(le_rank,re_rank,ls_rank,rs_rank,h,lam)
  m = size(lam,1)
  B = spzeros(m,m)
  rho = zeros(2)
  #convolution
  for i in le_rank:re_rank
    start = max(ls_rank,i+1)
    for j in start:rs_rank
        rho[1] += lam[i]*sum(h[1:(m-(rs_rank-(j-i)))])
        rho[2] += lam[i]*sum(h[(m-(rs_rank-(j-i))+1):end])
    end
  end
  rho ./= sum(rho)
  #expectation
  for i in le_rank:re_rank
    for j in (rs_rank+1):m
        B[i,j] += h[j-i]
    end
  end
  denb = sum(B)
  if denb > 0
    b = 0
    if zero(rho[1]) < rho[1] < one(rho[1])
      b = rand(Random.default_rng(), Geometric(rho[1]))
    end
    rind, cind, val = findnz(B)
    val ./= denb
    y = rand(Random.default_rng(), Multinomial(b, val))
    copy!(B, sparse(rind, cind, y, m, m))
  end
  return B
end

function Bsample!(b,le_rank,re_rank,ls_rank,rs_rank,h,lam)
  m = size(h,1)
  B = zero(b)
  for i in eachindex(re_rank)
    B .+= bsample(le_rank[i], re_rank[i], ls_rank[i], rs_rank[i], h, lam)
  end
  copy!(b, B)
end

function summingup!(α, β, A, B, h, lam)
  lp = 0.
  m = size(A,1)
  alpha = zero(α)
  beta = zero(β)
  for i in 1:(m-1)
    for j in (i+1):m
      alpha[j-i] += A[i,j] + B[i,j]
      beta[i] += A[i,j] + B[i,j]
      if h[j-i]>0
      lp += xlogy(A[i,j]+B[i,j],h[j-i])
      end
      if lam[i]>0
      lp += xlogy(A[i,j]+B[i,j],lam[i])
      end
    end
  end
  copy!(α,alpha)
  copy!(β,beta)
  return lp
end

function gibbssampler(LE, RE, LS, RS, Tmax, iter)
  ti = sort(unique([LE;RE;LS;RS;Tmax]))
  ti = ti[isfinite.(ti)]
  #n = size(LS,1)
  le_rank = indexin(LE, ti)
  re_rank = indexin(RE, ti)
  ls_rank = indexin(LS, ti)
  rs_rank = indexin(RS, ti)
  N = length(LE)
  m = size(ti,1)
  ccdf = zeros(iter, m)
  intensity = zeros(iter, m)
  alpha = ones(m)
  beta = ones(m)
  A = spzeros(m, m)
  B = spzeros(m, m)
  bt = zeros(iter, m)
  sumB = zeros(iter)
  logprob = zeros(iter)
  for i in 1:iter
    h = prob_sample(alpha)
    lam = prob_sample(beta)
    Asample!(A, le_rank, re_rank, ls_rank, rs_rank, h, lam)
    Bsample!(B, le_rank, re_rank, ls_rank, rs_rank, h, lam)
    logprob[i] = summingup!(alpha, beta, A, B, h, lam)
    ccdf[i,:] = prob2ccdf(h)
    sumB[i] = sum(B)
    bt[i,:] = colmarginal(B)
    intensity[i,:] = freq_sample(beta)
  end
  return (value = ti, ccdf = ccdf, intensity = intensity,  b = sumB, bt = bt, logprob = logprob)
end

end