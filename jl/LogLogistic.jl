module survdist

using Random
using Distributions
using HypergeometricFunctions

struct LogLogistic{Ta, Tb} <: ContinuousUnivariateDistribution
    a::Ta
    b::Tb
end

#### Parameters

shape(d::LogLogistic) = d.a
scale(d::LogLogistic) = d.b

params(d::LogLogistic) = (d.a, d.b)

#### Statistics

function mean(d::LogLogistic)
    a, b = params(d)
    if a>one(a)
        return (b*pi/a)/sin(pi/a)
    else
        return missing
    end
end

function var(d::Logistic)
    a, b = params(d)
    if a > 2.0
        return b^2.0 * (2*a/sin(2*a) - a^2.0 / sin(a)^2)
    else
        return missing
    end
end

function std(d::LogLogistic)
 return sqrt(var(d))
end


median(d::Logistic) = d.b

function pdf(d::LogLogistic, x)
    a, b = params(d)
    return ((a/b)*(x/b )^(a-1))/((1+(x/b)^a))^2
end

function logpdf(d::LogLogistic, x)
    a, b = params(d)
    return log(a)-log(b)+(a-1)*log(x/b) - 2*log1p((x/b)^a)
end

function cdf(d::LogLogistic, x)
    a, b = params(d)
    return (x^a) / (b^a +x^a)
end

function ccdf(d::LogLogistic, x)
    a, b = params(d)
    return (b^a) / (b^a +x^a)
end

function eqcdf(d::LogLogistic, x::Real)
    out = zero(x)
    if x >= zero(x)
        a, b = params(d)
        out = x*_₂F₁(1, 1/a, 1 + 1/a, -(x/b)^a)*sin(pi/a)/(b*pi/a)
    end
    return out
end

function quantile(d::LogLogistic,p)
    a, b = params(d)
    return b * ((p)/(1-p))^inv(a)
end

function rand(rng::AbstractRNG, d::LogLogistic)
    p = rand(rng)
    return quantile(d, p)
end

end