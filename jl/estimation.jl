using Turing, MCMCChains
using Distributions
using Random


@model function nucleation_prob(
    y::Vector{Float64},
    time_points::Vector{Float64}, 
    m::Float64, 
    theta::Float64, 
    exp_decay::Float64=10., 
    norm_sig::Float64=5.
)
    # priors for nucleation rate n and growth speed g
    n_p ~ truncated(Exponential(exp_decay), 0., 1.)
    shape_p ~ truncated(Normal(pi, norm_sig), 0., Inf)
    
    # computaion of beta and probbility p
    beta = (shape_p .* n_p).^(1 ./ m)
    p = repair_fraction(time_points, m, beta, theta)

    # The number of observations and probabilities
    for i in 1:length(y)
        y[i] ~ Bernoulli(p[i])
    end
end

function mcmc_sample(
    observation::Vector{Float64},
    t_time::Vector{Float64}, 
    m::Float64, 
    theta::Float64;
    exp_decay::Float64=1., 
    norm_sig::Float64=1., 
    iterations::Int=1000,
    show_progress::Bool=false
)
    # Create Gibbs sampler for each variable
    alg = Gibbs(MH(:n_p), MH(:shape_p))
    # Start sampling.
    chain = sample(
        nucleation_prob(
            observation, 
            t_time,
            m,
            theta,
            exp_decay,
            norm_sig
        ),
        alg,
        iterations,
        progress=show_progress,
    )
    return chain
end

function fetch_data(chain::MCMCChains.Chains, v_type::String="argmin_n")::Tuple{Float64, Float64, Float64}
    n_p, g_p, sig_p = nothing, nothing, nothing
    if lowercase(v_type) == "mean"
        n_p = mean(chain[:n_p])
        g_p = 1.
        sig_p = mean(chain[:shape_p])
    elseif lowercase(v_type) == "argmin_n"
        amin = argmin(chain[:n_p])
        n_p = chain[:n_p][amin]
        g_p = 1.
        sig_p = chain[:shape_p][amin]
    elseif lowercase(v_type) == "argmax_g"
        amax = argmax(chain[:g_p])
        n_p = chain[:n_p][amax]
        g_p = 1.
        sig_p = chain[:shape_p][amax]
    elseif lowercase(v_type) == "argmin_sig"
        amin = argmin(chain[:shape_p])
        n_p = chain[:n_p][amin]
        g_p = 1.
        sig_p = chain[:shape_p][amin]
    else
        throw(ArgumentError("Passed v_type %s is not accepted. 
        Try mean|argmin_n|argmin_g|argmin_sig" % v_type))
    end
    return n_p, g_p, sig_p
end
        