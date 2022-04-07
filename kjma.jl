using Random
using Distributions


function repair_fraction(time, m, beta, theta)
    return (1 .- exp.(-(beta .* time).^m)) .* theta
end

function nucleation_growth(time, n, g, shape, m, theta)
    beta =  ((shape .* n .* g.^(m-1)) ./ m).^(1 ./ m)
    return repair_fraction(time, m, beta, theta)
end

function repair_fraction_over_time(to_time, m, beta, theta)
    if typeof(beta) == Float64
        rf_t = Vector{Float64}(undef, convert(Int, to_time))
    elseif typeof(beta) == Array{Float64, 1}
        rf_t = Array{Float64, 2}(undef, convert(Int, to_time), length(beta))
    else
        throw(ArgumentError("Invalid data type of m"))
    end
        
    @inbounds for t in 1:to_time
        if typeof(rf_t) == Array{Float64, 2}
            rf_t[convert(Int, t), :] = repair_fraction(t, m, beta, theta)
        else
            rf_t[convert(Int, t)] = repair_fraction(t, m, beta, theta)
        end
    end
    return rf_t
end

function nucleation_growth_over_time(to_time, n, g, shape, m, theta)
    beta =  ((shape .* n .* g.^(m-1)) ./ m).^(1 ./ m)
    return repair_fraction_over_time(to_time, m, beta, theta)
end

function create_data(to_time::Float64, params::Vector{Float64}, n_points::Int=100)
    # Repair fraction as probabilities
    # Require the following order per row in params: m, beta, theta
    obs_fraction = repair_fraction_over_time(to_time, params[1], params[2], params[3])

    # sample data points
    observation = Array{Float64, 2}(undef, (convert(Int, to_time), n_points))
    @inbounds for (num, p) in enumerate(obs_fraction)
        observation[num, :] = rand(Binomial(1, minimum([maximum([p, 0]), 1])), n_points)
    end

    observation = collect(Iterators.flatten(transpose(observation)))  # transpose because column based
    t_time = repeat([1: 1: to_time;], inner=n_points)
    return observation, t_time
end