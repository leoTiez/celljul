using LinearAlgebra 
using Plots
using SpecialFunctions

const min_repair_t = 3.0

# Define class like nucleation structure
Base.@kwdef mutable struct Nucleation
    x::Int32
    y::Int32
    timestamp::Float32
    m::Float32
    d_vol::Float32
    volume::Float32=0.
end

function get_coordinates(self::Nucleation)::Vector{Float32}
    return [convert(Float32, self.x), convert(Float32, self.y)]
end

function update_radius!(self::Nucleation)::Nothing
    self.volume += self.d_vol
    return nothing
end

function do_dissociate(self::Nucleation, current_time::Float32, decay::Float32=.5)::Bool
    if current_time - self.timestamp > self.min_repair_t &
         rand() > exp(-decay * (current_time - self.timestamp))
        return true
    else
        return false
    end
end

function calc_circular_radius(self::Nucleation)::Float32
    # calculate volume fraction that changes w/ the dimension of the space
    factor = pi.^(self.m / 2.) / (gamma((self.m / 2.) + 1.) / pi)
    return (self.volume / factor)^(1. / self.m)
end

function circle!(
    state::Array{Float32, 2},
    time_state::Array{Float32, 2},
    centre::Vector{Vector{Float32}},
    radius::Vector{Float32},
    t::Float32
)::Nothing
    coordinates = findall(time_state .== -1)
    # coordinates = findall(iszero, state)
    @inbounds for coords in coordinates
        x = coords[1]
        y = coords[2]
        @inbounds for (c, r) in zip(centre, radius)
            if norm(c - [x, y]) <= r
                state[x, y] = 1
                time_state[x, y] = t
                break
            end
        end
    end
    return nothing
end

function simulation_step!(
    t::Float32,
    state::Array{Float32, 2},
    time_state::Array{Float32, 2},
    nucleation_list::Vector{Nucleation},
    m::Float32,
    n_p::Float64,
    d_vol::Float32,
    sim_protein::Bool=true,
    verbosity::Int32=1
)::Nothing
    # Growth
    update_radius!.(nucleation_list)
    circ_radius = calc_circular_radius.(nucleation_list)
    positions = get_coordinates.(nucleation_list)
    if verbosity > 1
        @time circle!(state, time_state, positions, circ_radius, t)
    else
        circle!(state, time_state, positions, circ_radius, t)
    end
    
    if sim_protein
        mask = ((t .- time_state) .> (min_repair_t + rand(Exponential(1.)))) .& (time_state != -1.)
        state[mask] .= 0
    end

    new_nucleation = findall((rand(Binomial(1, n_p), size(state)) .== 1.) .& (time_state .== -1.))
    @inbounds for pos in new_nucleation
        time_state[pos] = t
        state[pos] = 1.
        push!(
            nucleation_list, 
            Nucleation(
                x=pos[1],
                y=pos[2],
                timestamp=convert(Float64, t), 
                m=m,
                d_vol=d_vol
            )
        )
    end
    return nothing
end

function run_simulation(
    m::Float64,
    beta::Float64,
    theta::Float64;
    sim_protein::Bool=true,
    to_time::Float64=120.,
    dims::Tuple{Int64, Int64}=(100, 100),
    n_p=nothing,
    g_p=nothing,
    sig_p=nothing,
    verbosity::Int=2,
    to_gif::Bool=false
)::Vector{Float64}

    if n_p == nothing
        n_p = beta^m
    end

    if sig_p == nothing || g_p == nothing
        d_vol = 1.
    else
        d_vol = sig_p * g_p^m
    end

    state = zeros(Float32, dims)
    # Conversion times
    time_state = -ones(Float32, dims)

    # calculate scaling
    do_sort = false
    nucleation_list = Vector{Nucleation}()
    sim_rep_frac = Vector{Float32}()
    if verbosity > 0
        fig = heatmap(state, title="Time 0 min")
    end
    
    m, d_vol, theta = Float32(m), Float32(d_vol), Float32(theta)
    verbosity = Int32(verbosity)
    if to_gif
        anim = @animate for t in 1.:Float32(to_time)
            simulation_step!(
                t,
                state,
                time_state, 
                nucleation_list, 
                m, 
                n_p,
                d_vol,
                sim_protein,
                verbosity
            )
            push!(sim_rep_frac, theta * sum(state) / convert(Float64, length(state)))
            title_string = "Time $(round(Int, t)) min"
            display(heatmap!(fig, state, title=title_string))
        end
    else
        @inbounds for t in 1:Float32(to_time)
            simulation_step!(
                t,
                state,
                time_state, 
                nucleation_list, 
                m, 
                n_p,
                d_vol,
                sim_protein,
                verbosity
            )
            push!(sim_rep_frac, theta * sum(state) / convert(Float64, length(state)))
            if verbosity > 0
                title_string = "Time $(round(Int, t)) min"
                display(heatmap!(fig, state, title=title_string))
            end
        end
    end
    if to_gif
        anim_type = sim_protein ? "protein" : "repair"
        save_str = "figures/gif/simulation_$(anim_type)_$(size(state)).gif"
        gif(anim, save_str, fps = 10)
    end
    return sim_rep_frac
end