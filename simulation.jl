using LinearAlgebra 
using Plots

# Define class like nucleation structure
Base.@kwdef mutable struct Nucleation
    x::Int
    y::Int
    timestamp::Float64
    m::Float64
    g::Float64
    sig_p::Float64
    radius=0.
    min_repair_t=3.
end

function get_coordinates(self::Nucleation)::Vector{Int}
    return [convert(Float64, self.x), convert(Float64, self.y)]
end

function update_radius!(self::Nucleation, scaling::Float64=1.)::nothing
    self.radius += self.g^self.m * scaling
    return nothing
end

function do_dissociate(self::Nucleation, current_time::Float64, decay::Float64=.5)::Bool
    if current_time - self.timestamp > self.min_repair_t &
         rand() > exp(-decay * (current_time - self.timestamp))
        return true
    else
        return false
    end
end

function calc_circular_radius(self::Nucleation)::Float64
    surface = self.sig_p * self.radius
    return sqrt(surface / (4 * pi))
end

function circle!(
    state::Array{Float64, 2},
    coordinates::Array{Int, 2},
    centre::Vector{Int},
    radius::Float64
)::nothing
    coords = coordinates[[norm(c .- centre) for c in coordinates] .<= radius]
    for (x, y) in coords
        state[x, y] = 1
    end
    return nothing
end

function simulation_step!(
    t::Float64,
    state::Array{Float64, 2},
    coordinates::Array{Int, 2},
    nucleation_list::Vector{Nucleation},
    m::Float64,
    n_p::Float64,
    g_p::Float64,
    sig_p::Float64,
    scaling::Float64
)::nothing

    if coordinates == nothing
        xs = 1:size(state)[1]
        ys = 1:size(state)[2]
        coordinates = [convert.(Int, [x, y]) for x in xs for y in ys]
    end
    
    # Growth
    @inbounds for nucl in nucleation_list
        update_radius!(nucl, scaling)
        circ_radius = calc_circular_radius(nucl)
        circle!(state, coordinates, get_coordinates(nucl), circ_radius)
    end
        
    new_nucleation = findall((rand(Binomial(1, n_p), size(state)) .== 1.) .& (state .!= 1.))
    @inbounds for pos in new_nucleation
        push!(
            nucleation_list, 
            Nucleation(
                x=pos[1],
                y=pos[2], 
                timestamp=convert(Float64, t), 
                m=m,
                g=g_p,
                sig_p=sig_p
            )
        )
    end
    return nothing
end

function run_simulation(
    n_p::Float64,
    g_p::Float64,
    sig_p::Float64,
    m::Float64,
    theta::Float64;
    to_time::Float64=120.,
    dim::Tuple{Int64, Int64}=(100, 100)
)::Vector{Float64}

    state = zeros(dim)

    # Create coordinates
    xs = 1:size(state)[1]
    ys = 1:size(state)[2]
    coordinates = [convert.(Int, [x, y]) for x in xs for y in ys]

    # calculate scaling
    scaling = 2. ./ (sig_p .* g_p.^m)
    do_sort = false
    nucleation_list = Vector{Nucleation}
    sim_rep_frac = Vector{Float64}
    fig = heatmap(state)
    for t in 1:to_time
        simulation_step!(t, state, coordinates, nucleation_list, m, n_p, g_p, sig_p, scaling)
        display(heatmap!(fig, state))
        push!(sim_rep_frac, theta * sum(state) / convert(Float64, length(state)))
    end
    return sim_rep_frac
end