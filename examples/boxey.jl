# # [Box model](@id box_example)
# In this example we setup a [LOBSTER](@ref LOBSTER) biogeochemical model in a single box configuration.
# This example demonstrates:
# - How to setup OceanBioME's biogeochemical models as a stand-alone box model

# ## Install dependencies
# First we check we have the dependencies installed
# ```julia
# using Pkg
# pkg"add OceanBioME"
# ```

# ## Model setup
# Load the packages and setup the initial and forcing conditions
using OceanBioME, Oceananigans, Oceananigans.Units
using Oceananigans.Fields: FunctionField, ConstantField

const year = years = 365day
nothing #hide

# This is forced by a prescribed time-dependent photosynthetically available radiation (PAR)
PAR⁰(t) = 60 * (1 - cos((t + 15days) * 2π / year)) * (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2

z = -10 # specify the nominal depth of the box for the PAR profile
PAR(t) = PAR⁰(t) * exp(0.2z) # Modify the PAR based on the nominal depth and exponential decay

temp(t) = 15 # 2.4 * cos(t * 2π / year + 50days) + 10

nothing #hide

@info "preparing the model..."
#Δt = 15minutes
#stop_time = 6years
#last_Δt = 15minutes

#¤clock = Clock(; time = 0.0)
#¤ clock = Clock(; time = 0.0, last_Δt = Inf, iteration = 0, stage = 1)
#clock = Clock(; time = 0.0, last_Δt = Inf, iteration = 0)
# Set up the model. Here, first specify the biogeochemical model, followed by initial conditions and the start and end times
#model = BoxModel(biogeochemistry = LOBSTER(grid = BoxModelGrid, light_attenuation_model = nothing), forcing = (; PAR))
#model = BoxModel(biogeochemistry = NPZD(grid = BoxModelGrid, light_attenuation_model = nothing), forcing = (; PAR))
#model = BoxModel(biogeochemistry = NutrientPhytoplanktonZooplanktonDetritus(grid = BoxModelGrid, light_attenuation_model = nothing), forcing = (; PAR))
# !!!!!!!!!!!!!!!!!!!!!!!!!!! to recompile start from here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


model = BoxModel(biogeochemistry = OXYDEP(grid = BoxModelGrid, light_attenuation_model = nothing), forcing = (; PAR, temp))


#println(OXYDEP)
#set!(model, NO₃ = 5.0, NH₄ = 0.1, P = 0.01, Z = 0.01)
#set!(model, N = 7.0, P = 0.01, Z = 0.05)
set!(model, NUT = 10.0, PHY = 0.01, HET = 0.05, OXY = 150., DOM = 1.)  ## initial conditions??

simulation = Simulation(model; Δt = 15minutes, stop_time = 6years) #6

println(model)
simulation.output_writers[:fields] = JLD2OutputWriter(model, model.fields; filename = "boxey_od1.jld2", schedule = TimeInterval(10days), overwrite_existing = true)

prog(sim) = @info "$(prettytime(time(sim))) in $(prettytime(simulation.run_wall_time))"

#simulation.callbacks[:progress] = Callback(prog, IterationInterval(1000000))
simulation.callbacks[:progress] = Callback(prog, IterationInterval(9600)) # 1/100d for dt=15 min

# ## Run the model (should only take a few seconds)
@info "Running the model..."

run!(simulation)

# ## Load the output

times = FieldTimeSeries("boxey_od1.jld2", "PHY").times

timeseries = NamedTuple{keys(model.fields)}(FieldTimeSeries("boxey_od1.jld2", "$field")[1, 1, 1, :] for field in keys(model.fields))

# ## And plot
using CairoMakie

fig = Figure(size = (1200, 1200), fontsize = 24)

axs = []
for (name, tracer) in pairs(timeseries)
    idx = (length(axs))
    push!(axs, Axis(fig[floor(Int, idx/2), Int(idx%2)], ylabel = "$name", xlabel = "Year", xticks=(0:10)))
    lines!(axs[end], times / year, tracer, linewidth = 3)
end

fig
