# # [One-dimensional column example](@id OneD_column)
# In this example we setup a simple 1D column with the [LOBSTER](@ref LOBSTER) biogeochemical model and observe its evolution. The example demonstrates:
# - How to setup OceanBioME's biogeochemical models
# - How to visualise results
# This is forced by idealised mixing layer depth and surface photosynthetically available radiation (PAR) which are setup first.

# ## Install dependencies
# First we check we have the dependencies installed
# ```julia
# using Pkg
# pkg"add OceanBioME, Oceananigans, CairoMakie"
# ```

# ## Model setup
# We load the packages and choose the default LOBSTER parameter set
using OceanBioME, Oceananigans, Printf
using OceanBioME: Boundaries, GasExchange
using OceanBioME.Boundaries.Sediments: sinking_flux
using OceanBioME.SLatissimaModel: SLatissima
using Oceananigans.Fields: FunctionField, ConstantField
using Oceananigans.Units
using JLD2

import Oceananigans.Biogeochemistry: update_tendencies!
import Oceananigans.Biogeochemistry: required_biogeochemical_tracers,
                                     required_biogeochemical_auxiliary_fields,
                                     biogeochemical_drift_velocity

const year = years = 365days
nothing #hide

# ## Surface PAR and turbulent vertical diffusivity based on idealised mixed layer depth 
# Setting up idealised functions for PAR and diffusivity (details here can be ignored but these are typical of the North Atlantic)

@inline PAR⁰(x, y, t) = 60 * (1 - cos((t + 15days) * 2π / year)) * (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2
@inline H(t, t₀, t₁) = ifelse(t₀ < t < t₁, 1.0, 0.0)
@inline fmld1(t) = H(t, 50days, year) * (1 / (1 + exp(-(t - 100days) / 5days))) * (1 / (1 + exp((t - 330days) / 25days)))
@inline MLD(t) = - (10 + 340 * (1 - fmld1(year - eps(year)) * exp(-mod(t, year) / 25days) - fmld1(mod(t, year))))
@inline κₜ(x, y, z, t) = 1e-3 * (1 + tanh((z - MLD(t)) / 10)) / 2 + 0.5e-4    ### 1e-2 * (1 + tanh((z - MLD(t)) / 10)) / 2 + 1e-4
@inline temp(x, y, z, t) = 2.4 * cos(t * 2π / year + 50days) * (0.5 - 0.5 * tanh(0.25 * (abs(z)- 20)))  + 10
@inline salt(x, y, z, t) = (2.4 * cos(t * 2π / year + 50days)) * (0.5 - 0.5 * tanh(0.25 * (abs(z)- 20)))  + 33

nothing #hide

# ## Grid
# Define the grid.
depth_extent=40meters
grid = RectilinearGrid(size = (1, 1, 10), extent = (20meters, 20meters, depth_extent))

# ## Model
# First we define the biogeochemical model including carbonate chemistry (for which we also define temperature (``T``) and salinity (``S``) fields)
# and scaling of negative tracers(see discussion in the [positivity preservation](@ref pos-preservation))
# and then setup the Oceananigans model with the boundary condition for the DIC based on the air-sea CO₂ flux.

#¤ biogeochemistry = LOBSTER(; grid,
#¤                            surface_photosynthetically_active_radiation = PAR⁰,
#¤                            carbonates = true,
#¤                            scale_negatives = true)

biogeochemistry = OXYDEP(; grid, 
                          surface_photosynthetically_active_radiation = PAR⁰,
                          particles = nothing)

clock = Clock(; time = 0.0)
T = FunctionField{Center, Center, Center}(temp, grid; clock)
S = FunctionField{Center, Center, Center}(salt, grid; clock)
#---------------------------
# B O U N D A R Y   C O N D
#---------------------------
# for the bottom boundary
O2_suboxic = 30.0 
Trel = 10000. # relaxation time, s
b_ox = 15.0 #15.0    # dufference of DO conc. across SWI, uM
b_NUT = 20.
b_DOM_ox = 5.0
b_DOM_anox =20.0
@inline F_ox(conc,threshold)    = (0.5+0.5*tanh(conc-threshold))
@inline F_subox(conc,threshold) = (0.5-0.5*tanh(conc-threshold))
#---OXY----------------------
OXY_top = GasExchange(; gas = :OXY)
@inline OXY_bottom_cond(i, j, grid, clock, fields) =  
    @inbounds - (F_ox(fields.OXY[i, j, 1], O2_suboxic) * b_ox + F_subox(fields.OXY[i, j, 1], O2_suboxic) * (0.0 - fields.OXY[i, j, 1])) / Trel
OXY_bottom = FluxBoundaryCondition(OXY_bottom_cond,  discrete_form = true,)
#---DOM----------------------
@inline DOM_bottom_cond(i, j, grid, clock, fields ) =  
    @inbounds (F_ox(fields.OXY[i, j, 1], O2_suboxic) * (b_DOM_ox - fields.DOM[i, j, 1]) + F_subox(fields.OXY[i, j, 1], O2_suboxic) * 2.0 * (b_DOM_anox - fields.DOM[i, j, 1])) / Trel
DOM_bottom = FluxBoundaryCondition(DOM_bottom_cond, discrete_form = true) #, parameters = (; O2_suboxic, b_DOM_ox, Trel),)
#---NUT----------------------
@inline NUT_bottom_cond(i, j, grid, clock, fields ) =  
    @inbounds (F_ox(fields.OXY[i, j, 1], O2_suboxic) * (b_NUT - fields.NUT[i, j, 1]) + F_subox(fields.OXY[i, j, 1], O2_suboxic) * 2.0 * (b_DOM_anox - fields.NUT[i, j, 1])) / Trel
NUT_bottom = FluxBoundaryCondition(NUT_bottom_cond,  discrete_form = true,) #ValueBoundaryCondition(10.0)

#---------------------------
# Model instantiation
#---------------------------

model = NonhydrostaticModel(; grid,
                              clock,
                              closure = ScalarDiffusivity(ν = κₜ, κ = κₜ),
                              biogeochemistry,
                              boundary_conditions = (NUT = FieldBoundaryConditions(bottom = NUT_bottom ),
                                                     DOM = FieldBoundaryConditions(bottom = DOM_bottom),
                                                    # OXY = FieldBoundaryConditions(top = GasExchange(; gas = :OXY), bottom = OXY_bottom_flux),
                                                     OXY = FieldBoundaryConditions(top = OXY_top, bottom = OXY_bottom ),
                                                     ),
                              auxiliary_fields = (; T, S))
#¤                              boundary_conditions = (DIC = FieldBoundaryConditions(top = CO₂_flux), ),

#¤ set!(model, P = 0.03, Z = 0.03, NO₃ = 4.0, NH₄ = 0.05, DIC = 2239.8, Alk = 2409.0)

set!(model, NUT = 10.0, PHY = 0.01, HET = 0.05, OXY = 150., DOM = 1.)  ## initial conditions

# ## Simulation
# Next we setup a simulation and add some callbacks that:
# - Show the progress of the simulation
# - Store the model and particles output

stoptime = 730 #1095 #730 #1825

simulation = Simulation(model, Δt = 6minutes, stop_time = (stoptime)days) 

progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n",
                                iteration(sim),
                                prettytime(sim),
                                prettytime(sim.Δt),
                                prettytime(sim.run_wall_time))
                                                                  
simulation.callbacks[:progress] = Callback(progress_message, TimeInterval(10days))

filename = "columney_out"

model.tracers
NUT, PHY, HET, POM, DOM, OXY = model.tracers
T = model.auxiliary_fields.T
PAR = model.auxiliary_fields.PAR
S = model.auxiliary_fields.S

simulation.output_writers[:profiles] = JLD2OutputWriter(model, (; NUT, PHY, HET, POM, DOM, OXY, T, S, PAR),
                                                        filename = "$filename.jld2",
                                                        schedule = TimeInterval(1day),
                                                        overwrite_existing = true)
nothing #hide

# ## Run!
# We are ready to run the simulation

#---------------------------------------------------------------------------------------
run!(simulation)
#---------------------------------------------------------------------------------------

@info "Loading saved outputs..."

# ## Load saved output
PHY = FieldTimeSeries("$filename.jld2", "PHY")
NUT = FieldTimeSeries("$filename.jld2", "NUT")
HET = FieldTimeSeries("$filename.jld2", "HET")
POM = FieldTimeSeries("$filename.jld2", "POM")
DOM = FieldTimeSeries("$filename.jld2", "DOM")
OXY = FieldTimeSeries("$filename.jld2", "OXY")
T =   FieldTimeSeries("$filename.jld2", "T")
#S =   FieldTimeSeries("$filename.jld2", "S")
PAR = FieldTimeSeries("$filename.jld2", "PAR")


@info "Saved outputs loaded..."

z = jldopen("$filename.jld2")["grid"]["zᵃᵃᶜ"] #[3:29]  #[1:51] #[3:54]
times = T.times
nothing #hide

# ## Plot
# Finally, we plot!

using CairoMakie

fig = Figure(size = (1500, 1000), fontsize = 20)

#axis_kwargs = (xlabel = "Time (days)", ylabel = "z (m)", limits = ((0, times[end] / days), (-200meters, 0)))
#axis_kwargs = (ylabel = "z (m)", limits = ((0, times[end] / days), (-200meters, 0)))
axis_kwargs = (xlabel = "Time (days)", ylabel = "z (m)",
                limits = ((0, times[end] / days), (-depth_extent, 0)),
                xticks = collect(0:365:stoptime))

axPHY = Axis(fig[1, 3]; title = "PHY, mmolN/m³", axis_kwargs...)
hmPHY = heatmap!(times / days, z, interior(PHY, 1, 1, :, :)', colormap = Reverse(:davos100))
Colorbar(fig[1, 4], hmPHY)

axHET = Axis(fig[2, 3]; title = "HET, mmolN/m³", axis_kwargs...)
hmHET = heatmap!(times / days, z, interior(HET, 1, 1, :, :)', colormap = Reverse(:pink))
Colorbar(fig[2, 4], hmHET)

axPOM = Axis(fig[3, 3]; title = "POM, mmolN/m³", axis_kwargs...)
hmPOM = heatmap!(times / days, z, interior(POM, 1, 1, :, :)', colormap = Reverse(:bilbao25))
Colorbar(fig[3, 4], hmPOM)

axDOM = Axis(fig[3, 1]; title = "DOM, mmolN/m³", axis_kwargs...)
hmDOM = heatmap!(times / days, z, interior(DOM, 1, 1, :, :)', colormap = Reverse(:devon10))
Colorbar(fig[3, 2], hmDOM)

axNUT = Axis(fig[1, 1]; title = "NUT, mmolN/m³", axis_kwargs...)
hmNUT = heatmap!(times / days, z, interior(NUT, 1, 1, :, :)', colormap =  Reverse(:cherry))
Colorbar(fig[1, 2], hmNUT)


axOXY = Axis(fig[2, 1]; title = "OXY, mmol/m³", axis_kwargs...)
hmOXY = heatmap!(times / days, z, interior(OXY, 1, 1, :, :)', colormap = :turbo)
Colorbar(fig[2, 2], hmOXY)

#--------------------------------------

axT = Axis(fig[2, 5]; title = "T, oC", axis_kwargs...)
hmT = heatmap!(times / days, z, interior(T, 1, 1, :, :)', colormap = Reverse(:RdYlBu))
Colorbar(fig[2, 6], hmT)

#axS = Axis(fig[2, 3]; title = "S, psu", axis_kwargs...)
#hmS = heatmap!(times / days, z, interior(S, 1, 1, :, :)', colormap = Reverse(:speed))
#Colorbar(fig[2, 4], hmS)

axPAR = Axis(fig[1, 5]; title = "PAR  μE⋅m-2⋅s-1", axis_kwargs...)
hmPAR = heatmap!(times / days, z, interior(PAR, 1, 1, :, :)', colormap = :grayC10) # :linear_grey_0_100_c0_n256)
Colorbar(fig[1, 6], hmPAR)

fig
