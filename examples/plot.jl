using OceanBioME, Oceananigans, Oceananigans.Units

@info "before writing times.."

#times = FieldTimeSeries("boxey_od.jld2", "PHY").times
times = FieldTimeSeries("columney_out.jld2", "PHY").times

@info "after writing times.."
println(keys(model.fields))
#z = interior(FieldTimeSeries("boxey_od.jld2", "PHY"), 1,1,1,:)
z = interior(FieldTimeSeries("columney_out.jld2", "PHY"), 1,1,1,:)
# timeseries = NamedTuple{keys(model.fields)}(FieldTimeSeries("boxey.jld2", "$field")[1, 1, 1, :] for field in keys(model.fields))
@info "after timeseries.."
# ## And plot
using CairoMakie

fig = Figure(size = (1200, 1200), fontsize = 24)
lines(z)
fig

axs = []
for (name, tracer) in pairs(timeseries)
    idx = (length(axs))
    push!(axs, Axis(fig[floor(Int, idx/2), Int(idx%2)], ylabel = "$name", xlabel = "Year", xticks=(0:10)))
    lines!(axs[end], times / year, tracer, linewidth = 3)
end

fig