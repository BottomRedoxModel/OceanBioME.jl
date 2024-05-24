using OceanBioME, Oceananigans, Oceananigans.Units
# using Plots

# # Define the function
# my_function(x) = 1 - 0.5 * (1 + tanh(x - 30))

# # Generate x values
# x_values = range(0, stop=60, length=1000)  # Adjust the range as needed

# # Generate y values
# y_values = my_function.(x_values)

# # Plot the function
# plot(x_values, y_values, label="y = (1 - 0.5 * (1 + tanh(x - 30)))", xlabel="x", ylabel="y", title="Plot of y = (1 - 0.5 * (1 + tanh(x - 30)))")
"""
my_function(x) = 0.5 - 0.5 * tanh(0.25 * (abs(x)- 30))
x_values = range(0, stop=60, length=1000)  # Adjust the range as needed
y_values = my_function.(x_values)
plot(x_values, y_values, label="y = 0.5 - 0.5 * tanh(0.25 * (abs(x)- 30))", xlabel="x", ylabel="y", title="Plot")
"""


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