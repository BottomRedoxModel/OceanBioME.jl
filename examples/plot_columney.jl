using CairoMakie
using Oceananigans
using Oceananigans.Units

filename = "columney_out"

# ## Load saved output
PHY = FieldTimeSeries("$filename.jld2", "PHY")
NUT = FieldTimeSeries("$filename.jld2", "NUT")
HET = FieldTimeSeries("$filename.jld2", "HET")

# x, y, z = nodes(PHY)
times = PHY.times

times[end] / days

fig = Figure(size = (1000, 1500), fontsize = 20)
z = 0:-4:-200
axis_kwargs = (xlabel = "Time (days)", ylabel = "z (m)", limits = ((0, times[end] / days), (-100meters, 0)))

axPHY = Axis(fig[1, 1]; title = "Phytoplankton concentration (mmol N / m³)", axis_kwargs...)
hmPHY = heatmap!(times / days, z, interior(PHY, 1, 1, :, :)', colormap = :batlow)
Colorbar(fig[1, 2], hmPHY)

# axNUT = Axis(fig[1, 1]; title = "NUT concentration (mmol N / m³)")
# hmNUT = heatmap!(times / days, z, interior(NUT, 1, 1, :, :)', colormap = :batlow)
# Colorbar(fig[1, 2], hmNUT)

fig