

module itp_insuline
using CSV, DataFrames, Plots, Distributions, Statistics, Catalyst,LsqFit,BSplineKit

export Insuline
data_insuline = CSV.read("P_RT_I_preResecA.csv",DataFrame)
dataInsulineTime    = data_insuline[2:end,1]
dataInsulineMean = mean(Array(data_insuline[2:end,2:end]),dims = 2)

Time_read = [0, 2, 4, 9, 14, 29,180 ]
data_read = [39.97545454545455; 44.44909090909092; 44.88045454545455; 51.48863636363637; 58.06318181818183; 36.06818181818182; 7.458636363636364]
#insulineTimes = data_insuline[2:end,1]

phi(x) = -1490.6 ./ ((-932.3 ./ (x .+ 3.8836)) .- x)
phi_2(x)= (x * exp(x * -0.029063)) + 8.6344
phi_3(x₁)   = log(exp((x₁ * -0.038806) + 18.215) + x₁)
xdata = Float64.([0, 15, 30, 60, 90, 120, 180, 240, 300, 360])


y_G_preRessec = interpolate(vec(Float64.(dataInsulineTime)), vec(Float64.(dataInsulineMean)), BSplineOrder(4))
x = 0:1:360
"""
using SymbolicRegression
x =Float64.(xdata)
X = reshape(x, 1, :)

options= Options(
    binary_operators=[+, -, *,^],
    unary_operators=[exp],
    maxsize=20
)

hall_of_fame = equation_search(
    X,
    ydata,
    niterations=1000,
    options=options
) 

calculate_pareto_frontier(hall_of_fame)


plot(xdata, ydata, seriestype=:scatter, label="Data")
plot!(x, phi.(x), label="phi")
plot!(x, phi_2.(x), label="phi 2")
#plot!(x, phi_3.(x), label="phi 3")

"""

Insuline(x) = y_G_preRessec(x)

x_fine = 0:1:360
y_fine = y_G_preRessec.(x_fine)

plot(dataInsulineTime, dataInsulineMean, seriestype=:scatter)
plot!(x_fine, y_fine, lw=2)

#plot(dataInsulineTime,dataInsulineMean,l=:scatter,label="insuline RT")
plot(dataInsulineTime, dataInsulineMean, seriestype=:scatter)
plot!(x_fine, y_fine, lw=2)
end



