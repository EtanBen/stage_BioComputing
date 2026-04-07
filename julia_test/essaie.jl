using Plots

x = 0:0.1:10
y1 = sin.(x)
y2 = exp.(0.2*x)

p1 = plot(x, y1, color=:blue, label="sin(x)")

p2 = twinx()
plot!(p2, x, y2, color=:red, label="exp(0.2x)")

display(p1)