using Plots

include("rotmag.jl")
import .Rotmag

star = Rotmag.Star(2.0, 0.5, 15)
mag = Rotmag.Magnetosphere(2.2, 3, 1e-7, 10)
r_ms, θs, v_ps, v_ts = Rotmag.solvemag(star, mag)
Rs = @. r_ms*sin(θs)^3


plot(Rs, v_ps/1e5, linecolor = :red)
plot!(Rs, v_ts/1e4, linecolor = :blue)


