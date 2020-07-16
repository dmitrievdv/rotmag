include("rotmag.jl")
import .Rotmag

star = Rotmag.Star(2.0, 0.5, 15)
mag = Rotmag.Magnetosphere(2.2, 3, 1e-7, 10)
Rotmag.solvemag(star, mag)