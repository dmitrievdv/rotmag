using Plots

include("rotmag.jl")
import .Rotmag

star_fast = Rotmag.Star(2.0, 0.5, 35)
star_slow = Rotmag.Star(2.0, 0.5, 1)
mag = Rotmag.Magnetosphere(2, 3, 1e-7, 10)

function testderivative(r_m, θ, star, mag; eps = 1e-3)
    r = r_m*sin(θ)^2
    v_po, v_to = Rotmag.nonsolidrotation(r, θ, star, mag)
    ∇v_p, ∇v_t = Rotmag.velocityjacobian(r, θ, v_po, v_to, star, mag)
    dr = eps; dθ = eps
    v_pr, v_tr = Rotmag.nonsolidrotation(r + dr, θ, star, mag)
    v_pθ, v_tθ = Rotmag.nonsolidrotation(r, θ + dθ, star, mag)
    ∇v_pcalc = [(v_pr - v_po)/dr, (v_pθ - v_po)/dθ]
    ∇v_tcalc = [(v_tr - v_to)/dr, (v_tθ - v_to)/dθ]
    println("∇v_p: ", ∇v_p)
    println("∇v_p: ", ∇v_pcalc)
    println("∇v_t: ", ∇v_t)
    println("∇v_t: ", ∇v_tcalc)
end

# testderivative(2.5, π/4, star, mag, eps = 1e-5)

r_ms, θs, v_ps, v_ts, fast_source_s = Rotmag.solvemag(star_fast, mag, n_r = 3)
Rs = @. r_ms*sin(θs)^3
rs = @. r_ms*sin(θs)^2
plot(rs, fast_source_s, linecolor = :red)
r_ms, θs, v_ps, v_ts, slow_source_s = Rotmag.solvemag(star_slow, mag, n_r = 3)
plot!(rs, slow_source_s, linecolor = :blue)
# plot(Rs, v_ps/1e5, linecolor = :red)
# plot!(Rs, v_ts/1e4, linecolor = :blue)

