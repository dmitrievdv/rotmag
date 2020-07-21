# Specific file for testing purposes. Test Rotmag however you want!

using Plots

include("rotmag.jl")
import .Rotmag

star_fast = Rotmag.Star(v_eq = 35)
star_slow = Rotmag.Star(v_eq = 1)
mag = Rotmag.Magnetosphere(2, 3, 1e-7, 10)

function testderivative(r_m, star, mag; eps = 1e-6, n_θ = 20)
    θ_s = asin(1/r_m)
    dr = eps; dθ = eps
    Δp = zeros(2); Δt = zeros(2)
    for j=1:n_θ
        θ = π/2 - (j)*(π/2 - θ_s)/(n_θ)
        r = r_m*sin(θ)^2
        v_po, v_to = Rotmag.nonsolidrotation(r, θ, star, mag)
        ∇v_p, ∇v_t = Rotmag.velocityjacobian(r, θ, v_po, v_to, star, mag)
        v_pr, v_tr = Rotmag.nonsolidrotation(r + dr, θ, star, mag)
        v_pθ, v_tθ = Rotmag.nonsolidrotation(r, θ + dθ, star, mag)
        ∇v_pcalc = [(v_pr - v_po)/dr, (v_pθ - v_po)/dθ]
        ∇v_tcalc = [(v_tr - v_to)/dr, (v_tθ - v_to)/dθ]
        Δp = @. Δp + (∇v_p - ∇v_pcalc)^2/∇v_p^2
        Δt = @. Δt + (∇v_t - ∇v_tcalc)^2/∇v_t^2
    end
    Δp = @. √(Δp)
    Δt = @. √(Δt)
    return [Δp Δt]
end


# testderivative(2.5, star_fast, mag, eps = 1e-5)

r_ms, θs, v_ps, v_ts, n_u, n_l = Rotmag.solvemag(star_fast, mag, n_r = 3)
Rs = @. r_ms*sin(θs)^3
rs = @. r_ms*sin(θs)^2
plot(rs, n_u, linecolor = :red)

# r_ms, θs, v_ps, v_ts, slow_source_s = Rotmag.solvemag(star_slow, mag, n_r = 3)
# plot!(rs, slow_source_s, linecolor = :blue)
# r_ms, θs, fv_ps, fv_ts = Rotmag.magvel(star_fast, mag, n_r = 3)
# r_ms, θs, sv_ps, sv_ts = Rotmag.magvel(star_slow, mag, n_r = 3)
# Rs = @. r_ms*sin(θs)^3
# plot(Rs, fv_ts/1e5, linecolor = :red)
# plot!(Rs, sv_ts/1e5, linecolor = :blue)


