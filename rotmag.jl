module Rotmag
using Printf
using Polynomials

# export solvemag, solvepoint, test_export, R☉


const R☉ = 6.955e10
const M☉ = 1.989e33
const year_seconds = 3.155e7
const G = 6.67259e-8
const σ = 5.77e-5
const test_export = 0.223456789

struct Star{T<:Real}
	R::T
	M::T
	v_esc::T
	v_eq::T
	function Star(R, M, v_eq)
		_R = float(R)
		_M = float(M)
		_v_eq = float(v_eq) 
    	new{typeof(_R)}(_R*R☉, _M*M☉, √(2*G*_M/_R*M☉/R☉), _v_eq*1e5)
	end
end

struct Magnetosphere{T<:Real}
	r_mi::T
	r_mo::T
	M_dot::T
	v_start::T
	function Magnetosphere(r_mi, r_mo, M_dot, v_start)
			_r_mi = float(r_mi)
			_r_mo = float(r_mo)
			_M_dot = float(M_dot)
			_v_start = float(v_start)
    	new{typeof(r_mi)}(_r_mi, _r_mo, _M_dot, _v_start*1e5)
	end
end


function dipolermax(r, θ)
	return r/sin(θ)^2
end

function freefallvelocity(r, r_max, star)
	return star.v_esc*√(1/r - 1/r_max)
end

function surfmagfield(M_dot, R_in, star)
	return 4.2e2*(R_in/2.2)*(star.M/(0.5*M☉))^0.25*(M_dot/1e-8)^0.5*(star.R/(2*R☉))^(-3)
end

function dipolefield(surffield, r, r_max)
	return surffield*(4-3*r/r_max)^0.5/r^3
end

function eta(M_dot, r_mi, r_mo, star)
	B_star = surfmagfield(M_dot, r_mi, star)
	return 2.6e-9*(B_star/1e3)^(-1)*(M_dot/1e-8)*(star.R/(2*R☉))^-2*(1/r_mi-1/r_mo)^(-1)
end

function solvepoint(r_m :: Real, θ :: Real, star :: Star, mag :: Magnetosphere)
	r = r_m*sin(θ)^2
	R = r*sin(θ)
	η = eta(mag.M_dot, mag.r_mi, mag.r_mo, star)
	B_star = surfmagfield(mag.M_dot, mag.r_mi, star)
	# solve_start = [star.v_esc/B_star, 0]
	Bp = dipolefield(B_star, r, r_m)
	v_ff = freefallvelocity(r, r_m, star)
	v_eq = star.v_eq
	v_po = mag.v_start
	r_m2R2 = r_m^2 - R^2
	mech_E = v_po^2 + v_ff^2 - v_eq^2*r_m2R2
	# print(B_star, " ", Bp, " ", η, " ", 1e-5*v_ff, "\n")
	p = Polynomial([(r_m2R2)^2*mech_E, -2R*(r_m2R2)*mech_E, (R^2*mech_E - r_m2R2^2*v_eq^2 - (Bp*R/(4*π*η))^2), 
		            2*R*(r_m2R2)*v_eq^2, -R^2*v_eq^2])
	p = p/(R^2*v_eq^2)
	# println(Bp/η)
	x = roots(p)
	# println(x)
	xreal = real.(x[imag.(x) .≈ 0])
	# println(p.(xreal))
	# println(xreal)
	xreal = xreal*v_eq
	vt = sum(abs.(xreal))/2
	k = @. √(mech_E - xreal^2)/Bp
	Bϕ = @. xreal/k
	# println("k: ", k)
	# println("Bϕ: ", Bϕ)
	# println(mech_E, " ", xreal.^2, " ", Bp^2)
	ang_conservation = @. abs(v_eq*r_m2R2/(Bϕ*R*(k - (4*π*η)^(-1)))-1)
	# eng_conservation = @. (k^2*Bp*2 + xreal^2)-(mech_E)
	# println("ang: ", ang_conservation)
	# println("eng: ", eng_conservation)
	vt = xreal[ang_conservation .< 1e-10][1]
    # println(R, " ", r_m, " ", vt/1e5, " ", k*Bp/1e5)
	return vt
end

function solvemag(star :: Star, mag :: Magnetosphere; n_r = 5 :: Int, n_θ = 20 :: Int, file = "rotmag.dat")
	out = open(file, "w")
	B_star = surfmagfield(mag.M_dot, mag.r_mi, star)		
	for i=1:n_r
		r_m = mag.r_mi + (i-1)*(mag.r_mo-mag.r_mi)/(n_r-1)
		θ_s = asin(√(1/r_m))
		# solve_start = [star.v_esc/B_star, 0]
		for j=1:n_θ
			θ = π/2 - (j)*(π/2 - θ_s)/(n_θ)
			B_p = dipolefield(B_star, r_m*sin(θ)^2, r_m)
			v_t = solvepoint(r_m, θ, star, mag)
			v_ff = freefallvelocity(r_m*sin(θ)^2, r_m, star)
			v_p = √(mag.v_start^2 + v_ff^2 - v_t^2)
			v_rot = star.v_eq*r_m*sin(θ)^3
			solve_start = [v_p/B_p, v_t/v_p*B_p]
			@printf(out, "%8.f %8.f %8.f %8.f %8.f %8.f\n", r_m, θ, v_p*1e-5, v_t*1e-5, v_ff*1e-5, v_rot*1e-5)
		end
		print(out, "\n")
	end
	close(out)
	return 0
end
end


