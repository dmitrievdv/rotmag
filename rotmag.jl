module Rotmag
using Printf
using Polynomials
using LinearAlgebra

# export solvemag, solvepoint, test_export, R☉

const R☉ = 6.955e10
const M☉ = 1.989e33
const year_seconds = 3.155e7
const G = 6.67259e-8
const σ = 5.77e-5
const c = 2.99792458e10
const h = 6.626176e-27
const kB = 1.380649e-16

struct NoVelocitySolution <: Exception end
struct OuterRadiusTooSmall <: Exception end

Base.showerror(io::IO, e::NoVelocitySolution) = print(io, "can't find velocity solution")
Base.showerror(io::IO, e::NoVelocitySolution) = print(io, "r_mo ≤ r_mi")

struct Star{T<:Real}
	R::T
	M::T
	T::T
	v_esc::T
	v_eq::T
	function Star(R, M, T, v_eq)
		_R = float(R)
		_M = float(M)
		_T = float(T)
		_v_eq = float(v_eq) 
    	new{typeof(_R)}(_R, _M, _T, √(2*G*_M/_R*M☉/R☉), _v_eq)
	end
end

struct Magnetosphere{T<:Real}
	r_mi::T
	r_mo::T
	M_dot::T
	v_start::T
	function Magnetosphere(r_mi, r_mo, M_dot, v_start)
			if r_mo ≤ r_mi
				throw(OuterRadiusTooSmall)
			end
			_r_mi = float(r_mi)
			_r_mo = float(r_mo)
			_M_dot = float(M_dot)
			_v_start = float(v_start)
    	new{typeof(_r_mi)}(_r_mi, _r_mo, _M_dot, _v_start)
	end
end

function dipolermax(r, θ)
	return r/sin(θ)^2
end

function freefallvelocity(r, θ, star)
	return star.v_esc*√(1/r*(1 - sin(θ)^2))
end

function freefallvelocitygradient(r, θ, star)
	dvdθ = -star.v_esc*√(1/r)/√(1-sin(θ)^2)*sin(θ)*cos(θ)
	dvdr = -0.5*star.v_esc*√(1 - sin(θ)^2)*√(1/r^3)
	return [dvdr, dvdθ]
end

function surfmagfield(M_dot, R_in, star)
	return 4.2e2*(R_in/2.2)*(star.M/(0.5))^0.25*(M_dot/1e-8)^0.5*(star.R/(2))^(-3)
end

function dipolefield(surffield, r, θ)
	return surffield*(4-3*sin(θ)^2)^0.5/r^3
end

function  dipolefieldgradient(surffield, r, θ)
	dBdr = -3*surffield*(4-3*sin(θ)^2)^0.5/r^4
	dBdθ = -3*surffield/r^3/(4-3*sin(θ)^2)^0.5*sin(θ)*cos(θ)
	return [dBdr, dBdθ]
end

function eta(M_dot, r_mi, r_mo, star)
	B_star = surfmagfield(M_dot, r_mi, star)
	return 2.6e-9*(B_star/1e3)^(-1)*(M_dot/1e-8)*(star.R/(2))^-2*(1/r_mi-1/r_mo)^(-1)
end

function constructpolynomial(ρ, E, R, v_eq, Bp, η)
	P = Polynomial([ρ^2*E, -2R*ρ*E, (R^2*E - ρ^2*v_eq^2 - (Bp*R/(4*π*η))^2), 
					2*R*ρ*v_eq^2, -R^2*v_eq^2])
    return P
end

function polynomialgradient(ρ, ∇ρ, E, ∇E, R, ∇R, B_p, ∇B_p, v_eq, η)
	∇P = @. Polynomial(vcat(2*ρ*∇ρ*E + ρ^2*∇E, -2(∇R*E*ρ + R*∇E*ρ + R*E*∇ρ), 
						 (2R*∇R*E + R^2*∇E - 2v_eq^2*ρ*∇ρ - 2B_p*R/(4*π*η)^2*(∇B_p*R + B_p*∇R)),
						   2R*∇ρ*v_eq^2 + 2ρ*∇R*v_eq^2, -2R*∇R*v_eq^2))
	return ∇P 
end

function nonsolidrotation(r :: Real, θ :: Real, star :: Star, mag :: Magnetosphere)
	r_m = r/sin(θ)^2
	R = r*sin(θ)
	η = eta(mag.M_dot, mag.r_mi, mag.r_mo, star)
	B_star = surfmagfield(mag.M_dot, mag.r_mi, star)
	Bp = dipolefield(B_star, r, θ)
	v_ff = freefallvelocity(r, θ, star)
	v_eq = star.v_eq*1e5
	v_po = mag.v_start*1e5
	if star.v_eq ≈ 0.0 
		v_p = √(v_po^2 + v_ff^2)
		v_t = 0.0
		return [v_p, v_t]
	end
	ρ = r_m^2 - R^2
	E = v_po^2 + v_ff^2 - v_eq^2*ρ
	P = constructpolynomial(ρ, E, R, v_eq, Bp, η)
	# p = p/(R^2*v_eq^2)
	x = roots(P)
	xreal = real.(x[imag.(x) .≈ 0])
	xreal = xreal*v_eq
	# v_t = sum(abs.(xreal))/2
	k = @. √(E - xreal^2)/Bp
	Bϕ = @. xreal/k
	ang_conservation = @. abs(v_eq*ρ/(Bϕ*R*(k - (4*π*η)^(-1)))-1)
	
	v_t = try xreal[ang_conservation .< 1e-10][1]
	catch e
		if e isa(BoundsError)
			throw(NoVelocitySolution)
		end
	end
	v_p = √(E - v_t^2)
	return [v_p, v_t]
end

function velocityjacobian(r :: Real, θ :: Real, v_p :: Real, v_t :: Real, star :: Star, mag :: Magnetosphere)
	r_m = r/sin(θ)^2
	R = r*sin(θ)
	η = eta(mag.M_dot, mag.r_mi, mag.r_mo, star)
	B_star = surfmagfield(mag.M_dot, mag.r_mi, star)
	B_p = dipolefield(B_star, r, θ)
	∇B_p = dipolefieldgradient(B_star, r, θ)
	v_ff = freefallvelocity(r, θ, star)
	∇v_ff = freefallvelocitygradient(r, θ, star)
	v_eq = star.v_eq*1e5
	v_po = mag.v_start*1e5
	if star.v_eq ≈ 0.0
		∇v_t = [0, 0] 
		∇v_p = [-v_ff^2/v_p/(2r), -v_ff^2/v_p*tan(θ)]
		return ∇v_p, ∇v_t
	end
	∇r_m = [1/sin(θ)^2, -2r/sin(θ)^3*cos(θ)]
	∇R = [sin(θ), r*cos(θ)]
	ρ = r_m^2 - R^2
	∇ρ = @. 2r_m*∇r_m - 2R*∇R
	E = v_po^2 + v_ff^2 - v_eq^2*ρ
	∇E = @. 2*v_ff*∇v_ff - v_eq^2*∇ρ
	P = constructpolynomial(ρ, E, R, v_eq, B_p, η)
	dP = derivative(P)
	∇P = polynomialgradient(ρ, ∇ρ, E, ∇E, R, ∇R, B_p, ∇B_p, v_eq, η)
	x = v_t/v_eq
	∇v_t = -v_eq*[∇P[1](x)/dP(x), ∇P[2](x)/dP(x)]
	∇v_p = @. 1/v_p*(∇E/2 - v_t*∇v_t)
	return ∇v_p, ∇v_t
end

function covariantvelocityjacobian(r, θ, v_p, v_t, ∇v_p, ∇v_t)
	ξ = √(4 - 3sin(θ)^2)
	k = [-2cos(θ)/ξ, -r*sin(θ)/ξ, 0]
	t = [0, 0, r*sin(θ)]
	v = v_p*k + v_t*t

	∇k = [     0         -sin(θ)/ξ   0; 
		  2sin(θ)/ξ^3 -r*4cos(θ)/ξ^3 0;
		       0          0          0]
	
	∇t = [ 0 0  sin(θ);
		   0 0 r*cos(θ);
		   0 0    0]

	k∇v_p = [k[1]*∇v_p[1] k[2]*∇v_p[1] k[3]*∇v_p[1];
			 k[1]*∇v_p[2] k[2]*∇v_p[2] k[3]*∇v_p[2];
				 0            0            0]
				 
	t∇v_t = [0 0 t[3]*∇v_t[1];
			 0 0 t[3]*∇v_t[2];
			 0 0      0]

	∇v = v_p*∇k + k∇v_p + v_t*∇t + t∇v_t
	return v, ∇v
end

function directionalgradient(r, θ, α, β, v, ∇v)
	n = [cos(α), 1/r*sin(α)*cos(β),  1/r/sin(θ)*sin(α)*sin(β)]
	∇n = [      0             -n[2]/r              -n[3]/r;
		      r*n[2]          -n[1]/r           -n[3]/tan(θ);
		  r*sin(θ)^2*n[3] sin(θ)*cos(θ)*n[3] -n[1]/r - n[2]/tan(θ)]

	n∇vn = dot(n,(∇v*n + ∇n*v))
	return n∇vn
end

function integrategradalldirections(r, θ, v_p, v_t, ∇v_p, ∇v_t; n=100)
	αs = [0:π/n:(n-1)*π/n;] .+ π/(2*n)
	βs = [0:2*π/n:(n-1)*2*π/n;] .+ π/n
	I = 0
	v, ∇v = covariantvelocityjacobian(r, θ, v_p, v_t, ∇v_p, ∇v_t)
	for α ∈ αs, β ∈ βs
		I = I + abs(directionalgradient(r, θ, α, β, v, ∇v))*sin(α)
	end
	return I*π/n*2*π/n
end

function integrategradstar(r, θ, v_p, v_t, ∇v_p, ∇v_t; n=100)
	ϕ = 0
	try
		ϕ = asin(1/r)
	catch e
		if isa(e, DomainError)
			ϕ = π/2
		end
	end

	αs = -([0:ϕ/n:(n-1)*ϕ/n;] .- π) .+ ϕ/(2*n)
	βs = [0:2*π/n:(n-1)*2*π/n;] .+ π/n
	I = 0
	v, ∇v = covariantvelocityjacobian(r, θ, v_p, v_t, ∇v_p, ∇v_t)
	for α ∈ αs, β ∈ βs
		I = I + abs(directionalgradient(r, θ, α, β, v, ∇v))*sin(α)
	end
	return I*ϕ/n*2*π/n
end

function calcsource(r, θ, v_p, v_t, ∇v_p, ∇v_t)
	return (integrategradstar(r, θ, v_p, v_t, ∇v_p, ∇v_t)/
	        integrategradalldirections(r, θ, v_p, v_t, ∇v_p, ∇v_t))
end

function calcuplevel(source :: Real, T_star; n_l = 1e5, λ = 1.0830e4, g_u = 5, g_l = 3)
	ν = c/(λ*1e-8)
	# print(h*ν, "\n")
	n_u = g_u/g_l*n_l/((exp(h*ν/(kB*T_star))-1)/source + 1)
	return n_u
end

function solvemag(star :: Star, mag :: Magnetosphere; n_r = 5 :: Int, n_θ = 20 :: Int, modelname = "rotmag")
	mkpath("models/popul")
	mkpath("models/vel")
	mkpath("models/data")
	
	datafile = "models/data/"*modelname*"_data.dat"
    populfile = "models/popul/"*modelname*"_popul.dat"
	velfile = "models/vel/"*modelname*"_vel.dat"
	
    dataout = open(datafile, "w")
    println(dataout, "# Parameters for model "*modelname)
	println(dataout, star.M, " # Mstar")
	println(dataout, star.R, " # Rstar")
	println(dataout, star.T, " # Tstar")
	println(dataout, star.v_eq, " # equatorial rotation")
	println(dataout, mag.M_dot, " # Mdot")
	println(dataout, 7500, " # Tmag")
	println(dataout, "dipole # field type")
	println(dataout, "# borders")
	println(dataout, mag.r_mi, " # first")
	println(dataout, mag.r_mo, " # second")
	close(dataout)

	B_star = surfmagfield(mag.M_dot, mag.r_mi, star)
	θ_grid = zeros((n_θ, n_r))
	r_m_grid = zeros((n_θ, n_r))
	v_p_grid = zeros((n_θ, n_r))
	v_t_grid = zeros((n_θ, n_r))
	source_grid = zeros((n_θ, n_r))
	n_u_grid = zeros((n_θ, n_r))
	n_l_grid = fill(1e5, (n_θ, n_r))

	for i=1:n_r
		r_m = mag.r_mi + (i-1)*(mag.r_mo-mag.r_mi)/(n_r-1)
		θ_s = asin(√(1/r_m))
		# solve_start = [star.v_esc/B_star, 0]
		for j=1:n_θ
			θ = π/2 - (j)*(π/2 - θ_s)/(n_θ)
			r = r_m*sin(θ)^2
			r_m_grid[j,i] = r_m
			θ_grid[j,i] = θ
			B_p = dipolefield(B_star, r_m*sin(θ)^2, r_m)
			v_p, v_t = nonsolidrotation(r, θ, star, mag)
			v_t_grid[j,i] = v_t
			v_p_grid[j,i] = v_p
			∇v_p, ∇v_t = velocityjacobian(r, θ, v_p, v_t, star, mag)
			source_grid[j,i] = calcsource(r, θ, v_p, v_t, ∇v_p, ∇v_t)
			n_u_grid[j,i] = calcuplevel(source_grid[j,i], star.T) 
		end
	end

	velout = open(velfile, "w")
	println(velout, "# ", n_r, " ", n_θ)
	println(velout, "#")
	populout = open(populfile, "w")
	println(populout, "# ", n_r, " ", n_θ)
	println(populout, "#")
	for i=1:n_r
		r_m = mag.r_mi + (i-1)*(mag.r_mo-mag.r_mi)/(n_r-1)
		θ_s = asin(√(1/r_m))
		# solve_start = [star.v_esc/B_star, 0]
		for j=n_θ:-1:1
			@printf(populout, "%.2f %.7f %.4e %.4e %.4e %.4e \n", 
			    r_m_grid[j,i], θ_grid[j,i], 1e4, source_grid[j,i], n_l_grid[j,i], n_u_grid[j,i])
			@printf(velout, "%.2f %.7f %.4e %.4e\n", 
			    r_m_grid[j,i], θ_grid[j,i], v_p_grid[j,i], v_t_grid[j,i])
		end
		print(populout, "\n")
		print(velout, "\n")
	end
	close(velout)
	close(populout)
	return r_m_grid, θ_grid, v_p_grid, v_t_grid, n_u_grid, n_l_grid
end

function magvel(star :: Star, mag :: Magnetosphere; n_r = 5 :: Int, n_θ = 20 :: Int, file = "rotmag.dat")
	# out = open(file, "w")
	B_star = surfmagfield(mag.M_dot, mag.r_mi, star)
	θ_grid = zeros((n_θ, n_r))
	r_m_grid = zeros((n_θ, n_r))
	v_p_grid = zeros((n_θ, n_r))
	v_t_grid = zeros((n_θ, n_r))
	source_grid = zeros((n_θ, n_r))
	for i=1:n_r
		r_m = mag.r_mi + (i-1)*(mag.r_mo-mag.r_mi)/(n_r-1)
		θ_s = asin(√(1/r_m))
		# solve_start = [star.v_esc/B_star, 0]
		for j=1:n_θ
			θ = π/2 - (j)*(π/2 - θ_s)/(n_θ)
			r = r_m*sin(θ)^2
			r_m_grid[j,i] = r_m
			θ_grid[j,i] = θ
			B_p = dipolefield(B_star, r_m*sin(θ)^2, r_m)
			v_p, v_t = nonsolidrotation(r, θ, star, mag)
			v_t_grid[j,i] = v_t
			v_p_grid[j,i] = v_p
			v_rot = star.v_eq*r_m*sin(θ)^3
			# @printf(out, "%8.f %8.f %8.f %8.f %8.f %8.f\n", r_m, θ, v_p*1e-5, v_t*1e-5, v_ff*1e-5, v_rot*1e-5)
		end
		# print(out, "\n")
	end
	# close(out)
	return r_m_grid, θ_grid, v_p_grid, v_t_grid
end
end


