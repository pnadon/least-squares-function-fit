using Pkg; Pkg.add("QuadGK")
using QuadGK
using LinearAlgebra

function main()
	func( t) = 2*t^2 + 4
	bounds = [ 0, 2]
	num_digits = 3

	x = approx_func( func, bounds, num_digits)
	display( x)
end

function approx_func( func, bounds, num_digits)
	lower, upper = bounds

	coefs = [
		[ ( upper - lower) ((upper^2 - lower^2) / 2) (( upper^3 - lower^3) / 3)];
		[ ((upper^2 - lower^2) / 2) (( upper^3 - lower^3) / 3) (( upper^4 - lower^4) / 4)];
		[ ((upper^3 - lower^3) / 3) (( upper^4 - lower^4) / 4) (( upper^5 - lower^5) / 5)]]

	right_side = [ 
		quadgk( func, lower, upper, rtol=1e-5);
		quadgk( t -> t * func(t), lower, upper, rtol=1e-5);
		quadgk( t -> t^2 * func(t), lower, upper, rtol=1e-5)]
	right_side = first.(right_side)
	
	L,U = lu( coefs, Val(false))
	
	y = solve_lower( L, right_side)
	x = solve_upper( U, y)

	x
end

function solve_lower( A, b)
	x1 = b[1] / A[1,1]
	x2 = (b[2] - A[2,1]*x1) / A[2,2]
	x3 = (b[3] - A[3,1]*x1 - A[3,2]*x2) / A[3,3]

	[ x1; x2; x3]
end

function solve_upper( A, b)
	x3 = b[3] / A[3,3]
	x2 = (b[2] - A[2,3]*x3) / A[2,2]
	x1 = (b[1] - A[1,3]*x3 - A[1,2]*x2) / A[1,1]

	x = [ x1; x2; x3]
	x = map( x -> round(x; sigdigits = 3), x)
end

main()