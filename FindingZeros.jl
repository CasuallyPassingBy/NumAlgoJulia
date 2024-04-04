
"""
    bisection(f, a, b; tol=1e-6, max_iter=1000, bracketing=false)

Find a zero of the function `f` using the bisection method within the interval `[a, b]`.

# Arguments
- `f::Function`: The function whose zero is to be found.
- `a::Number`: The left endpoint of the interval.
- `b::Number`: The right endpoint of the interval.
- `tol::Real`: The tolerance for the zero. The iteration stops when `|f(c)| ≤ tol`, where `c` is the midpoint of the interval.
- `max_iter::Integer`: The maximum number of iterations allowed.
- `bracketing::Bool`: If `true`, return the bracketing interval instead of the approximate zero.

# Returns
- If `bracketing` is `false`, returns the approximate zero of the function within the specified tolerance.
- If `bracketing` is `true`, returns the bracketing interval `[a, b]`.

# Examples
```julia
f(x) = x^2 - 2
bisection(f, 1, 2)  # Output: 1.4142131805419922
bisection(f, 1, 2, bracketing=true)  # Output: [1.4142131805419922, 1.4142141342163086]
```
"""
function bisection(f, a, b; tol=1e-6, max_iter=1000, bracketing = false)
    # Check if f(a) and f(b) have opposite signs
    if sign(f(a)) * sign(f(b)) > 0
        error()
    end
    
    # Initialize iteration counter and midpoint
    iter = 0
    c = (a + b) / 2
    
    while abs(f(c)) > tol && iter < max_iter
        c = a + (b-a)/2
        fp = f(c)
        # Check if the root is at midpoint
        if fp == 0 || (b-a)/2 < tol
            if !bracketing
                return c   
            else
                if sign(fp) * sign(f(a)) < 0 # We use sign to avoid underflow/overflow
                    return sort!([a, c])
                else
                    return sort!([c, b])
                end
            end
        end
        # Update the interval
        if sign(fp) * sign(f(a)) < 0 # We use sign to avoid underflow/overflow
            b = c
        else
            a = c
        end
        iter += 1
    end
    error("Did not converge within the maximum iterations.")
end
"""
    newton_method(f, f_prime, x0, tol=1e-6, max_iter=1000)

Find a zero of the function `f` using Newton's method, starting from the initial approximation `x0`.

# Arguments
- `f::Function`: The function whose zero is to be found.
- `f_prime::Function`: The derivative of the function `f`.
- `x0::Number`: The initial approximation.
- `tol::Real`: The tolerance for convergence. The iteration stops when `|x_next - x_current| ≤ tol`.
- `max_iter::Integer`: The maximum number of iterations allowed.

# Returns
- `x_next::Number`: The approximate zero of the function within the specified tolerance.

# Notes
If the absolute value of the derivative at the current approximation is less than `tol`, the method terminates with a warning.

# Examples
```julia
f(x) = x^2 - 2
f_prime(x) = 2x
newton_method(f, f_prime, 1.5) # Output 1.4142135623730951
```
"""
function newton_method(f, f_prime, x0; tol=1e-6, max_iter=1000)
    # Initialize the iteration counter and the initial approximation
    x_current = x0
    
    for iter in 1:max_iter
        f_val = f(x_current)       # Evaluate the function
        f_prime_val = f_prime(x_current)  # Evaluate the derivative
        
				# Check if the derivative is close to zero
        if abs(f_prime_val) < tol
            println("Derivative close to zero. Newton's method cannot proceed.")
            return x_current
        end
				# Update the approximation using Newton's method
        x_next = x_current - f_val / f_prime_val

        # Check for convergence
        if abs(x_next - x_current) < tol # we can change this inequality for any of the other 3, 
            return x_next
        end
				
        x_current = x_next
    end
    
    error("Did not converge within the maximum iterations.")
    return x_current
end
"""
    function secant_method(f, x0, x1; tol = 1e-6, max_iter = 1000)

Find a zero of the function `f` using the Secant method, starting from the initial approximations `x0` and `x1`.

# Arguments
- `f::Function`: the fucntion whose zero is to be found.
- `x0::Number`: An initial approximation
- `x1::Number`: Another initial approximation
- `tol::Real`: The tolerance for convergence. The iteration stops when `|x_next - x_current| ≤ tol`.
- `max_iter::Integer`: The maximum number of iterations allowed.

# Returns
- `x_current::Number` The appriximate zero fo the function within the specified tolerance.

# Notes 
If the absolute value of the difference x_current and x_prev is less than the `tol`, the method terminates with a warning.

# Examples
```julia
f(x) = x^2 -2
secant_method(f, 1, 2) # Output 1.4142135623730954
```
"""
function secant_method(f, x0, x1; tol = 1e-6, max_iter = 1000)
    # Initialize the iteration counter and the initial approximations
    x_prev = x0
    x_current = x1
    
    for iter in  1:max_iter
        f_prev = f(x_prev)       # Evaluate the function at the previous approximation
        f_current = f(x_current) # Evaluate the function at the current approximation
        
        # Check for convergence
        if abs(x_prev - x_current) < tol  # this can change to any of the 3 inqualites to check for convergence
            return x_current
        end
        
        # Check if the difference is close to zero
        if abs(x_current - x_prev) < tol
            println("Difference between approximations close to zero. Secant method cannot proceed.")
            return x_current
        end
        
        # Update the approximation using the secant method
        x_next = x_current - f_current * (x_current - x_prev) / (f_current - f_prev)
        x_prev = x_current
        x_current = x_next
    end
    
    println("Did not converge within the maximum iterations.")
    return x_current
end

"""
    false_position(f, a, b; tol=1e-6, max_iter=1000, bracketing=false)

Find a zero of the function `f` using the false position method within the interval `[a, b]`.

# Arguments
- `f::Function`: The function whose zero is to be found.
- `a::Number`: The left endpoint of the interval.
- `b::Number`: The right endpoint of the interval.
- `tol::Real`: The tolerance for convergence. The iteration stops when `|f(x_next)| ≤ tol`.
- `max_iter::Integer`: The maximum number of iterations allowed.
- `bracketing::Bool`: If `true`, return the bracketing interval instead of the approximate zero.

# Returns
- If `bracketing` is `false`, returns the approximate zero of the function within the specified tolerance.
- If `bracketing` is `true`, returns the bracketing interval `[a, b]`.

# Examples
```julia
f(x) = x^2 - 2
false_position(f, 1.0, 2.0)  # Output: 1.4142131805419922
false_position(f, 1.0, 2.0; bracketing=true)  # Output: [1.41421143847487, 1.4142135620573204]
```
"""
function false_position(f, a, b; tol=1e-6, max_iter=1000, bracketing = false)
	# Check if f(a) and f(b) have opposite signs
	if f(a) * f(b) > 0
		error("Function has the same sign at both endpoints. False position method cannot proceed.")
	end
	
	# Initialize iteration counter and variables
	iter = 0
	x_prev = a
	x_current = b
	
	while iter < max_iter
		f_prev = f(x_prev)
		f_current = f(x_current)
		
		# Check for convergence
		if abs(f_current) < tol
			if bracketing
                bracket = [x_prev, x_current]
                return sort!(bracket)
            else
                return x_current
            end
		end
		
		# Check if the difference is close to zero
		if abs(x_current - x_prev) < tol
			println("Difference between approximations close to zero. False position method cannot proceed.")
			return x_current
		end
		
		# Update the approximation using the false position method
		x_next = x_current - f_current * (x_current - x_prev) / (f_current - f_prev)
		
		# Determine the new interval
		if f(x_next) * f(a) < 0
			b = x_next
		else
			a = x_next
		end
		
		# Update variables for the next iteration
		x_prev = x_current
		x_current = x_next
		
		iter += 1
	end
	
	println("Did not converge within the maximum iterations.")
	return x_current
end
"""
    modified_newton(f, f_prime, f_double_prime, x0, tol=1e-6, max_iter=1000)

Find a zero of the function `f` using the modified Newton's method, which converges quadratically regardless of the multiplicity of the zero.

# Arguments
- `f::Function`: The function whose zero is to be found.
- `f_prime::Function`: The first derivative of the function `f`.
- `f_double_prime::Function`: The second derivative of the function `f`.
- `x0::Number`: The initial approximation.
- `tol::Real`: The tolerance for convergence. The iteration stops when `|f(x_next)| ≤ tol`.
- `max_iter::Integer`: The maximum number of iterations allowed.

# Returns
- `x_next::Number`: The approximate zero of the function within the specified tolerance.

# Notes
If the denominator in the iterative formula becomes close to zero, the method terminates with a warning.

# Examples
``` julia
f(x) = x^2 - 2
f_prime(x) = 2x
f_double_prime(x) = 2
modified_newton(f, f_prime, f_double_prime, 1.5)  # Output: 1.4142135623746899
```
"""
function modified_newton(f, f_prime, f_double_prime, x0, tol=1e-6, max_iter=1000)
    # Initialize the iteration counter and the initial approximation
    iter = 0
    x_current = x0
    
    while iter < max_iter
        f_val = f(x_current)       # Evaluate the function
        f_prime_val = f_prime(x_current)  # Evaluate the first derivative
        f_double_prime_val = f_double_prime(x_current)  # Evaluate the second derivative
        
        # Check for convergence
        if abs(x_current - x_next) < tol # This can be changed for different inequalities
            return x_current
        end
        
        # Check if the derivative is close to zero
        if abs(f_double_prime_val * f_prime_val + f_val) < tol
            println("The denominator is close to zero. Modified Newton's method cannot proceed.")
            return x_current
        end
        
        # Update the approximation using the modified Newton's method
        x_next = x_current - f_val * (f_prime_val / (f_double_prime_val * f_prime_val + f_val))
        
        x_current = x_next  # Update the approximation
        
        iter += 1
    end
    
    println("Did not converge within the maximum iterations.")
    return x_current
end
"""
    steffensen_acceleration(f, x0; tol=1e-6, max_iter=1000)

Find a fixed point of the function `f` using Steffensen's method, which accelerates convergence by applying quadratic interpolation.

# Arguments
- `f::Function`: The function whose zero is to be found.
- `x0::Number`: The initial approximation.
- `tol::Real`: The tolerance for convergence. The iteration stops when `|x_accelerated - x_current| ≤ tol`.
- `max_iter::Integer`: The maximum number of iterations allowed.

# Returns
- `x_accelerated::Number`: The approximate zero of the function within the specified tolerance.

# Examples
```julia
f(x) = x^2 - 2
steffensen_acceleration(f, 1)  # Output: -1
```
"""
function steffensen_acceleration(f, x0; tol=1e-6, max_iter=1000)
    x_current = x0
    
    for iter in 1:max_iter
        x_next = f(x_current)
        x_next_next = f(x_next)
        
        if abs(x_next_next - 2x_next + x_current) < tol
            return x_current
        end
        # Apply Steffensen's acceleration
        x_accelerated = x_current - ((x_next - x_current)^2) / (x_next_next - 2x_next + x_current)
        
        # Check for convergence
        if abs(x_accelerated - x_current) < tol
            return x_accelerated
        end
        
        x_current = x_accelerated  # Update the approximation
    end
    
    println("Did not converge within the maximum iterations.")
    return x_current
end
"""
    accelerated_newton(f, f_prime, x0, tol=1e-6, max_iter=1000)

Find a zero of the function `f` using Newton's method with Steffensen's acceleration.

# Arguments
- `f::Function`: The function whose zero is to be found.
- `f_prime::Function`: The first derivative of the function `f`.
- `x0::Number`: The initial approximation.
- `tol::Real`: The tolerance for convergence. The iteration stops when `|x_accelerated - x_current| ≤ tol`.
- `max_iter::Integer`: The maximum number of iterations allowed.

# Returns
- `x_accelerated::Number`: The approximate zero of the function within the specified tolerance.

# Notes
- If the derivative of the function is too close to zero at any point during the iteration, Newton's method cannot proceed, and the function returns the current approximation.
- If convergence is not achieved within the maximum iterations, the function throws an error.

# Examples
```julia
f(x) = x^2 - 2
f_prime(x) = 2x
accelerated_newton(f, f_prime, 1.5)  # Output: 1.4142135623746899
```
"""
function accelerated_newton(f, f_prime, x0, tol=1e-6, max_iter=1000)
    # Initialize the current approximation
    x_current = x0
    
    # Define a nested function for Newton's method
    function newton_function(x0)
        # Evaluate the function and its derivative at the current approximation
        f_val = f(x0)
        f_prime_val = f_prime(x0)
        
        # Check if the derivative is too close to zero, which may cause division by zero
        if abs(f_prime_val) < tol
            println("Derivative close to zero. Newton's method cannot proceed.")
            error()
        end
        
        # Compute the next approximation using Newton's method
        return x0 - f_val / f_prime_val
    end
    
    # Initialize variables for storing next two iterations
    x_next = x_next_next = 0
    
    # Iterate until convergence or reaching the maximum number of iterations
    for iter in 1:max_iter
        try
            # Apply Newton's method to compute the next two iterations
            x_next = newton_function(x_current)
            x_next_next = newton_function(x_next)
        catch
            # Handle potential errors (e.g., division by zero) by returning the current approximation
            return x_current
        end
        
        # Apply Steffensen's acceleration to the current approximation
        if abs(x_next_next - 2x_next + x_current) < tol
            return x_current
        end
        x_accelerated = x_current - ((x_next - x_current)^2) / (x_next_next - 2x_next + x_current)
        # Check for convergence based on the difference between consecutive approximations
        if abs(x_accelerated - x_current) < tol
            return x_accelerated
        end
        
        # Update the current approximation
        x_current = x_accelerated
    end
    
    # If convergence is not achieved within the maximum iterations, print a message and return the current approximation
    error("Did not converge within the maximum iterations.")
end

"""
    fixed_point_iteration(f, x0, tol=1e-6, max_iter=1000)

Find a fixed point of the function `f` using the fixed-point iteration method.

# Arguments
- `f::Function`: The function for which the fixed point is to be found.
- `x0::Number`: The initial approximation.
- `tol::Real`: The tolerance for convergence. The iteration stops when `|x_next - x_current| ≤ tol`.
- `max_iter::Integer`: The maximum number of iterations allowed.

# Returns
- `x_fixed::Number`: The approximate fixed point of the function within the specified tolerance.

# Notes
- If convergence is not achieved within the maximum iterations, the function throws an error.

# Examples
```julia
f(x) = x^2 - 2
fixed_point_iteration(f, 1.5; max_iter = 10000)  # Output: 1.9999992794691956
```
"""
function fixed_point_iteration(f, x0; tol=1e-6, max_iter=1000)
    # Initialize the current approximation
    x_current = x0
    # Iterate until convergence or reaching the maximum number of iterations
    for iter in 1:max_iter
        # Compute the next approximation using fixed-point iteration
        x_next = f(x_current)
        
        # Check for convergence
        if abs(x_next - x_current) < tol
            return x_next
        end
        
        # Update the current approximation
        x_current = x_next
    end

    # If convergence is not achieved within the maximum iterations, throw an error
    error("Did not converge within the maximum iterations.")
end

"""
    mullers_method(f, x0, x1, x2; tol=1e-10, max_iter=1000)

Find a root of the function `f` using Müller's method.

# Arguments
- `f::Function`: The function for which to find the root.
- `x0::Number`: The first initial guess.
- `x1::Number`: The second initial guess.
- `x2::Number`: The third initial guess.
- `tol::Real`: The tolerance for convergence. The iteration stops when `|x - x2| ≤ tol`.
- `max_iter::Integer`: The maximum number of iterations allowed.

# Returns
- `x::Number`: The approximate root of the function within the specified tolerance.
- `iter::Integer`: The number of iterations performed.

# Notes
- Müller's method is an iterative numerical root-finding algorithm that utilizes quadratic interpolation.
- The initial guesses `x0`, `x1`, and `x2` should be chosen such that `x1` is closer to the root than `x0` and `x2`.
- If convergence is not achieved within the maximum iterations, the function throws an error.

# Examples
```julia
f(x) = x^2 - 2
x0, x1, x2 = 1.0, 1.5, 2.0
root, iterations = mullers_method(f, x0, x1, x2)
```
"""
function mullers_method(f, x0, x1, x2; tol=1e-10, max_iter=1000)
    # Initialize iteration counter
    iter = 0
    
   for iter in 1:max_iter
        # Compute differences and slopes
        h1 = x1 - x0
        h2 = x2 - x1
        delta1 = (f(x1) - f(x0)) / h1
        delta2 = (f(x2) - f(x1)) / h2
        
        # Compute coefficients of the quadratic interpolation
        a = (delta2 - delta1) / (h2 + h1)
        b = a * h2 + delta2
        c = f(x2)
        
        # Calculate the discriminant
        radicand = b^2 - 4*a*c
        
        # Choose the root with smaller denominator
        if abs(b + sqrt(radicand)) > abs(b - sqrt(radicand))
            den = b + sqrt(radicand)
        else
            den = b - sqrt(radicand)
        end
        
        # Update root estimate
        x = x2 - 2c / den
        
        # Check for convergence
        if abs(x - x2) < tol
            return x
        end
        
        # Update points for the next iteration
        x0, x1, x2 = x1, x2, x
    end
    
    error("Müller's method did not converge within the specified number of iterations.")
end

"""
    iterated_inverse_interpolation(x_values, y_values; target_y = 0)

Compute the interpolated x value corresponding to a given y value using the iterated inverse interpolation method.

# Arguments
- `x_values::Vector`: A vector of x values.
- `y_values::Vector`: A vector of corresponding y values.
- `target_y::Number`: The target y value for which the corresponding x value is to be computed. Default is 0.

# Returns
- `x_interpolated::Number`: The interpolated x value corresponding to the target y value.

# Examples
```julia
x_values = [0, 1, 2, 3, 4]
y_values = [0, 1, 4, 9, 16]
target_y = 8
iterated_inverse_interpolation(x_values, y_values; target_y) # Output
```
"""
function iterated_inverse_interpolation(x_values, y_values; target_y = 0)
    n = length(x_values)
    Q = similar(x_values, Float64) # Initialize the first column of the table with y-values
    for j in 2:n
        for i in 1:n-j+1
            Q[i] = ((target_y - y_values[i]) * Q[i+1] - (target_y - y_values[i+j-1]) * Q[i]) / (y_values[i+j-1] - y_values[i])
        end
    end
    return Q[1] # The interpolated value is in the top-right corner of the table
end