module Interpolating
    """
        lagrange_polynomial(x_values::Vector{Float64}, y_values::Vector{Float64})::Function

    Compute the Lagrange interpolation polynomial based on the provided data points.

    # Arguments
    - `x_values::Vector`: Vector of x values representing the data points.
    - `y_values::Vector`: Vector of corresponding y values.

    # Returns
    - `lagrange_interpolator::Function`: A function that computes the Lagrange interpolation polynomial for a given x value.

    # Example
    ```julia
    x_values = [0, 1, 2, 3, 4]
    y_values = [0, 1, 4, 9, 16]
    lagrange_polynomial = lagrange_polynomial(x_values, y_values)
    lagrange_polynomial(2.5) # Output: 6.25
    ```
    """
    function lagrange_polynomial(x_values::Vector{Float64}, y_values::Vector{Float64})::Function
        n = length(x_values)
        if n != length(y_values)
            error("The input data x_values and y_values have different lengths")
        end

        function lagrange_interpolator(x)
            if x in x_values
                # If x is one of the data points, return the corresponding y value
                return y_values[findfirst(isequal(x), x_values)]
            end

            result = 0.0
            for i in 1:n
                term = y_values[i]
                for j in 1:n
                    if i != j
                        term *= (x - x_values[j]) / (x_values[i] - x_values[j])
                    end
                end
                result += term
            end
            return result
        end
        return lagrange_interpolator
    end

    """
        nevilles_method(x_values::Vector{Float64}, y_values::Vector{Float64}, x::Float64)

    Compute the interpolated y value corresponding to a given x value using Neville's method.

    # Arguments
    - `x_values::Vector`: A vector of x values.
    - `y_values::Vector`: A vector of corresponding y values.
    - `x::Number`: The x value for which the corresponding y value is to be computed.

    # Returns
    - `y_interpolated::Number`: The interpolated y value corresponding to the input x value.

    # Errors
    - Throws an error if the input data `x_values` and `y_values` have different lengths.

    # Examples
    ```julia
    x_values = [0, 1, 2, 3, 4]
    y_values = [0, 1, 4, 9, 16]
    x = 2.5
    nevilles_method(x_values, y_values, x) # Output: 6.25
    ```
    """
    function nevilles_method(x_values::Vector{Float64}, y_values::Vector{Float64}, x::Float64)::Float64
        n = length(x_values)
        if n != length(y_values)
            error("the input data x_values and y_values have different legnths")
        end
        Q = map(Float64, y_values) # Initialize the first column of the table with y-values
        for j in 2:n
            for i in 1:n-j+1
                Q[i] = ((x - x_values[i]) * Q[i+1] - (x - x_values[i+j-1]) * Q[i]) / (x_values[i+j-1] - x_values[i])
            end
        end
        
        return Q[1]  # The interpolated value is in the top-right corner of the table
    end

    """
        calculate_divided_differences(x::Vector{Float64}, y::Vector{Float64})

    Calculate the difference table using Newton's divided differences method.

    # Arguments
    - `x::Vector`: The vector of x values.
    - `y::Vector`: The vector of corresponding y values.

    # Returns
    - `F::Matrix`: The difference table computed using Newton's divided differences.

    # Example
    ```julia
    x_values = [0, 1, 2, 3, 4]
    y_values = [0, 1, 4, 9, 16]
    calculate_divided_differences(x_values, y_values) # Output
    # 0.0  1.0  1.0  0.0  0.0
    # 1.0  3.0  1.0  0.0  0.0
    # 4.0  5.0  1.0  0.0  0.0
    # 9.0  7.0  0.0  0.0  0.0
    # 16.0  0.0  0.0  0.0  0.0
    ```
    """
    function calculate_divided_differences(x::Vector{Float64}, y::Vector{Float64})::Matrix{Float64}
        n = length(x)
        F = zeros(n, n)
        if length(y) != n
            error("the input data x_values and y_values have different legnths")
        end
        # Initialize the first column with y values
        F[:, 1] .= y

        # Compute divided differences
        for j in 2:n
            for i in 1:n-j+1
                F[i, j] = (F[i+1, j-1] - F[i, j-1]) / (x[i+j-1] - x[i])
            end
        end
        return F
    end

    """
        update_divided_differences(F::Matrix{Float64}, x::Vector{Float64}, y::Vector{Float64})

    Update the difference table `F` with new data points `x` and `y`.

    # Arguments
    - `F::Matrix`: The difference table computed using Newton's divided differences.
    - `x::Vector`: The new x-values to be added.
    - `y::Vector`: The corresponding y-values for the new x-values.

    # Returns
    - `new_F::Matrix`: The updated difference table with the combined data points.

    # Example
    ```julia
    x_values = [0, 1, 2]
    y_values = [0, 1, 4]
    F = calculate_divided_differences(x_values, y_values)
    x_new_values = [0, 1, 2, 3, 4]
    y_new_values = [0, 1, 4, 9, 16]
    F = update_divided_differences(F, x_new_values, y_new_values)# Output
    # 0.0  1.0  1.0  0.0  0.0
    # 1.0  3.0  1.0  0.0  0.0
    # 4.0  5.0  1.0  0.0  0.0
    # 9.0  7.0  0.0  0.0  0.0
    # 16.0  0.0  0.0  0.0  0.0
    ```
    """
    function update_divided_differences(F::Matrix{Float64}, x::Vector{Float64}, y::Vector{Float64})::Matrix{Float64}
        n = size(F)[1]

        # Determine the number of new data points
        num_new_points = length(x)

        # Extend the difference table if needed
        if num_new_points > n
            F = vcat(F, zeros(num_new_points - n, n))
            F = hcat(F, zeros(num_new_points, num_new_points - n))
        end

        # Add new data points
        for (i, yi) in enumerate(y)
            F[i, 1] = yi  # Update the first column with new y-values
        end

        # Update the last columns of F
        N = num_new_points
        for j in 2:N
            for i in 1:N-j+1
                if (i < n) && (j < n)
                    continue
                end
                F[i, j] = (F[i+1, j-1] - F[i, j-1]) / (x[i+j-1] - x[i])
            end
        end

        return F
    end

    """
        newton_interpolation_from_table(F::Matrix{Float64}, x_values::Vector{Float64}, x::Float64)

    Calculate the interpolation given by the difference table of Newton's divided differences.

    # Arguments
    - `F::Matrix`: The difference table computed using Newton's divided differences.
    - `x_values::Vector`: The vector of x values corresponding to the difference table.
    - `x::Number`: The x value at which to interpolate.

    # Returns
    - `interpolated_value::Number`: The interpolated value at the given x value.

    # Example
    ```julia
    x_values = [0, 1, 2, 3, 4]
    y_values = [0, 1, 4, 9, 16]
    F = calculate_divided_differences(x_values, y_values) 
    newton_interpolation_from_table(F, x_values, 2.5) # Output 6.25
    ```
    """
    function newton_interpolation_from_table(F::Matrix{Float64}, x_values::Vector{Float64}, x::Float64)::Number
        # Extract the first row of the difference table
        coefficients = F[1, :]
        
        n = length(x_values)
        interpolated_value = 0
        
        # Iterate over each coefficient and corresponding x value
        for i in 1:n
            product = 1
            # Compute the product of (x - x_values[j]) for j = 1 to i-1
            for j in 1:i-1
                product *= (x - x_values[j])
            end
            # Add the product multiplied by the coefficient to the total
            interpolated_value += coefficients[i] * product
        end
        
        return interpolated_value
    end

    """
        linear_interpolation(x_values::Vector{Float64}, y_values::Vector{Float64}, x::Float64; sorted=true::Bool)

    Compute linear interpolation of `y` for a given `x` value using the provided data points.

    # Arguments
    - `x_values::Vector`: Vector of x values.
    - `y_values::Vector`: Vector of corresponding y values.
    - `x::Number`: The x value for which the corresponding y value is to be interpolated.
    - `sorted::Bool (optional)`: Boolean indicating whether the input data `x_values` is sorted. Default is `true`.

    # Returns
    - `y::Number`: Interpolated y value corresponding to the input `x` value.

    # Errors
    - Throws an error if the input data `x_values` and `y_values` have different lengths.

    # Example
    ```julia
    x_values = [0, 1, 2, 3, 4]
    y_values = [0, 1, 4, 9, 16]
    x = 2.5
    linear_interpolation(x_values, y_values, x) # Output: 6.5
    ```
    The function performs linear interpolation between adjacent data points 
    (x_values[i], y_values[i]) and (x_values[i+1], y_values[i+1]) to estimate 
    the y value corresponding to the given x value. If x is outside the range of x_values, 
    the function extrapolates using the nearest data points.
    """
    function linear_interpolation(x_values::Vector{Float64}, y_values::Vector{Float64}, x::Float64; sorted=false::Bool)::Number
        n = length(x_values)
        if n != length(y_values)
            error("The input data x_values and y_values have different lengths")
        end

        # Ensure x_values is sorted
        if !sorted
            pairs_of_values = sort(hcat(x_values, y_values), dims=1)
        else
            pairs_of_values = hcat(x_values, y_values)
        end

        # Find the index of the interval where x lies
        index = searchsortedfirst(pairs_of_values[:, 1], x)

        # Perform linear interpolation
        if index == 1
            y = pairs_of_values[1, 2]
        elseif index == n + 1
            y = pairs_of_values[n, 2]
        else
            x0, y0 = pairs_of_values[index, :]
            x1, y1 = pairs_of_values[index + 1, :]
            y = y0 + (y1 - y0) * (x - x0) / (x1 - x0)
        end

        return y
    end
end
