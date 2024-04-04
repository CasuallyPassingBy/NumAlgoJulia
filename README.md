# Numerical Analysis Algorithms in Julia

Welcome to the Numerical Analysis Algorithms repository in Julia! This repository hosts a collection of efficient algorithms for numerical analysis tasks, specifically focusing on finding zeros/fixed points of functions and data interpolation.

## Features

- **Zero Finding Algorithms**: Implementations of various methods for finding zeros and fixed points of functions, including but not limited to:
  - Bisection method
  - Newton-Raphson method
  - Secant method
  - Fixed-point iteration method
- **Interpolation Techniques**: Tools for interpolating data using methods such as:
  - Linear interpolation
  - Polynomial interpolation (Lagrange, Newton)
  - Cubic spline interpolation

## Getting Started

To get started with using the algorithms in this repository, follow these steps:

1. Clone the repository to your local machine:

   ```bash
   git clone https://github.com/CasuallyPassingBy/NumAlgoJulia.git
   ```
1. Navigate to the repository directory:
    ```bash
     cd NumAlgoJulia
     ```
1. Explore the src directory to find the implementations of the numerical analysis algorithms.



Use the algorithms in your Julia projects by importing the necessary modules.

## Usage
Here's a basic example demonstrating how to use the zero finding algorithms:
```julia
using .FindingZeros
# Function you want to find the zero of
f(x) = x^2 - 2
#Find the approximate zero between 1 and 2 using the Bisection method
approx_zero = bisection(f, 1, 2)
```
And here's how to perform linear interpolation on a dataset:
```julia
using .Interpolating
# Sample data points
x_values = [0, 1, 2, 3, 4]
y_values = [0, 1, 4, 9, 16]

# Interpolate at x = 2.5 using Neville's method
interp_value_neville = nevilles_method(x_values, y_values, 2.5)
```
## Contributing

Contributions to this repository are welcome! If you have implemented additional numerical analysis algorithms or have suggestions for improvements, feel free to open an issue or submit a pull request.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

