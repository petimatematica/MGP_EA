# Gradient Projected Method with GPA1 and GPA2

## GPA1.jl
This file contains the implementation of the Projected Gradient Method equipped with Armijo's linesearch along feasible directions.

### Function GPA1
Armijo linesearch along feasible directions.

#### Arguments

- **$\boldsymbol{x}$ (Vector):** Point where the linear search will applied
- **$f$ (Function):** The objective function
- **$\nabla f$ (Function):** The gradient of the objective function
- **$\sigma$ (FLoat64):** The parameter of Armijo linesearch
- **min_step (FLoat64):** The minimum stepsize parameter acceptable
- **$\gamma$ _start (FLoat64):** Initial stepsize
- **$zk$ (Vector):** $P_C(x_k - \beta_k\nabla f(x_k))$

#### Outputs

- **$\boldsymbol{\gamma}$ (Float64):** The stepsize parameter in iterarion k
- **ierror (Int64):** The error of linesearch (0 - no error, 1 - step length too small)
- **evalf (Int64):** Number of function evaluations

### Function method1
Gradient Projected Method with Armijo linesearch along feasible directions.

#### Arguments

- **$\boldsymbol{x}$ (Vector):** Initial point
- **$f$ (Function):** The objective function
- **$\nabla f$ (Function):** The gradient of the objective function
- **$\varepsilon$ (FLoat64):** The convergence tolerance parameter
- **max_iter (Int64):** Maximum number of iterations allowed 
- **GPA1 (Function):** The linesearch function
- **projection (Function):** The projection function

#### Outputs

- **$\boldsymbol{x}$ (Vector):** Minimizer point
- **$\boldsymbol{f(x)}$ (Float64):** Objective function evaluated in the minimizer
- **info (DataFrame):** Data frame with information about the minimization process
- **et (Float64):** Time spent solving the problem
- **ierror (Int64):** Error returned (0 - no error, 1 - step length too small, 2 - Maximum of iterations was achieved)
- **seqx )Matrix{Vector{Float64}}):** Matrix with the terms of the sequence $\{x_k\}$
- **sum(evalf_\gamma) (Int64):** Total number of function evaluations 

## GPA2.jl
This file contains the implementation of the Projected Gradient Method equipped with Armijo linesearch along the projection arc.

### Function GPA2
Armijo linesearch along the projection arc.

#### Arguments

- **$\boldsymbol{x}$ (Vector):** Point where the linear search will be carried out
- **$f$ (Function):** The objective function
- **$\nabla f$ (Function):** The gradient of the objective function
- **$\sigma$ (FLoat64):** The parameter of Armijo linesearch
- **min_step (FLoat64):** The minimum stepsize parameter acceptable
- **$\beta$ _start (FLoat64):** Initial stepsize

#### Outputs

- **$\boldsymbol{\beta}$ (Float64):** The stepsize parameter in iterarion k
- **zkj (Vector):** $P_C(x_k - \beta_k\nabla f(x_k))$
- **ierror (Int64):** The error of linesearch (0 - no error, 1 - step length too small)
- **evalf (Int64):** Number of function evaluations

### Function method2
Gradient Projected Method with Armijo linesearch along the projection arc.

#### Arguments

- **$\boldsymbol{x}$ (Vector):** Initial point
- **$f$ (Function):** The objective function
- **$\nabla f$ (Function):** The gradient of the objective function
- **$\varepsilon$ (FLoat64):** The convergence tolerance parameter
- **max_iter (Int64):** Maximum number of iterations allowed 
- **GPA2 (Function):** The linesearch function
- **projection (Function):** The projection function

#### Outputs

- **$\boldsymbol{x}$ (Vector):** Minimizer point
- **$\boldsymbol{f(x)}$ (Float64):** Objective function evaluated in the minimizer
- **info (DataFrame):** Data frame with information about the minimization process
- **et (Float64):** Time spent solving the problem
- **ierror (Int64):** Error returned (0 - no error, 1 - step length too small, 2 - Maximum of iterations was achieved)
- **seqx )Matrix{Vector{Float64}}):** Matrix with the terms of the sequence $\{x_k\}$
- **sum(evalf_\gamma) (Int64):** Total number of function evaluations 


## main.jl
This Julia script performs a comprehensive test on the Dixon-Price function using the projected gradient method and generates a data frame with the data obtained during the optimization process.

### Code overview
The script is organized into the following sections:

1. **Imports and Inclusions:** The script imports the necessary packages (`LinearALgebra`, `DataFrames` and `Random`) and includes the essential Julia files: `GPA1.jl` and `GPA2.jl`.

2. **Dixon-Price function:**
   - **$f$ (Function):** The Dixon-Price function
   - **$\nabla f$ (Function):** The gradient of the Dixon-Price function

3. **Projection functions for diferent feasible sets:**
   - **projection1:** $C = \{x \in \mathbb{R}^n; |x-c| \le 20 \}$, onde $c = (-50,0,0,0,\ldots)$.
   - **projection2:** $C = \{x \in \mathbb{R}^n ; \langle a, x \rangle = 0\}$, onde $a = (1,1,\cdots,1)$
   - **projection3:** $C = \{x \in \mathbb{R}^n ; \langle a, x \rangle \le 0\}$, onde $a = (-1,-1,\cdots,-1)$
   - **projection4:** $C = \{x \in \mathbb{R}^n ; Ax = (0,0)\}$, onde $A \in \mathbb{R}(2,n)$ tem posto $2$

4. **Parameters:**
Definition of all parameters to be used in the optimization process.

6. **Optimization Cycle:**
The script performs the gradient method optimization for the chosen case and show the informations at the end of process.
