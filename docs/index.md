# What is DSP?

DSP (**D**ecomposition of **S**tructured **P**rograms) is an open-source and parallel package that implements decomposition
methods for **structured mixed-integer linear programming problems**.

Examples of the structured programming problem include:

- distributionally robust optimization problem (with Wasserstein ambiguity set)
- stochastic programming problem
- long-term planning problem
- network optimization problem

## Algorithms in DSP

DSP provides **serial** and **parallel** implementations for the four types of algorithms.

### Dantzig-Wolfe Decomposition

This solves the dual of the dual decomposition.
The key difference from the dual decomposition is the implementation of branch-and-bound method on top of the decomposition,
which allows to find a global optiaml solution.

For more technical details, please refer the following papers.

1. Kibaek Kim and Brian Dandurand. "Scalable Branching on Dual Decomposition of Stochastic Mixed-Integer Programming Problems" Mathematical Programming Computation (to appear), 2020
1. Kibaek Kim, Audun Botterud, and Feng Qiu. "Temporal Decomposition for Improved Unit Commitment in Power System Production Cost Modeling" IEEE Transactions on Power Systems. 2018

### Dual Decomposition

This solves the Lagrangian relaxation of the constraints that couple blocks of the structure.
Hence, the algorithm finds a Lagrangian dual bound that can be strictly lower than the optimal objective value for nonconvex problems (e.g., mixed-integer programs).
For two-stage stochastic programs, each block is represented by a scenario subproblem.
Stochastic programs can be read from [SMPS](https://ieeexplore.ieee.org/abstract/document/8142546) files or our C API functions (or Julia package [DSPopt.jl](https://github.com/kibaekkim/DSPopt.jl)).
For generic optimization problems, each block can be specified by [mps-dec](https://gcg.or.rwth-aachen.de/doc/reader__dec_8h.html) files or our C API functions (or Julia package DSPopt.jl).

> NOTE: A **distributionally robust** variant of dual decomposition is available for the problems given in SMPS files and `.dro` file.

For more technical details, please refer the following papers.

1. Kibaek Kim. "Dual Decomposition of Two-Stage Distributionally Robust Mixed-Integer Programming under the Wasserstein Ambiguity Set" Optimization Online, 2020
1. Kibaek Kim, Cosmin G. Petra, and Victor M. Zavala. "An Asynchronous Bundle-Trust-Region Method for Dual Decomposition of Stochastic Mixed-Integer Programming" SIAM Journal on Optimization, 2019
1. Kibaek Kim, Mihai Anitescu, and Victor M. Zavala. "An Asynchronous Decomposition Algorithm for Security Constrained Unit Commitment under Contingency Events" Proceedings in Power System Computation Conference, 2018
1. Kibaek Kim and Victor M. Zavala. "Algorithmic Innovations and Software for the Dual Decomposition Method applied to Stochastic Mixed-Integer Programs" Mathematical Programming Computation, 2017

### Integer Benders Decomposition

This solves the integer Benders decomposition of stochastic programming problems with the lower bound initialized by the dual decomposition.
Integer Benders decomposition runs when the first stage is a pure-binary program and the second stage is a mixed-integer program.

> NOTE: A **distributionally robust** variant of Benders decomposition is available for the problems given in SMPS files and `.dro` file.

> Note: When the first stage is not a pure-binary program, any second-stage integer variables will be relaxed.

> Note: This algorithm requires SCIP and is available only for stochastic programs.

### Extensive Form Solver

This reformulates the problem into one large optimization problem and uses available external solver.
Parallel computing is available only through the external solver (i.e., multi-threading).

> Note: There is no MPI parallelism for this algorithm.
