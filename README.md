DSP
===
DSP is an open-source optimization solver that implements decomposition methods for stochastic mixed-integer programming (SMIP) problems. These are structured optimization problems considering uncertain scenario realization s with probability p_s in the following form:

        minimize   c^T x + \sum_{s=1}^S p_s q_s^T y_s
        subject to   A x                              = b
                   T_s x +                    W_s y_s = h_s for s = 1, .., S
                   some x, y_s are integers

where x and y_s are decision variable vectors with dimensions n_1 and n_2, respectively, A, T_s and W_s are matrices of dimensions m_1 by n_1, m_2 by n_1 and m_2 by n_2, respectively, and c, q_s, b, and h_s are vectors of appropriate dimensions.

The problem structures allow decomposition approaches for solving the problem in parallel computing system. DSP provides parallel implementations for decomposition methods (Benders decomposition and dual decomposition). The methods can be run on cluster as well as desktop computers.

Getting started
---------------
The file INSTALL.md, included in this directory, is the installation guide for the DSP package.

Julia interface
---------------
DSP provides an interface to an algebraic modeling package StochJuMP based on Julia. The interface package can be installed by the Julia command
```julia
Pkg.clone("https://github.com/kibaekkim/DSP.jl.git");
```
and the examples are located at
```
./examples/julia
```
Once the Julia interface has been properly installed, users can model problems by using StochJuMP, a algebraic modeling package in Julia for stochastic programming, which can be installed by the following Julia command:
```julia
Pkg.clone("https://github.com/kibaekkim/StochJuMP.jl.git");
```
Then, users can call a decomposition solver for the solution. This also  can be run on cluster by using MPI.

Examples
--------
The examples/ directory contains examples to demonstrate the use of DSP. 

Contents of the distribution
----------------------------
The DSP distribution contains the following top-level subdirectories:
* examples/ -- Example problems and programs.
* extra/ -- External packages required and used for DSP. These include MA27, OOQP, SCIP and SMI.
* src/ -- Source files for the DSP distribution.

Credits
-------
DSP is maintained by:
* Kibaek Kim, Mathematics and Computer Science Division, Argonne National Laboratory.
* Victor M. Zavala, Mathematics and Computer Science Division, Argonne National Laboratory.

