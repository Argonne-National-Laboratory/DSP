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
If you have all the prerequisite packages installed on your system (see the [INSTALL.md](./INSTALL.md)), then you need to go to the root directory of DSP and type
```cmake
cmake .
```
to configure OOQP. If you wish to install the package in a more permanent location, you may then type
```
make install
```
External packages (MA27, OOQP, SCIP Optimization Suite, Smi) used in DSP are built automatically. A shared object is installed in ./lib directory. Once the installation has been successfully done, you need to set an environment variable DSP_INC and (DY)LD_LIBRARY_PATH.
For Linux,
```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<DSP source path>/lib
```
For Mac,
```
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:<DSP source path>/lib
```

Julia interface
---------------
If Julia is available on your machine, you can model problems by using StochJuMP, a algebraic modeling package in Julia for stochastic programming, which can be installed by the following Julia command:
```julia
Pkg.clone("https://github.com/kibaekkim/StochJuMP.jl.git");
```
DSP provides an interface to StochJuMP. The DSP interface package can be installed by the Julia command
```julia
Pkg.clone("https://github.com/kibaekkim/DSPsolver.jl.git");
```
and the examples are located at
```
./examples/julia
```
Then, users can call a decomposition solver for the solution. This also can be run on cluster by using MPI.

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

