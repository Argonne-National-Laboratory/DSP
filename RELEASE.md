DSP Release Notes
=================
August 2, 2017
--------------

* DSP can compile with CPLEX

Version 0.3.0 (August 22, 2016)
----------------------------

* New design of the software package
* Implemented the **serial** versions of the methods
* Working with [JuMP](https://github.com/JuliaOpt/JuMP.jl) via new Julia interface [Dsp.jl](https://github.com/kibaekkim/Dsp.jl)

Version 0.2.1 (November 30, 2015)
---------------------------------

* Implemented Benders decomposition with MPI
* Implemented a doubly stabilized bundle method for dual decomposition

Version 0.2.0 (July 20, 2015)
-----------------------------

* Parallelized for cut generation and upper bound evaluataion in dual decomposition
* Reorganized the source code directories
* Released documentation (http://dsp.readthedocs.org/)

Version 0.1.0 (June 18, 2015)
-----------------------------

* Parallel dual decomposition via MPI
* Benders decomposition via OpenMP
* Extensive form solution
* Interfaces for StochJuMP, SMPS and C
