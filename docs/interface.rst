Interfaces
----------

Stand-alone binary
^^^^^^^^^^^^^^^^^^

DSP will be compiled as a binary file runDsp that can read 

* SMPS file format (https://ieeexplore.ieee.org/abstract/document/8142546)
* MPS file format with additional DEC file and solve the problems. See more about MPS with DEC file format in http://www.or.rwth-aachen.de/gcg/doc/reader__dec_8h.html and also example in https://github.com/Argonne-National-Laboratory/DSP/tree/dev-coin/examples/mps-dec.

Julia Interface
^^^^^^^^^^^^^^^

DSP can use Julia as a modeling interface for the computational experiments. 
Julia is a high level dynamic programming language for technical computing, 
with syntax that is familiar to users of other technical computing environments such as MATLAB and Python. 
Julia can be downloaded from http://julialang.org/downloads/. 
If Julia is successfully installed on your machine, then you can start the Julia command-line tool by typing::

    julia

Now we need to install four Julia packages required to run DSP in the Julia environment. 
The packages should be installed in the Julia command-line tool. Please update the package list of Julia by typing::

    Pkg.update();

The following packages should be installed.

* JuMP.jl is a algebraic modeling package in Julia for mathematical programming, which can be installed by the following Julia command::

    Pkg.add("JuMP");

* The Dsp.jl package provides an interface to JuMP.jl. Dsp.jl can be installed by the Julia command::

    Pkg.clone("https://github.com/kibaekkim/Dsp.jl.git");

MPI.jl is an *optional* package to run DSP in parallel on high-performance computing machines using MPI library. 
This is an MPI interface package, which can be installed by the following Julia command. MPICH is required for this package::

    Pkg.add("MPI");