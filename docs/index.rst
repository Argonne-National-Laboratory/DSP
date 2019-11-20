.. DSP documentation master file, created by
   sphinx-quickstart on Sat Jul 11 14:57:57 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

DSP --- Decompositions for Structured Programming
=================================================

DSP is an open-source and parallel package that implements decomposition methods for **structured mixed-integer programming problems**. These are structured optimization problems in the following form:

.. math::

   \min \quad & c^T x + \sum_{s=1}^S q_s^T y_s \\
   \text{s.t.} \quad
   & A x = b \\
   & T_s x + W_s y_s = h_s \quad \forall s = 1, .., S \\
   & \text{mixed-integer } x, y_s

The current branch ``dsp-bb`` provides **parallel dual decomposition** methods with **branch-and-bound** procedure.

Interfaces
----------

Stand-alone binary
^^^^^^^^^^^^^^^^^^

DSP will be compiled as a binary file runDsp that can read SMPS, MPS with DEC files and solve the problems. 
See more about MPS with DEC file format in http://www.or.rwth-aachen.de/gcg/doc/reader__dec_8h.html and also example in https://github.com/Argonne-National-Laboratory/DSP/tree/dev-coin/examples/mps-dec.

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

JuMP.jl is a algebraic modeling package in Julia for mathematical programming, which can be installed by the following Julia command::

    Pkg.add("JuMP");

The Dsp.jl package provides an interface to JuMP.jl. Dsp.jl can be installed by the Julia command::

    Pkg.clone("https://github.com/kibaekkim/Dsp.jl.git");

MPI.jl is an optional package to run DSP in parallel on high-performance computing machines using MPI library. 
This is an MPI interface package, which can be installed by the following Julia command. MPICH is required for this package::

    Pkg.add("MPI");

Credits
-------

DSP has been developed and is maintained by:

* `Kibaek Kim <http://mcs.anl.gov/~kibaekkim/>`_, Mathematics and Computer Science Division, Argonne National Laboratory.
* `Victor M. Zavala <http://zavalab.engr.wisc.edu/>`_, Department of Chemical and Biological Engineering, University of Wisconsin-Madison.

Publication
-----------

* Kibaek Kim and Brian Dandurand. Scalable Branching on Dual Decomposition of Stochastic Mixed-Integer Programming Problems. Optimization Online, 2018
* Kibaek Kim, Audun Botterud, and Feng Qiu. Temporal Decomposition for Improved Unit Commitment in Power System Production Cost Modeling. IEEE Transactions on Power Systems 33(5), 5276-5287, 2018
* Kibaek Kim and Victor M. Zavala. Algorithmic innovations and software for the dual decomposition method applied to stochastic mixed-integer programs. Mathematical Programming Computation 10(2), 225-266, 2018
* Kibaek Kim, Fan Yang, Victor M. Zavala, and Andrew A. Chien. Data Centers as Dispatchable Loads to Harness Stranded Power. IEEE Transactions on Sustainable Energy 8(1), 208-218, 2017
* Kibaek Kim and Victor M. Zavala. Large-Scale Stochastic Mixed-Integer Programming Algorithms for Power Generation Scheduling. In Alternative Energy Sources and Technologies, 493-512, 2016

Acknowledgements
----------------

This material is based upon work supported by the U.S. Department of Energy, Office of Science, under contract number DE-AC02-06CH11357. 
We thank E. Michael Gertz and Stephen Wright for providing the `OOQP <http://pages.cs.wisc.edu/~swright/ooqp/>`_ software package.

Contents:

.. toctree::
   :maxdepth: 2

   installation.rst
   quickstart.rst

