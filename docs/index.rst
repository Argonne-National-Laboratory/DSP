.. DSP documentation master file, created by
   sphinx-quickstart on Sat Jul 11 14:57:57 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

DSP --- Decompositions for Structured Programming
=================================================

DSP is an object-oriented open-source software package written in C++ for solving structured programming problems such as stochastic mixed-integer programming (SMIP). The current version of DSP has implemented decomposition methods for solving SMIP problems:

* Dual decomposition -- DSP has implemented new dual decomposition methods

  * **Interior-point cutting-plane method with early termination criteria**
  * **Eliminating infeasible first-stage solutions**
  * Cutting-plane method
  * Subgradient method

* Benders decomposition
* Extensive form solution

**Credits**

DSP has been developed and is maintained by:

* `Kibaek Kim <http://mcs.anl.gov/~kibaekkim/>`_, Mathematics and Computer Science Division, Argonne National Laboratory.
* `Victor M. Zavala <http://zavalab.engr.wisc.edu/>`_, Department of Chemical and Biological Engineering, University of Wisconsin-Madison.

**Publication**

* Kibaek Kim and Victor M. Zavala. `Algorithmic innovations and software for the dual decomposition method applied to stochastic mixed-integer programs`_ Mathematical Programming Computation 10(2), 225-266, 2018
* Kibaek Kim, Fan Yang, Victor M. Zavala, and Andrew A. Chien. `Data Centers as Dispatchable Loads to Harness Stranded Power'_ IEEE Transactions on Sustainable Energy 8(1), 208-218, 2017
* Kibaek Kim and Victor M. Zavala. `Large-Scale Stochastic Mixed-Integer Programming Algorithms for Power Generation Scheduling'_ In Alternative Energy Sources and Technologies, 493-512, 2016

**Acknowledgements**

This material is based upon work supported by the U.S. Department of Energy, Office of Science, under contract number DE-AC02-06CH11357. We gratefully acknowledge the computing resources provided on *Blues*, a high-performance computing cluster operated by the Laboratory Computing Resource Center at Argonne National Laboratory. We thank E. Michael Gertz and Stephen Wright for providing the `OOQP <http://pages.cs.wisc.edu/~swright/ooqp/>`_ software package.

Contents:

.. toctree::
   :maxdepth: 2

   installation.rst
   quickstart.rst

