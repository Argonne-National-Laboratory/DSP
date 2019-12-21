# DSP
[![DOI](https://zenodo.org/badge/26612881.svg)](https://zenodo.org/badge/latestdoi/26612881)
[![Documentation Status](https://readthedocs.org/projects/dsp/badge/?version=master)](https://dsp.readthedocs.io/?badge=master)

DSP is an open-source and parallel package that implements decomposition methods for **structured mixed-integer programming** problems. These are structured optimization problems in the following form:

        minimize   c^T x + \sum_{s=1}^S q_s^T y_s
        subject to   A x                              = b
                   T_s x +                    W_s y_s = h_s for s = 1, .., S
                   some x, y_s are integers

where x and y_s are decision variable vectors with dimensions n_1 and n_2, respectively, A, T_s and W_s are matrices of dimensions m_1 by n_1, m_2 by n_1 and m_2 by n_2, respectively, and c, q_s, b, and h_s are vectors of appropriate dimensions.

**DSP Solution Methods:**
* Extensive form solver (global solver)
* Serial/parallel dual decomposition (dual bounding solver)
* Serial/parallel Dantzig-Wolfe decomposition (global solver)
* Serial/parallel Benders decomposition

**Problem Input Formats:**
* SMPS file format for stochastic programs
* MPS and DEC files for generic block-structured optimization problems
* Julia modeling package [Dsp.jl](https://github.com/kibaekkim/Dsp.jl)

## Download

You can clone this repository in your preferred directory by typing:
```bash
git clone --recursive https://github.com/Argonne-National-Laboratory/DSP.git
```
or
```bash
git clone https://github.com/Argonne-National-Laboratory/DSP.git
cd DSP
git submodule update --init --recursive
```

## Documentation

The package documentation is available in [Readthedocs](https://dsp.readthedocs.io/?badge=master).

## Credits

DSP has been developed and is maintained by:
* [Kibaek Kim](http://mcs.anl.gov/~kibaekkim/), Mathematics and Computer Science Division, Argonne National Laboratory.
* [Victor M. Zavala](http://zavalab.engr.wisc.edu/), Department of Chemical and Biological Engineering, University of Wisconsin-Madison.


## Key Publications

* Kibaek Kim and Briand Dandurand. "[Scalable Branching on Dual Decomposition of Stochastic Mixed-Integer Programming Problems](http://www.optimization-online.org/DB_HTML/2018/10/6867.html)" Optimization Online, 2018
* Kibaek Kim, Cosmin Petra, and Victor Zavala. "[An Asynchronous Bundle-Trust-Region Method for Dual Decomposition of Stochastic Mixed-Integer Programming](https://epubs.siam.org/doi/abs/10.1137/17M1148189)" SIAM Journal on Optimization 29(1), 2019
* Kibaek Kim and Victor M. Zavala. "[Algorithmic innovations and software for the dual decomposition method applied to stochastic mixed-integer programs](https://link.springer.com/article/10.1007/s12532-017-0128-z)" Mathematical Programming Computation 10(2), 2017


## Acknowledgements

This material is based upon work supported by the U.S. Department of Energy, Office of Science, under contract number DE-AC02-06CH11357. We gratefully acknowledge the computing resources provided on *Blues*, a high-performance computing cluster operated by the Laboratory Computing Resource Center at Argonne National Laboratory. We thank E. Michael Gertz and Stephen Wright for providing the [OOQP](http://pages.cs.wisc.edu/~swright/ooqp/) software package.


[![Analytics](https://ga-beacon.appspot.com/UA-64449559-2/DSP/README.md)](https://github.com/igrigorik/ga-beacon)
