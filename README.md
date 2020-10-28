# DSP

Release: ![GitHub tag (latest SemVer)](https://img.shields.io/github/v/tag/Argonne-National-Laboratory/DSP?label=release&sort=semver)
[![DOI](https://zenodo.org/badge/26612881.svg)](https://zenodo.org/badge/latestdoi/26612881)

Documentation: [![Documentation Status](https://readthedocs.org/projects/dsp/badge/?version=master)](https://dsp.readthedocs.io/?badge=master)

Status: [![Build Status](https://travis-ci.org/Argonne-National-Laboratory/DSP.svg?branch=master)](https://travis-ci.org/Argonne-National-Laboratory/DSP)
[![codecov](https://codecov.io/gh/Argonne-National-Laboratory/DSP/branch/master/graph/badge.svg)](https://codecov.io/gh/Argonne-National-Laboratory/DSP)

--------

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

**Problem Types:**
* Two-stage stochastic mixed-integer linear programs
* Distributionally robust stochastic mixed-integer linear programs
* Structured mixed-integer linear programs

**Problem Input Formats:**
* SMPS file format for stochastic programs (`.dro` optionally for distributionally robust)
* MPS and DEC files for generic block-structured optimization problems
* Julia modeling package [DSPopt.jl](https://github.com/kibaekkim/DSPopt.jl)

## Installation

```
git clone --recursive https://github.com/Argonne-National-Laboratory/DSP.git
```

## Contributors

* [Kibaek Kim](https://kibaekkim.github.io/), Mathematics and Computer Science Division, Argonne National Laboratory.
* [Victor M. Zavala](http://zavalab.engr.wisc.edu/), Department of Chemical and Biological Engineering, University of Wisconsin-Madison.
* Christian Tjandraatmadja, Google Research.
* Yingqiu Zhang, Industrial and Systems Engineering, Virginia Tech.

## Key Publications

* Kibaek Kim. "[Dual Decomposition of Two-Stage Distributionally Robust Mixed-Integer Programming under the Wasserstein Ambiguity Set](http://www.optimization-online.org/DB_HTML/2020/04/7723.pdf)" Optimization Online, 2020
* Kibaek Kim and Briand Dandurand. "[Scalable Branching on Dual Decomposition of Stochastic Mixed-Integer Programming Problems](http://www.optimization-online.org/DB_HTML/2018/10/6867.html)" Mathematical Programming Computation (to appear), 2020
* Kibaek Kim, Cosmin Petra, and Victor Zavala. "[An Asynchronous Bundle-Trust-Region Method for Dual Decomposition of Stochastic Mixed-Integer Programming](https://epubs.siam.org/doi/abs/10.1137/17M1148189)" SIAM Journal on Optimization 29(1), 2019
* Kibaek Kim and Victor M. Zavala. "[Algorithmic innovations and software for the dual decomposition method applied to stochastic mixed-integer programs](https://link.springer.com/article/10.1007/s12532-017-0128-z)" Mathematical Programming Computation 10(2), 2017


## Acknowledgements

This material is based upon work supported by the U.S. Department of Energy, Office of Science, under contract number DE-AC02-06CH11357.
