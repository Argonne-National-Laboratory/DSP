.. DSP documentation master file, created by
   sphinx-quickstart on Sat Jul 11 14:57:57 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

DSP --- Decompositions for Structured Programming
=================================================

DSP is an open-source and parallel package that implements decomposition
methods for **structured mixed-integer linear programming problems**.
These are structured optimization problems in the following form:

.. math::

   \min \quad & c^T x + \sum_{s=1}^S q_s^T y_s \\
   \text{s.t.} \quad
   & A x = b \\
   & T_s x + W_s y_s = h_s \quad \forall s = 1, .., S \\
   & \text{mixed-integer } x, y_s

Note that x and y_s are decision variable vectors with dimensions n_1 and n_2,
respectively, A, T_s and W_s are matrices of dimensions m_1 by n_1, m_2 by n_1
and m_2 by n_2, respectively, and c, q_s, b, and h_s are vectors of appropriate
dimensions.

Algorithms in DSP
-----------------

DSP provides **serial** and **parallel** implementations for the four types of algorithms.

Extensive Form Solver
^^^^^^^^^^^^^^^^^^^^^

Dual Decomposition
^^^^^^^^^^^^^^^^^^

Dantzig-Wolfe Decomposition
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Benders Decomposition
^^^^^^^^^^^^^^^^^^^^^

* Dual decomposition (with subgradient method and several bundle methods) with
  **branch-and-bound** procedure
* Benders decomposition

.. toctree::
   :caption: Table of Contents
   :maxdepth: 1

   prerequisites.rst
   installation.rst
   interface.rst
   examples.rst
   acknowledgement.rst
