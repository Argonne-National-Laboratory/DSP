.. DSP documentation master file, created by
   sphinx-quickstart on Sat Jul 11 14:57:57 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

DSP --- Decompositions for Structured Programming
=================================================

DSP is an open-source and parallel package that implements decomposition methods for **structured mixed-integer programming problems**. 
These are structured optimization problems in the following form:

.. math::

   \min \quad & c^T x + \sum_{s=1}^S q_s^T y_s \\
   \text{s.t.} \quad
   & A x = b \\
   & T_s x + W_s y_s = h_s \quad \forall s = 1, .., S \\
   & \text{mixed-integer } x, y_s

.. note::
    This version comes from the ``dsp-bb`` branch of the Github repository for ``DSP``. 
    This version only implements **parallel dual decomposition** methods with **branch-and-bound** procedure.
    The ``master`` branch of the Github repository for ``DSP``  supports more algorithms.

Contents:

.. toctree::
   :caption: Table of Contents
   :maxdepth: 1

   installation.rst
   interface.rst
   examples.rst
   acknowledgement.rst
