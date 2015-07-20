DSP Overview
------------

This section provides an introduction to DSP: Decompositions for Structured Programming. DSP is an object-oriented open-source software package written in C++ for solving structured programming problems. DSP implements different solution methods for solving *stochastic mixed-integer programming (SMIP)* problems. The constraint matrix of SMIP problem has a *block-angular structure* as follows:

.. math::

   \left[\begin{matrix} A & & & & \\
   T_1 & W_1 & & & \\
   T_2 & & W_2 & & \\
   \vdots & & & \ddots & \\
   T_S & & & & W_S\end{matrix}\right] \left[\begin{matrix} x \\
   y_1 \\
   y_2 \\
   \vdots \\
   y_S \end{matrix}\right] = \left[\begin{matrix} b \\
   h_1 \\
   h_2 \\
   \vdots \\
   h_S \end{matrix}\right]

The block-angular structure leads to different decomposition methods that enable large-scale optimization problems (e.g., stochastic programming problmes) to be efficiently solved or approximated particularly on high-performance computing systems. DPS exploits the block-angular structure of SMIP and provides the following **methods**:

* Parallel dual decomposition methods

  * Different dual update methods:

    * Interior-point cutting-plane method with early termination criteria
    * Cutting-plane method
    * Subgradient method

  * Effectively eliminating infeasible first-stage solution and thus aids to obtain upper bounds.
  * Significant reductions in number of iterations and solution time, as compared with standard dual decomposition methods.
  * Parallel runs on high performance computing systems using MPI
  * Scaling in number of scenarios with many computing cores

* Parallel Benders decomposition methods

  * Solving stochastic program with first-stage mixed-integer variables
  * Parallel runs with many cores using OpenMP

* Subgradient method
* Extensive form solution

.. note:: The methods implemented in DSP are formally described in :cite:`kim2015algorithmic`.

**Additional features** of DSP include:

* Interface for the StochJuMP.jl package in Julia

  * `StochJuMP.jl <https://github.com/kibaekkim/StochJuMP.jl>`_ -- a scalable algebraic modeling package for stochastic programming in Julia
  * Extension of the JuMP.jl package, one of the fastest optimization modeling package in Julia
  * User can run DSP in parallel via Julia without coding any MPI codes.

* Solver independence

  * Any external solvers can be used in DSP through the SolverInterface class
  * Current support: SCIP, SOPLEX, CLP, OOQP

Decomposition of SMIP
^^^^^^^^^^^^^^^^^^^^^

Stochastic mixed-integer programming (SMIP) is an optimization modeling framework for problems that involve mixed-integer decision variables and uncertain parameters. We considers the following two-stage SMIP problem:

.. math::

   z := \min_x \left\{ c^T x + \mathbb{E}_\xi\left[Q(x,\xi)\right] :\; A x = b,\; x \in X \right\},
   
where the scenario recourse function

.. math::

   Q(x,\xi) := \min_y\; \{q(\xi)^T y :\; W(\xi) y = h(\xi) - T(\xi) x, \; y \in Y\}.

We assume that the random parameter :math:`\xi` follows a discrete distribution with finite support :math:`\{\xi^1, \dots, \xi^S\}` and corresponding probabilities :math:`p_1,\dots,p_S` (continuous distributions can be handled by using a sample-average approximation).

.. note:: The sets :math:`X \subseteq \mathbb{R}_+^{n_1}` and :math:`Y \subseteq \mathbb{R}_+^{n_2}` represent integer restrictions on a subset of the decision variables :math:`x` and :math:`y`, respectively.

The first-stage problem data comprises :math:`A \in \mathbb{R}^{m_1\times n_1}, b\in \mathbb{R}^{m_1}`, and :math:`c\in \mathbb{R}^{n_1}`.  The second-stage data are given by :math:`T(\xi^s)\in \mathbb{R}^{m_2 \times n_1}, W(\xi^s)\in \mathbb{R}^{m_2 \times n_2}, h(\xi^s) \in \mathbb{R}^{m_2}`, and :math:`q(\xi^s) \in \mathbb{R}^{n_2}`. For simplicity, we use the notation :math:`(T_s, W_s, h_s, q_s)` for :math:`s\in \mathcal{S} := \{1,\dots,S\}`.

Dual Decomposition
******************

Dual decomposition considers SMIP of the extensive form

.. math::
   :label: DDEF

   z = \min_{x_s,y_s} \quad & \sum_{s\in\mathcal{S}} p_s \left( c^T x_s + q_s^T y_s \right) \\
   \text{s.t.} \quad
   & \sum_{s\in\mathcal{S}} H_s x_s = 0 \quad \textit{(Nonanticipativity constraint)} \\
   & (x_s, y_s) \in G_s, \quad \forall s\in\mathcal{S}

where the scenario feasibility set is defined as

.. math::

  G_s := \{(x_s, y_s) \;:\; A x_s = b, \; T_s x_s + W_s y_s = h_s, \; x_s \in X, \; y_s \in Y\},

the *nonanticipativity* constraints represents the equations :math:`x_1 = x_S` and :math:`x_s = x_{s-1}` for :math:`s=2,\dots,S`, and :math:`H_s` is a suitable :math:`S\cdot n_1 \times n_1` matrix. 

.. note:: SMIP may not have relatively complete recourse. Without this property, there can exist :math:`(\hat x,\hat y)` such that :math:`(\hat x,\hat y) \in G_s` and :math:`(\hat x,\hat y) \notin G_{s'}` for :math:`s\neq s'`. 

We apply a Lagrangian relaxation of these constraints to obtain the Lagrangian dual function of :eq:`DDEF`: 

.. math::

   D(\lambda) := \min_{x_s,y_s} \left\{ \sum_{s\in\mathcal{S}} L_s(x_s,y_s,\lambda) : (x_s,y_s) \in G_s,\; \forall s\in\mathcal{S} \right\},

where 

.. math::

   L_s(x_s,y_s,\lambda) := p_s \left(c^T x_s + q_s^T y_s\right) + \lambda^T (H_s x_s). 

For fixed :math:`\lambda`, the Lagrangian dual function can be decomposed as

.. math::

   D(\lambda) = \displaystyle\sum_{s\in\mathcal{S}} D_s(\lambda),

where

.. math::

   D_s(\lambda) := \min_{x_s,y_s} \left\{ L_s(x_s,y_s,\lambda) : (x_s,y_s) \in G_s \right\}.

We thus seek to obtain the best lower bound for :eq:`DDEF` by solving the maximization problem (the Lagrangian dual problem):

.. math::

  z_\text{LD} := \max_{\lambda} \sum_{s\in\mathcal{S}} D_s(\lambda).


Benders Decomposition
*********************

Benders decomposition considers SMIP of the form

.. math::
  \min \quad & c^T x + \sum_{s\in \widetilde{\mathcal{S}}} p_s q_s^T y_s + \theta\\
   \text{s.t.} \quad
   & A x = b, \; x \in X \\
   & T_s x + W_s y_s = h_s, \; y_s \in Y, \; \forall s\in \widetilde{\mathcal{S}}, \\
   & \theta \geq \sum_{s\in S\backslash\widetilde{\mathcal{S}}} p_s Q(x,\omega_s),

where :math:`\widetilde{\mathcal{S}}` is a subset of :math:`\mathcal{S}` given by user. The method performs outer-approximation of the recourse function :math:`Q(x,\omega_s)` by iteratively adding a set of linear inequalities. DSP implements a standard Benders decomposition method for solving SMIP problems with first-stage mixed-integer variables.

.. warning:: The second-stage integrality is relaxed in DSP solution.

Design of the DSP Development Framework
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The software design is object-oriented and implemented in C++. It consists of *Model* classes and *Solver* classes for handling optimization models and scenario data.

Model Classes
*************

An abstract *Model* class is designed to define a generic optimization model data structure. The *StoModel* class defines the data structure for generic stochastic programs, including two-stage stochastic programs and multistage stochastic programs. The underlying data structure of *StoModel* partially follows the SMPS format. The class also defines core functions for problem decomposition. The *TssModel* class derived defines the member variables and functions specific to two-stage stochastic programs and decompositions. Following the design of the model classes, users are able to derive new classes for their own purposes and efficiently manage model structure provided from several interfaces (e.g., StochJuMP and SMPS).

Solver Classes
**************

An abstract *Solver* class is designed to provide different algorithms for solving stochastic programming problems defined in the *Model* class. DSP implements the *TssSolver* class to define solvers specific to two-stage stochastic programs. From the *TssSolver* class, three classes are derived for each method: *TssDe*, *TssBd*, and *TssDd*.

* The *TssDe* class implements a wrapper of external solvers to solve the extensive form of two-stage stochastic programs. The extensive form is constructed and provided by the *TssModel* class.
* The *TssBd* class implements a Benders decomposition method for solving two-stage stochastic programs with continuous recourse. A proper decomposition of the model is performed and provided by the *TssModel* class, while the second-stage integrality restriction is automatically relaxed. Depending on parameters provided, *TssModel* can make a different form of the problem decomposition for *TssBd*. For example, the user can specify the number of cuts added per iteration, which determines the number of auxiliary variables in the master problem of Benders decomposition. Moreover, the Benders master can be augmented for a subset :math:`\widetilde{\mathcal{S}}` of scenarios.
* The *TssDd* class implements the proposed dual decomposition method for solving two-stage stochastic programs with mixed-integer recourse. For this method, an abstract *TssDdMaster* class is designed to implement methods for updating the dual variables. The subgradient method and the cutting-plane method are implemented in such derived classes. Moreover, a subclass derived from the *TssBd* is reused for implementing the Benders-type cutting-plane procedure for the subproblems. An :math:`l_\infty`-norm trust region is also applied in order to stabilize the cutting-plane method. The rule of updating the trust region follows that proposed in :cite:`linderoth2003decomposition`. Users can also implement their own method for updating the dual variables.

External Solver Interface Classes
*********************************

DSP uses external MIP solvers to solve subproblems under different decomposition methods. The *SolverInterface* class is an abstract class to create interfaces to the decomposition methods implemented. Several classes are derived from the abstract class in order to support specific external solvers. The current implementation supports the following external optimization solvers:

* Linear programming solvers

  * Clp :cite:`clp`
  * SoPlex :cite:`Wunderling1996`
  * OOQP :cite:`gertz2003object`

* Mixed-integer programming solver

  * SCIP :cite:`achterberg2009scip`

Users familiar with the COIN-OR Open Solver Interface :cite:`saltzman2004coin` should easily be able to use the *SolverInterfaceOsi* class to derive classes for other solvers (e.g., CPLEX :cite:`cplex`, Gurobi :cite:`gurobi`).

Parallelization
***************

The proposed dual decomposition method can be run on distributed memory and on shared memory computing systems with multiple cores. The implementation protocol is MPI.  In a distributed memory environment, the scenario data and corresponding Lagrangian subproblems are distributed to multiple processors based on scenario indices. The root processor updates the Lagrangian multipliers and solves a subset of the subproblems. When solving the subproblems in distributed computing nodes, subproblem solutions and the dual variables must be communicated with the root processor. In addition, each computing node communicates the primal first-stage solutions and the valid inequalities generated for a subproblem with the rest of the nodes.

.. bibliography:: overview.bib

