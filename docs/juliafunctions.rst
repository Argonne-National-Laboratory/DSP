Julia Interface Functions
-------------------------

We describe Julia interface functions.

Model and Solution Functions
++++++++++++++++++++++++++++

.. function:: readSmps(filename)

   Read SMPS files from files ``filename``.cor, ``filename``.tim and ``filename``.sto.

.. function:: loadProblem(model::JuMP.Model)

   Read model from stochastic model ``model::JuMP.Model``.

.. function:: freeModel()

   Release all the internal model data. You may call this only if you need to read a new model.

.. function:: solve([solver]) 

   Solve the model by using solver type ``solver``. Possible values of ``solver`` are ``DSP_SOLVER_DD``, ``DSP_SOLVER_BD`` AND ``DSP_SOLVER_DE``. Default value is ``DSP_SOLVER_DD``.

.. function:: evaluateSolution(solution::Vector{Cdouble})

   Evaluate the objective value of the first-stage solution specified in ``solution``.

Get Functions
+++++++++++++

Common functions
****************

.. function:: getNumScenarios()

   Get number of scenarios.

.. function:: getNumRows()

   Get number of rows in the extensive form.

.. function:: getNumRows(stage)

   :param stage: stage of the stochastic model. Possible values are ``DSP_FIRST_STAGE`` and ``DSP_SECOND_STAGE``.

   Get number of rows in ``stage``.

.. function:: getNumCols()

   Get number of columns in the extensive form.

.. function:: getNumCols(stage)

   :param stage: stage of the stochastic model. Possible values are ``DSP_FIRST_STAGE`` and ``DSP_SECOND_STAGE``.

   Get number of columns in ``stage``.

.. function:: getObjCoef()

   Get coefficient vector of the objective function.

.. function:: getSolutionStatus()

   Get solution status. Possible values are:

   ============================= ===============================================
   Return Value                  Status
   ============================= ===============================================
   ``DSP_STAT_OPTIMAL``          Optimal
   ``DSP_STAT_PRIM_INFEASIBLE``  Primal infeasible
   ``DSP_STAT_DUAL_INFEASIBLE``  Dual infeasible
   ``DSP_STAT_LIM_ITERorTIME``   Interation limit or solution time limit exceeds
   ``DSP_STAT_STOPPED_GAP``      Tolerance gap reached
   ``DSP_STAT_STOPPED_NODE``     Node limit reached
   ``DSP_STAT_STOPPED_TIME``     Solution time exceeds
   ``DSP_STAT_STOPPED_USER``     Stopped by user
   ``DSP_STAT_STOPPED_SOLUTION`` Stopped with solution issue
   ``DSP_STAT_STOPPED_ITER``     Iteration limit reached
   ``DSP_STAT_STOPPED_UNKNOWN``  Stopped with unknown reason
   ``DSP_STAT_STOPPED_MPI``      MPI error
   ``DSP_STAT_ABORT``            Error
   ``DSP_STAT_LIM_PRIM_OBJ``     Pirmal objective function limit reached
   ``DSP_STAT_LIM_DUAL_OBJ``     Dual objective function limit reached
   ``DSP_STAT_UNKNOWN``          Unknown
   ============================= ===============================================

.. function:: getSolutionTime()

   Get solution time in seconds.

.. function:: getObjValue()

   Get objective function value.

.. function:: getPrimalBound()

   Get the best primal bound of the objective function value.

.. function:: getDualBound()

   Get the best dual bound of the objective function value.

.. function:: getSolution([num::Integer])

   :param num: number of columns to retrieve solution values. Default to ``getNumCols()``

   Get solution for a given size of the solution vector.

.. function:: getNumIterations()

   Get number of iterations.

.. function:: getNumNodes()

   Get number of branch-and-bound nodes. This may be called for solver type ``DSP_SOLVER_BD`` or ``DSP_SOLVER_DE``.

Dual decomposition functions
****************************

.. function:: getDdNumInfeasSolutions()

   Get number of infeasible solutions detected.

.. function:: getDdIterTime()

   Get array of solution time per iteration.

.. function:: getDdMasterTime() 

   Get array of master problem solution time per iteration.

.. function:: getDdSubprobTime()

   Get array of subproblem solution time per iteration.

.. function:: getDdMasterObjValues()

   Get array of master problem objective value per iteration.

.. function:: getDdSubproblemObjValues()

   Get array of subproblem objective value per iteration.

.. function:: getDdPrimalBounds()

   Get array of primal bound per iteration.

.. function:: getDdDualBounds()

   Get array of dual bound per iteration.

.. function:: getDdCpuTime()

   Get solution time in CPU seconds.

.. function:: getDdChangesOfMultiplier()

   Get array of Euclidean distance of dual variable values between two consecutive iterations. This requires to set ``setDdDualVarsLog(DSP_YES)``.

Set Functions
+++++++++++++

Common functions
****************

.. function:: setLogLevel(level::Integer)

   Set display ``level`` from 0 (none) to 5 (verbose).

.. function:: setNumCores(number::Integer)

   Set ``number`` of cores used in OpenMP. This option works for solver type ``DSP_SOLVER_BD`` only.

.. function:: setNodeLimit(number::Integer)

   Set maximum ``number`` of branch-and-bound nodes for solver type ``DSP_SOLVER_BD`` and ``DSP_SOLVER_DE``.

.. function:: setIterLimit(number::Integer)

   Set maximum ``number`` of iterations for solver types ``DSP_SOLVER_DD`` and ``DSP_SOLVER_BD``.

.. function:: setWallLimit(t::Number)

   Set maximum solution time ``t`` in wall clock.

.. function:: setIntRelax(stage)

   Indicate to relax integrality in ``stage``.

Dual decomposition functions
****************************

.. function:: setDdMasterSolver(solverType)

   Set solver type for the master problem. Possible values of ``solverType`` are ``DSP_SOLVER_DD``, ``DSP_SOLVER_BD`` and ``DSP_SOLVER_DE``.

.. function:: setDdMasterNumCutsPerIter(num::Integer)

   Set number of outer-approximation cuts added to the master problem per iteration.

.. function:: setDdAddFeasCuts(freq::Integer)

   Set frequency of adding feasibility cuts for the Lagrangian subproblems. Possible values of ``freq`` are:

   ================= ==================================
   Value             Description
   ================= ==================================
   -1                Disable
    0                Enable only at the first iteration
    :math:`n \geq 1` Enable at every :math:`n` node(s)
   ================= ==================================

.. function:: setDdAddOptCuts(freq::Integer)

   Set frequency of adding optimality cuts for the Lagrangian subproblems. Possible values of ``freq`` are same as for ``setDdAddFeasCuts(freq:Integer)``.

.. function:: setDdEvalUb(freq::Integer)

   Set frequency of evaluating upper bounds for the original problem. Possible values of ``freq`` are same as for ``setDdAddFeasCuts(freq:Integer)``.

.. function:: setDdStoppingTolerance(tol::Number)

   Set stopping tolerance for the relative gap between the best upper bound and the best lower bound. Default value of ``tol`` is ``1.0e-5``.

.. function:: setDdDualVarsLog(yesNo)

   Indicate if the dual variable values are stored for every iteration. Possible values of ``yesNo`` are ``DSP_YES`` and ``DSP_NO``.

.. function:: setDdTrustRegionSize(num::Integer)

   Set the initial trust region size for the master problem.

Benders decomposition functions
*******************************

.. function:: setBdAugScenarios(num::Integer, scenarios::Array{Int,1})

   Set scenario indices that are extended to the first-stage problem. The scenarios specified in ``scenarios`` are considered in the first stage as an extensive form. The number of elements in ``scenarios`` should also be specified in ``num``.

SCIP functions
**************

.. function:: setScipDisplayFreq(freq::Integer)

   Set frequency of displaying node information lines. Default value of ``freq`` is ``100``.

.. function:: setScipLimitsGap(gap::Number)

   Set stopping tolerance for the relative gap between primal and dual bounds. Default value of ``gap`` is ``0``.

.. function:: setScipLimitsTime(time::Number)

   Set maxima time in seconds to run. Default value of ``time`` is ``1e+20``.

General Decomposition Functions
+++++++++++++++++++++++++++++++

.. note:: Calling any function in this section signals to DSP to apply general decomposition to the problem. These functions must be called before calling ``loadProblem``.

.. function:: addCouplingConstraint(m::JuMP.Model, constr::JuMP.LinearConstraint)

   Add a coupling constraint for general decomposition.

.. function:: setVarSubproblem(m::JuMP.Model, var::JuMP.Variable, subproblem::Int)

   Indicate to which subproblem a variable should belong in general decomposition.