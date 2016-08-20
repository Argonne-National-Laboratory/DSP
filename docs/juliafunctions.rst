Julia Interface Functions
-------------------------

We describe Julia interface functions.

.. function:: readSmps(filename)

   Read SMPS files from files ``filename``.cor, ``filename``.tim and ``filename``.sto.

.. function:: getblocksolution(model::JuMP.Model)

   Return second-stage solutions from stochastic model ``model::JuMP.Model``.

.. function:: JuMP.solve(model::JuMP.Model; suppress_warnings = false, comm = nothing, options...)

   Solve the model and returns solution status. Argument ``comm`` may have a communicator from MPI.jl. This can also have arguments ``solve_type`` and ``param``. Possible values of ``solve_type`` are ``:Dual``, ``:Benders``, and ``:Extensive``. ``param`` may have a parameter file name.

.. function:: optimize(;suppress_warnings = false, comm = nothing, options...)

   Same as ``JuMP.solve``, but this is used when model is read from SMPS files.

.. function:: getprimobjval()

   Returns primal objective value. The value can also be available from ``Dsp.model.primVal``.

.. function:: getdualobjval()

   Returns dual objective value. The value can also be available from ``Dsp.model.dualVal``.

.. function:: getsolutiontime()

   Returns solution time in seconds.

.. function:: getblockids()

   Returns block IDs required for the current processor when MPI.jl is initialized.