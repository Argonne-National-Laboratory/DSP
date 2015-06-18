#### Example: Stochastic Unit Commitment with IEEE 118-Bus System

This contains an example for stochastic unit commitment with IEEE 118-bus system. The model considers uncertain wind power generation and quick-start generating units (introducing second-stage integer variables).

* suc_mod.jl defines the algebraic model using the StochJuMP.jl package.
* suc_run.jl solves the model by dual decomposition
```bash
    mpirun -n 1 julia suc_run.jl
```
* suc_run_de.jl solves the model by extensive form solution
```bash
    julia suc_run_de.jl
```

