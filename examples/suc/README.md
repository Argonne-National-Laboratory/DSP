This contains an example for two-stage stochastic unit commitment using IEEE 118-bus system.

* suc_mod.jl defines the algebraic model using the StochJuMP.jl package.
* suc_run.jl solves the model by dual decomposition
```bash
    mpirun -n 1 julia suc_run.jl
```
* suc_run_de.jl solves the model by extensive form solution
```bash
    julia suc_run_de.jl
```

