This contains Farmer example.

* farmer_data.jl defines the model data.
* farmer_model.jl defines the algebraic model using the StochJuMP.jl package.
* farmer_run.jl runs DSP. You can run farmer_run.jl
  * dual decomposition
```bash
    mpirun -n 1 julia farmer_run.jl DD
```
  * Benders decomposition
```bash
    julia farmer_run.jl BD
```
  * Extensive form solution
```bash
    julia farmer_run.jl DE
```

