This contains Farmer example modeled in extensive form for illustration purposes.

* ext_farmer_model.jl defines the algebraic model using the StochJuMP.jl package.
* ext_farmer_run.jl runs DSP. You can run ext_farmer_run.jl
  * dual decomposition
```bash
    mpirun -n 1 julia ext_farmer_run.jl DD
```
  * Benders decomposition
```bash
    julia ext_farmer_run.jl BD
```
  * Extensive form solution
```bash
    julia ext_farmer_run.jl DE
```

