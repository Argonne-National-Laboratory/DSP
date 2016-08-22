This contains Farmer example.

* farmer_model.jl defines the algebraic model using the JuMP.jl and Dsp.jl packages.
* farmer_run.jl runs DSP. You can run farmer_run.jl
```bash
    julia farmer_run.jl
```
* farmer_run_mpi.jl runs parallel DSP using MPI library. You can run farmer_run_mpi.jl
```bash
    mpirun -n 2 julia farmer_run_mpi.jl
```
