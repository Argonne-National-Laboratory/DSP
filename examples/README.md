This contains examples for using DSP with Julia interface.

* smps.jl reads SMPS files given as an argument and solves the problem using serial DSP.
```bash
    julia smps.jl ./smps/sslp_15_45_5
```
* smps.jl reads SMPS files given as an argument and solves the problem using parallel DSP.
```bash
    mpirun -n 2 julia smps_mpi.jl ./smps/sslp_15_45_5
```
