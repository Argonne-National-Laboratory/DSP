This contains examples for using DSP with Julia interface.

* smps_dd.jl reads SMPS files in subdirectory ./smps and solves the problem using dual decomposition.
```bash
    mpirun -n 1 julia smps_dd.jl
```
* smps_bd.jl reads SMPS files in subdirectory ./smps and solves the problem using Benders decomposition.
```bash
    julia smps_bd.jl
```
