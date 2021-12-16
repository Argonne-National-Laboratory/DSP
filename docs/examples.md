# Examples

Here we give some examples for running `DSP` with the executable file.
The examples are given with `--algo dd` (i.e., using the dual decomposition), but can be easily modified to use the other algorithms (e.g., `bd`, `dw`) with argument `--algo`.

!!! Example "Example: Wasserstein DRO"

    This example solves the distributionall robust optimization (DRO) of `farmer` example, where the Wasserstein ambiguity set is of order `2` with the size `0.1`.

    === "Serial Run"

        ```
        ./bin/runDsp --algo drdd \
                     --smps ../examples/farmer \
                     --wassnorm 2.0 \
                     --wasseps 0.1 \
                     --soln mysoln \
                     --param myparam.txt
        ```

    === "Parallel Run"

        ```
        mpiexec -np 3 ./bin/runDsp --algo drdd \
                                   --smps ../examples/farmer \
                                   --wassnorm 2.0 \
                                   --wasseps 0.1 \
                                   --soln mysoln \
                                   --param myparam.txt
        ```

!!! Example "Example: Stochastic MIQCP"

    This examples solves the two-stage stochastic mixed-integer quadratically constrained program (MIQCP) of `farmer` example.
    The quadratic constraints are defined in `../example/farmer.txt`.

    === "Serial Run"

        ```
        ./bin/runDsp --algo dd \
                     --smps ../examples/farmer \
                     --quad ../example/farmer \
                     --soln mysoln \
                     --param myparam.txt
        ```

    === "Parallel Run"

        ```
        mpiexec -np 3 ./bin/runDsp --algo dd \
                                   --smps ../examples/farmer \
                                   --quad ../example/farmer \
                                   --soln mysoln \
                                   --param myparam.txt
        ```

    !!! info
        As a special case, the stochastic mixed-integer linear program can be set by omitting argument `--quad`.

!!! Example "Example: Generic structured program"

    This example solves a generic block-angular structure problem written in files `noswot.mps` and `noswot.dec`.

    === "Serial Run"

        ```
        ./bin/runDsp --algo dw \
                     --mps ../examples/mps-dec/noswot.mps \
                     --dec ../examples/mps-dec/noswot.dec \
                     --soln mysoln \
                     --param myparam.txt
        ```

    === "Parallel Run"

        ```
        mpiexec -np 3 ./bin/runDsp --algo dw \
                                   --mps ../examples/mps-dec/noswot.mps \
                                   --dec ../examples/mps-dec/noswot.dec \
                                   --soln mysoln \
                                   --param myparam.txt
        ```