# How to Use DSP

There are two ways to use algorithms in `DSP`:

- Executable binary file `runDsp`
- Callable C functions in shared library `libDsp.so` (or `libDsp.dylib` for Mac)

Building the DSP source files creates an executable file `./bin/runDsp` and a shared library file `./lib/libDsp.*`.

## Executable file

The executable file `./bin/runDsp` can be run as a command-line tool with the following arguments:

| Argument     | Description |
| ------------ | ----------- |
| `--algo`     | A required argument for the choice of algorithms from one of the following:<ul><li>`de`: deterministic equivalent form</li><li>`bd`: (integer) Benders decomposition</li><li>`dd`: dual decomposition</li><li>`drbd`: distributionally robust (integer) Bender decomposition</li><li>`drdd`: distributionally robust dual decomposition</li><li>`dw`: Dantzig-Wolfe decomposition with branch-and-bound</li></ul> |
| `--wassnorm` | Wasserstein distance norm (>= 1.0). This argument should be used with `--wasseps` argument when `--algo drbd` or `--algo drdd`.|
| `--wasseps`  | Wasserstein distance limit (>= 0.0). This argument should be used with `--wassnorm` argument when `--algo drbd` or `--algo drdd`. |
| `--smps`     | [SMPS](https://ieeexplore.ieee.org/abstract/document/8142546) file name without extensions.<br>For example, if your SMPS files are `../test/farmer.cor`, `../test/farmer.sto`, and `../test/farmer.tim`, this value should be `../test/farmer` |
| `--mps`      | [MPS](https://neos-guide.org/content/mps-format) file name. This argument should be used with `--dec` argument. |
| `--dec`      | [DEC](https://gcg.or.rwth-aachen.de/doc/reader__dec_8h.html) file name. This argument should be used with `--mps` argument. |
| `--quad`     | Quadratic file name without extension. The Quadratic file should extend the parsed SMPS file and its suffix needs to be `.txt` in the current version. For example, if your Quadratic file is `../test/farmer.txt`, this value should be `../test/farmer`. |
| `--soln`     | An optional argument for solution file prefix.<br>For example, For example, if `--soln mysoln` is given, `runDsp` will write the following solution files:<ul><li>`mysoln.primobj.txt` for primal objective value</li><li>`mysoln.dualobj.txt` for dual objective value</li><li>`mysoln.primal.txt` for primal variable values in the order of variables defined in the input file, if the primal objective value is less than `1e+20`.</li><li>`mysoln.dual.txt` for dual variable values if `--algo dd` is given.</li></ul> |
| `--param`    | An optional paramater for parameter file name (see [Parameters](./parameters.md) section)|


!!! Example "Your first example!"

    This examples solves the two-stage stochastic mixed-integer linear program (MILP) of `farmer` example by using the (parallel) dual decomposition.

    === "Serial Run"

        ```
        ./bin/runDsp --algo dd --smps ../examples/smps/farmer
        ```

    === "Parallel Run"

        ```
        mpiexec -np 3 ./bin/runDsp --algo dd --smps ../examples/smps/farmer
        ```

    !!! attention

        Of course, `Parallel Run` is available only if you build `DSP` with MPI library.

## Shared library

The shared library provides access to C API functions.
The library needs to be placed in the searchable path.
For example,

=== "Linux"

    ```
    export LD_LIBRARY_PATH=/path/to/lib/libDsp.so:$LD_LIBRARY_PATH
    ```

=== "Mac"

    ```
    export DYLD_LIBRARY_PATH=/path/to/lib/libDsp.dylib:$DYLD_LIBRARY_PATH
    ```

Using the library can be useful to interface with other (high-level) langauges such as Julia.
We strongly recomment to check our Julia interface [DSPopt.jl](https://github.com/kibaekkim/DSPopt.jl).