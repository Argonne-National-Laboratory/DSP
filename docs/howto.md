# How to Use DSP

Building the DSP source files creates an executable file `./bin/runDsp` and a shared library file `./lib/libDsp.so` (or `./lib/libDsp.dylib` for Mac).

## Executable file

The executable file can be run as a command-line tool.
Running the file without any argument will display available arguments for the running.

```shell
./runDsp
```

```shell
Usage: --algo <de,bd,dd,drbd,drdd,dw> --smps <smps file> --mps <mps file> --dec <dec file> [--soln <solution file prefix> --param <param file> --test <benchmark objective value>]

       --algo	choice of algorithms.
                de: deterministic equivalent form
                bd: Benders decomposition
                dd: dual decomposition
                drbd: distributionally robust bd
                drdd: distributionally robust dd
                dw: Dantzig-Wolfe decomposition with branch-and-bound
       --smps	SMPS file name without extensions. For example, if your SMPS files are ../test/farmer.cor, ../test/farmer.sto, and ../test/farmer.tim, this value should be ../test/farmer
       --mps	MPS file name
       --dec	DEC file name
       --soln	optional argument for solution file prefix. For example, if the prefix is given as MySol, then two files MySol.primal.txt and MySol.dual.txt will be written for primal and dual solutions, respectively.
       --param	optional paramater for parameter file name
```

### Input files

`runDsp` can read problem from

- [SMPS](https://ieeexplore.ieee.org/abstract/document/8142546) files: This input should consist of three files with extensions: `.cor`, `.tim`, and `sto`.
- `.dro` file: For `--algo dro`, `runDsp` requires `.dro` file, in addition to SMPS files, that defines the Wasserstein uncertainty set.
- [MPS and DEC files](https://gcg.or.rwth-aachen.de/doc/reader__dec_8h.html): This input should consist of two files, which can be given to options `--mps` and `--dec`.

DSP has a number of parameters for algorithms.
Users can set parameter values by giving parameter file with `--param`.
See more details in [Parameters](parameters.md).
This argument is optional.


### Output files

If `--soln` is givne, `runDsp` will write solution files.
For example, if `--soln mysoln` is given, `runDsp` will write

- `mysoln.primobj.txt` for primal objective value
- `mysoln.dualobj.txt` for dual objective value
- `mysoln.primal.txt` for primal variable values in the order of variables defined in the input file, if the primal objective value is less than `1e+20`.
- `mysoln.dual.txt` for dual variable values if `--algo dd` is given.

## Shared library

The shared library provides access to C API functions.
The library needs to be placed in the searchable path (e.g., `LD_LIBRARY_PATH` for linux or `DYLD_LIBRARY_PATH` for Mac).

Using the library can be useful to interface with other (high-level) langauges such as Julia.
We strongly recomment to check our Julia interface [DSPopt.jl](https://github.com/kibaekkim/DSPopt.jl).