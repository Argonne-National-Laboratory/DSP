# Numerical experiments for DRO

This file describes how to reproduce the numerical results reported in the paper [Dual Decomposition of Two-Stage Distributionally Robust Mixed-Integer Programming under the Wasserstein Ambiguity Set](https://optimization-online.org/2020/04/7723/).

## Problem instances

The 80 problem instances used in the experiments are available in `examples/dro/experiments` of the DSP repository.
You can also generate the same (or other) instances by using `generate_instances.jl`.

## Solve an instance

You can solve instance `drdcap_233_10_20` with the Wasserstein limit of 1.0 and outputs the results into `outputs/drdcap_233_20_1.txt` by the command below.

```bash
../../../build/bin/runDsp --algo drdd --wassnorm 2.0 --wasseps 1.0 --smps drdcap_233_10_20 --param params.txt > outputs/drdcap_233_20_1.txt
```

The bash script used on a cluster machine is available in `run_script.sh`.

## Create figures and tables

Once all the output files are generated in `outputs`, we can generate figures and tables used in the paper.

1. Run `output_file_to_csv.jl`, reading the output files and collecting them into `tmp.csv` file.
2. Run `csv_to_figs.jl`, reading `tmp.csv` and generating figures and tables.
