# Parameters

`DSP` provides a number of parameters that changes the algorithmic settings.

## Parameter file format

Users can write and pass parameter file to `./runDsp` with `--param` option.

In parameter file, a parameter value can be set in the following format:

```
<type> <name> <value>
```

For example, a parameter file containing

```
int LOG_LEVEL 2
double DW/TIME_LIM 300
#bool DW/HEURISTICS false
bool DW/HEURISTICS/SMIP false
```

changes the value of parameter `LOG_LEVEL` of type `int` to `2`.
Similarly for the other parameters.

> NOTE: With `#` the third line will be commented out.

## List of Parameters

This will be available soon.
Meanwhile, please refer `./src/Utility/DspParams.cpp` for the parameters available.
