# What is DSP?

DSP (**D**ecomposition of **S**tructured **P**rograms) is an open-source and parallel package that implements decomposition
methods for structured **Mixed-Integer Quadratically Constrained Quadratic Programming (MIQCQP)** problems.
Structured programming problems refer to the class of problems that embed decomposable structures (e.g., block-angular matrices).
A well-known example is a two-stage stochastic program with a finite number of scenarios of the form
$$
\min_{x,y_i} f^x(x) + \sum_{i=1}^N p_i f_i^y(y_i)
$$
subject to
$$
(x,y_i) \in \mathcal{G}_i := \\left\\{(x,y_i):g_i^x(x) + g_i^y(y_i) \leq h_i \\right\\}, \quad \forall i=1,\dots,N,
$$
where $x$ and $y_i$ are the first- and second-stage variables, respectively, for each scenario index $i$, $p_i$ is the probability of scenario $i$, and $f^x$, $f_i^y$, $g_i^x$, and $g_i^y$ are (convex) quadratic functions.
Multiple decomposition methods can effectively utilize such structures in order to accelerate the solutions.

!!! check "Key Features"

    - Implementation of parallel Benders and dual decompositions of stochastic MIQCQP problems
    - Support for Wasserstein-based distributionally robust solutions of the form

    $$
    \min_{(x,y_i)\in\mathcal{G}_i, \; i=1,\dots,N} \left\{ f^x(x) + \max_{p\in\mathbb{P}_r^W(\epsilon)} \sum_{i=1}^N p_i f_i^y(y_i) \right\},
    $$

      where $\mathbb{P}_r^W(\epsilon)$ is the Wasserstein ambiguity set of order $r$ with size $\epsilon$.

    - Support for generic decomposable structures (e.g., network topology, time horizon)
    - Parallel solve on a cluster by using MPI library
    - Run in parallel with Julia modeling interface `DSPopt.jl`. A toy example can be simply written in the following few lines:

    === "Serial Run"

        ```julia
        using StructJuMP
        using DSPopt

        # scenarios
        xi = [[7,7] [11,11] [13,13]]

        # create StructuredModel with number of scenarios
        m = StructuredModel(num_scenarios = 3)

        @variable(m, 0 <= x[i=1:2] <= 5, Int)
        @objective(m, Min, -1.5 * x[1] - 4 * x[2])
        for s = 1:3
            # create a StructuredModel linked to m with id s and probability 1/3
            blk = StructuredModel(parent = m, id = s, prob = 1/3)
            @variable(blk, y[j=1:4], Bin)
            @objective(blk, Min, -16 * y[1] + 19 * y[2] + 23 * y[3] + 28 * y[4])
            @constraint(blk, 2 * y[1] + 3 * y[2] + 4 * y[3] + 5 * y[4] <= xi[1,s] - x[1])
            @constraint(blk, 6 * y[1] + y[2] + 3 * y[3] + 2 * y[4] <= xi[2,s] - x[2])
        end

        # The Wasserstein ambiguity set of order-2 can be imposed with the size limit of 1.0.
        DSPopt.set(WassersteinSet, 2, 1.0)

        status = optimize!(m, 
            is_stochastic = true, # Needs to indicate that the model is a stochastic program.
            solve_type = DSPopt.DW, # see instances(DSPopt.Methods) for other methods
        )
        ```

    === "Paralell Run"

        ```julia
        using MPI
        using StructJuMP
        using DSPopt

        # Initialize MPI communication
        MPI.Init()

        # Initialize DSPopt.jl with the communicator.
        DSPopt.parallelize(MPI.COMM_WORLD)

        # scenarios
        xi = [[7,7] [11,11] [13,13]]

        # create StructuredModel with number of scenarios
        m = StructuredModel(num_scenarios = 3)

        @variable(m, 0 <= x[i=1:2] <= 5, Int)
        @objective(m, Min, -1.5 * x[1] - 4 * x[2])
        for s = 1:3
            # create a StructuredModel linked to m with id s and probability 1/3
            blk = StructuredModel(parent = m, id = s, prob = 1/3)
            @variable(blk, y[j=1:4], Bin)
            @objective(blk, Min, -16 * y[1] + 19 * y[2] + 23 * y[3] + 28 * y[4])
            @constraint(blk, 2 * y[1] + 3 * y[2] + 4 * y[3] + 5 * y[4] <= xi[1,s] - x[1])
            @constraint(blk, 6 * y[1] + y[2] + 3 * y[3] + 2 * y[4] <= xi[2,s] - x[2])
        end

        # The Wasserstein ambiguity set of order-2 can be imposed with the size limit of 1.0.
        DSPopt.set(WassersteinSet, 2, 1.0)

        status = optimize!(m, 
            is_stochastic = true, # Needs to indicate that the model is a stochastic program.
            solve_type = DSPopt.DW, # see instances(DSPopt.Methods) for other methods
        )

        # Finalize MPI communication
        MPI.Finalize()
        ```

## Acknowledgements

This material is based upon work supported by the U.S. Department of Energy, Office of Science, under contract number DE-AC02-06CH11357.