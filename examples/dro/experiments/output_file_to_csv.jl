"""
A Julia script to read output files and create a latex table format.
"""

using Printf

instance_numbers = [233, 243, 332, 342]
scenarios = [20, 50, 100, 200, 300]
epsilons = [1, 100, 500, 1000]

mycsv = open("tmp.csv", "w")
println(mycsv, "instance,K,eps,Iterations,UB,LB,Time")

for i in instance_numbers, s in scenarios
    for e in epsilons
        local filename = "outputs/drdcap_$(i)_$(s)_$(e).txt"
        primal = dual = gap = time = 0.0
        iter = 0
        # print("dcap$(i)\\_$(s)\\_$(e) & ")
        open(filename, "r") do io
            while !eof(io)
                word = readline(io)
                if length(word) < 5
                    continue
                end
                if word[begin:5] == "Statu"
                    if word[9:end] != "3000"
                        @error "Unexpected solution status"
                    end
                elseif word[begin:5] == "Prima"
                    primal = parse(Float64, word[15:end])
                    # @show primal
                elseif word[begin:5] == "Dual "
                    dual = parse(Float64, word[15:end])
                    # @show dual
                elseif word[begin:5] == "Gap ("
                    gap = parse(Float64, word[15:end])
                    # @show gap
                elseif word[begin:5] == "Itera"
                    iter = parse(Float64, word[15:end])
                    # @show iter
                elseif word[begin:5] == "Time "
                    time = parse(Float64, word[15:end])
                    # @show time
                end
            end
        end
        @printf(mycsv, "%d,%d,%d,%d,%f,%f,%f\n", i, s, e, iter, primal, dual, time)
    end
    # println("\\hline")
end

close(mycsv)
