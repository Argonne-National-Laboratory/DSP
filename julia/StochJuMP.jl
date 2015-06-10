#module StochJuMP

if Pkg.installed("JuMP") == nothing
	Pkg.add("JuMP")
end

# import modules
import JuMP
import MathProgBase
import MathProgBase.MathProgSolverInterface

export StochasticModel, getStochastic, getparent, getchildren, 
       num_scenarios, StochasticBlock, @second_stage

#############################################################################
# JuMP exports
export
# Objects
    Model, Variable, AffExpr, QuadExpr, LinearConstraint, QuadConstraint, MultivarDict,
    ConstraintRef,
# Functions
    # Model related
    getNumVars, getNumConstraints, getObjectiveValue, getObjective,
    getObjectiveSense, setObjectiveSense, writeLP, writeMPS, setObjective,
    addConstraint, addSOS1, addSOS2, solve,
    getInternalModel, setPresolve, buildInternalModel,
    # Variable
    setName, getName, setLower, setUpper, getLower, getUpper,
    getValue, setValue, getDual,
    # Expressions and constraints
    affToStr, quadToStr, conToStr, chgConstrRHS,
# Macros and support functions
    @addConstraint, @addConstraints, @defVar, 
    @defConstrRef, @setObjective, addToExpression,
    @setNLObjective, @addNLConstraint, @gendict


# Define stochastic model
type StochasticData
	probability::Vector{Float64}
    children::Vector{JuMP.Model}
    parent
    num_scen::Int
end

# constructor with no argument
StochasticData() = StochasticData(Float64[],JuMP.Model[],nothing,0)
StochasticData(children, parent, nscen) = StochasticData(Float64[],children, parent, nscen)

# constructor with the number of scenarios
function StochasticModel(numScen::Int)
    m = JuMP.Model()
    m.ext[:Stochastic] = StochasticData(JuMP.Model[],nothing,numScen)
    return m
end

# constructor with child and parent models
StochasticModel(children,parent) = StochasticModel(children,parent,0)
function StochasticModel(children, parent, nscen)
    m = JuMP.Model(solver=parent.solver)
    m.ext[:Stochastic] = StochasticData(children, parent, nscen)
    m.ext[:Skip] = false;
    return m
end

function getStochastic(m::JuMP.Model)
    if haskey(m.ext, :Stochastic)
        return m.ext[:Stochastic]
    else
        error("This functionality is only available for StochasticModels")
    end
end

getparent(m::JuMP.Model)      = getStochastic(m).parent
getchildren(m::JuMP.Model)    = getStochastic(m).children
getprobability(m::JuMP.Model) = getStochastic(m).probability
num_scenarios(m::JuMP.Model)  = getStochastic(m).num_scen

StochasticBlock(m::JuMP.Model) = StochasticBlock(m::JuMP.Model, 1.0/num_scenarios(m))
function StochasticBlock(m::JuMP.Model, probability::Float64)
    stoch = getStochastic(m)
    #stoch.parent = m
    ch = StochasticModel(JuMP.Model[], m)
    push!(stoch.children, ch)
    push!(stoch.probability, probability)
    return ch
end

# create sparse matrix in [T W] form, where T is technology matrix and W is recourse matrix.
# taken and modified from benders.jl
function prepConstrMatrix(m::JuMP.Model)
    stoch = getStochastic(m)
    if stoch.parent == nothing
    	return JuMP.prepConstrMatrix(m)
    else
        rind = Int[]
        cind = Int[]
        value = Float64[]
        linconstr = deepcopy(m.linconstr)
        for (nrow,con) in enumerate(linconstr)
            aff = con.terms
            for (var,id) in zip(reverse(aff.vars), length(aff.vars):-1:1)
                push!(rind, nrow)
                if m.linconstr[nrow].terms.vars[id].m == stoch.parent
                    push!(cind, var.col)
                elseif m.linconstr[nrow].terms.vars[id].m == m
                    push!(cind, stoch.parent.numCols + var.col)
                end
                push!(value, aff.coeffs[id])
                splice!(aff.vars, id)
                splice!(aff.coeffs, id)
            end
        end
    end
    return sparse(rind, cind, value, length(m.linconstr), stoch.parent.numCols + m.numCols)
end


#end
