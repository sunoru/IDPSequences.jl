module IDPSequences

export Feature, Record
export sequence
include("./types.jl")

export DisorderScoringMatrix
include("./data.jl")

export ELM
include("./ELM.jl")

export msa
include("./msa.jl")

end # module
