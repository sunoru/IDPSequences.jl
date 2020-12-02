module IDPSequences

export Feature, Record, FeaturedSequence
export sequence, identifier, description, motifs
include("./types.jl")

export DisorderScoringMatrix
include("./data.jl")

export ELM
include("./ELM.jl")

export Profile, ScoreProfile, FeatureProfile, DisorderMSAModel
include("./profile.jl")

include("./pairwise.jl")

export load_data
export load_fasta
include("./io.jl")

export msa
include("./msa.jl")

end # module
