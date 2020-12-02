import BioAlignments: SubstitutionMatrix, parse_ncbi_submat
using BioSequences

const Alphabets = parse.([AminoAcid], collect("ARNDCQEGHILKMFPSTWYV"))
const AlphabetCount = length(Alphabets)

# Scoring Matrix
const DisorderScoringMatrix = parse_ncbi_submat(
    AminoAcid, joinpath(dirname(@__FILE__), "data", "submat", "DISORDER")
)

index(s::AminoAcid) = reinterpret(UInt8, s) + 0x01

const ∞ = 2007012811
