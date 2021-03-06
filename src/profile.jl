import BioAlignments: AbstractSubstitutionMatrix, AbstractScoreModel
using BioSequences

abstract type Profile end
struct ScoreProfile <: Profile
    # AlphabetCount × n
    data::Matrix{Float64}
end

@inline Base.getindex(p::ScoreProfile, aa::AminoAcid, i) = let m = index(aa)
    if m ≤ AlphabetCount
        p.data[m, i]
    elseif aa ≡ AA_B
        (p.data[index(AA_D), i] + p.data[index(AA_N), i]) * 0.5
    elseif aa ≡ AA_Z
        (p.data[index(AA_E), i] + p.data[index(AA_Q), i]) * 0.5
    else
        sum(p.data[:, i] * 0.05)
    end
end
@inline Base.getindex(p::ScoreProfile, m, i) = p.data[m, i]
Base.setindex!(p::ScoreProfile, val, m, i) = p.data[m, i] = val

function create_profile(sequences::SequenceList)
    n = length(sequences[1])
    profile = zeros(AlphabetCount, n)
    for seq in sequences
        @assert length(seq) ≤ n
        @inbounds for i in 1:n
            r = seq[i]
            if r ≡ AA_B
                profile[index(AA_D), i] += 0.5
                profile[index(AA_N), i] += 0.5
            elseif r ≡ AA_Z
                profile[index(AA_E), i] += 0.5
                profile[index(AA_Q), i] += 0.5
            elseif r ∉ Alphabets
                profile[:, i] .+= 0.05
            else
                profile[index(r), i] += 1.0
            end
        end
    end
    ScoreProfile(profile)
end

function ScoreProfile(
    sequences::SequenceList, submat::AbstractSubstitutionMatrix
)
    profile = create_profile(sequences)
    nseq = length(sequences)
    profile.data ./= nseq
    n = size(profile.data, 2)
    for i in 1:n
        profile[:, i] = submat.data[1:AlphabetCount, 1:AlphabetCount] * profile[:, i]
    end
    # TODO: add B / Z / X support
    profile
end


struct FeatureProfile <: Profile
    data::Dict{String, Vector{Float64}}
    motifs::Set{String}
    motif_weight::Float64
end
@inline Base.getindex(p::FeatureProfile, feature) = p.data[feature]
Base.setindex!(p::FeatureProfile, val, feature) = p.data[feature] = val

function FeatureProfile(
    records::RecordList,
    motif_weight
) 
    motifs = Set{String}()
    for record in records
        for motif in record.motifs
            push!(motifs, motif.id)
        end
    end
    FeatureProfile(Dict(), motifs, motif_weight)
end

function get_occurences(fp::FeatureProfile, sequences::SequenceList)
    n = length(sequences[1])
    occurences = Dict(motif => zeros(n) for motif in fp.motifs)
    for i in 1:n
        for seq in sequences
            @assert length(seq) ≤ n
            for motif in seq.features[i]
                occurences[motif][i] += 1.0
            end
        end
    end
    occurences
end

function score_motif(fp::FeatureProfile, occurences, i, motif)
    # Assuming all motifs' probabilities are 1
    occurences[motif][i] * fp.motif_weight
end

function update!(fp::FeatureProfile, sequences::SequenceList) 
    occurences = get_occurences(fp, sequences)
    nseq = length(sequences)
    for (_, occ) in occurences
        occ ./= nseq
    end
    n = length(sequences[1])
    for motif in fp.motifs
        fp.data[motif] = [
            score_motif(fp, occurences, i, motif)
            for i in 1:n
        ]
    end
end


struct DisorderMSAModel{
    S <: AbstractSubstitutionMatrix
} <: AbstractScoreModel{Float64}
    score_profile::Ref{ScoreProfile}
    feature_profile::FeatureProfile
    submat::S
    gap_opening_penalty::Float64
    gap_extension_penalty::Float64
    ending_gaps_penalty::Float64
end
