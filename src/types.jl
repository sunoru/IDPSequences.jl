using FASTX
import FASTX: identifier, description, sequence
using BioSequences

struct Feature
    id::String
    range::UnitRange{Int64}
end

struct Record
    identifier::String
    description::String
    sequence::LongAminoAcidSeq
    motifs::Vector{Feature}
end
Record(fasta::FASTA.Record, motifs) = let d = description(fasta)
    Record(
        identifier(fasta),
        isnothing(d) ? "" : d,
        sequence(fasta),
        motifs
    )
end
const RecordList = AbstractVector{Record}

identifier(r::Record) = r.identifier
description(r::Record) = r.description
sequence(r::Record) = r.sequence
motifs(r::Record) = r.motifs
Base.length(r::Record) = length(sequence(r))

struct FeaturedSequence
    sequence::LongAminoAcidSeq
    features::Vector{Set{String}}
    lowered::Vector{Bool}
end
const SequenceList = AbstractVector{FeaturedSequence}

Base.length(s::FeaturedSequence) = length(s.sequence)
@inline Base.getindex(s::FeaturedSequence, i) = s.sequence[i]
Base.setindex(s::FeaturedSequence, val, i) = Base.setindex(s.sequence, val, i)
sequence(s::FeaturedSequence) = s.sequence

FeaturedSequence(s, f) = FeaturedSequence(s, f, zeros(Bool, length(s)))

function FeaturedSequence(s::LongAminoAcidSeq, motifs::AbstractVector{Feature} = [])
    n = length(s)
    features = [Set{String}() for _ in 1:n]
    for motif in motifs
        for i in motif.range
            push!(features[i], motif.id)
        end
    end
    FeaturedSequence(s, features, zeros(Bool, n))
end

FeaturedSequence(r::Record) = FeaturedSequence(sequence(r), motifs(r))
