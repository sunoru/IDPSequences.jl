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

function _get_unipro_id(id)
    for each in split(id, "|")
        isnothing(match(r"^[0-9A-Z]{6}$", each)) || return each
    end
    id
end

identifier(r::Record) = r.identifier
description(r::Record) = r.description
sequence(r::Record) = r.sequence
motifs(r::Record) = r.motifs
unipro_id(r::Record) = _get_unipro_id(r.identifier)
Base.length(r::Record) = length(sequence(r))

struct FeaturedSequence
    sequence::LongAminoAcidSeq
    features::Vector{Set{String}}
    indices::Vector{Int}
    length::Int
end
const SequenceList = AbstractVector{FeaturedSequence}

Base.length(s::FeaturedSequence) = length(s.sequence)
@inline Base.getindex(s::FeaturedSequence, i) = s.sequence[i]
Base.setindex(s::FeaturedSequence, val, i) = Base.setindex(s.sequence, val, i)
sequence(s::FeaturedSequence) = s.sequence

FeaturedSequence(s, f) = FeaturedSequence(s, f, collect(1:length(s)), length(s))

function FeaturedSequence(s::LongAminoAcidSeq, motifs::AbstractVector{Feature} = [])
    n = length(s)
    features = [Set{String}() for _ in 1:n]
    for motif in motifs
        for i in motif.range
            push!(features[i], motif.id)
        end
    end
    FeaturedSequence(s, features, collect(1:n), n)
end

FeaturedSequence(r::Record) = FeaturedSequence(sequence(r), motifs(r))
