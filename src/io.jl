using Printf
using FASTX

load_fasta(fasta_file) = collect(open(FASTA.Reader, fasta_file))

function load_data(
    fasta_file,
    motif_file = fasta_file[1:first(findlast(".", fasta_file))] * "yaml"
)
    motifs = ELM.load_yaml(motif_file)
    records = Record[]
    open(fasta_file) do fi
        for fasta in load_fasta(fasta_file)
            push!(records, Record(
                fasta, get(motifs, _get_unipro_id(identifier(fasta)), [])
            ))
        end
    end
    records
end

function Base.convert(::Type{S}, seq::FeaturedSequence) where {S <: AbstractString}
    chars = Char[]
    last_index = 0
    for (i, x) in enumerate(seq.sequence)
        c = Char(x)
        if x ≠ AA_Gap
            if seq.indices[i] ≠ last_index + 1
                if length(chars) > 0
                    chars[end] = lowercase(chars[end])
                end
                c = lowercase(c)
            end
            last_index = seq.indices[i]
        end
        push!(chars, c)
    end
    if seq.indices[end] < seq.length
        chars[end] = lowercase(chars[end])
    end
    S(chars)
end
Base.String(seq::FeaturedSequence) = convert(String, seq)

Base.print(io::IO, seq::FeaturedSequence) = print(io, String(seq))
Base.print(io::IO, seqs::SequenceList) = print(io, join(
    (String(seq) for seq in seqs),
    '\n'
))

function get_next_index(seq, start)
    for i in start + 1:length(seq)
        seq.indices[i] > 0 && return i
    end
    -∞
end

function save_msf(filename::AbstractString, sequences::SequenceList, records::RecordList)
    nseqs = length(sequences)
    @assert nseqs > 0 && nseqs ≡ length(records)
    originals = sequence.(records)
    original_lengths = length.(originals)
    # TODO: regenerate the gaps
    # aa_seqs = [AminoAcid[] for _ in 1:nseqs]
    aa_seqs = collect.(String.(sequences))

    open(filename, "w") do fo
        n = length(aa_seqs[1])
        # TODO: checksum?
        @printf fo "PileUp\n\n\n\n   MSF:%5d  Type: P    Check:  0000   ..\n\n" n
        for record in records
            @printf fo "Name: %s oo  Len:%5d  Check:  0000  Weight:  10.0\n" unipro_id(record) n
        end
        @printf fo "\n//\n\n"
        line_startings = [["$(unipro_id(record))     "] for record in records]
        i = 0
        while i < n
            @printf fo "\n\n"
            lines = deepcopy(line_startings)
            for _ in 1:5
                for (j, seq) in enumerate(aa_seqs)
                    push!(lines[j], join(c ≡ '-' ? '.' : c for c in seq[i+1:min(i+10, n)]))
                end
                i += 10
                i ≥ n && break
            end
            for line in lines
                println(fo, join(line, " "))
            end
        end
    end
end
