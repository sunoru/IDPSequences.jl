module ELM

using HTTP, YAML
using FASTX

using ..IDPSequences: Feature, Record, RecordList

function get_data(unipro_id; read_fasta = true)
    r = HTTP.request("GET", "http://elm.eu.org/instances.gff?q=$unipro_id")
    @assert r.status == 200
    motifs = Feature[]
    fasta_lines = []
    reading_fasta = false
    for line in split(String(r.body), "\n")
        if reading_fasta
            if startswith(line, ">")
                continue
            end
            push!(fasta_lines, line)
        elseif startswith(line, "##FASTA")
            read_fasta || break
            reading_fasta = true
        elseif contains(line, "sequence_feature")
            tmp = split(line)
            start, ending = parse.([Int], tmp[4:5])
            id = split(tmp[9], "=")[2]
            push!(motifs, Feature(id, start:ending))
        end
    end
    fasta = read_fasta ? FASTA.Record(unipro_id, join(fasta_lines)) : FASTA.Record
    Record(fasta, motifs)
end

function save_yaml(filename, records::RecordList)
    data = Dict(
        record.identifier => map(record.motifs) do motif
            t = Dict()
            t["id"] = motif.id
            t["range"] = [first(motif.range), last(motif.range)]
            t
        end for record in records
    )
    YAML.write_file(filename, data)
end

function load_yaml(filename)
    data = YAML.load_file(filename)
    Dict(
        id => isnothing(motifs) ? [] : [
            Feature(t["id"], t["range"][1]:t["range"][2])
            for t in motifs
        ] for (id, motifs) in data
    )
end

end
