using BioSequences

update_score(
    (score, step), (new_score, new_step)
) = if new_score > score
    new_score, new_step
else
    score, step
end
max_score(scorings...) = reduce(update_score, scorings)

function find_best_score(V, ending_penalty)
    n, m = size(V)
    max_coords = n, m
    maxv = V[n, m]
    for i = 1:n
        v = V[i, m] + ending_penalty * (n - i)
        if v > maxv
            maxv = v
            max_coords = i, m
        end
    end
    for j = 1:m
        v = V[n, j] + ending_penalty * (m - j)
        if v > maxv
            maxv = v
            max_coords = n, j
        end
    end
    max_coords
end

function backtrace_alignment(
    seq, model,
    V, H, G,
    traces_V, traces_H, traces_G
)
    n, m = size(V)
    s1 = AminoAcid[]
    s2 = AminoAcid[]
    features = Set{String}[]
    indices = Int[]
    best_i, best_j = find_best_score(V, model.ending_gaps_penalty)
    i, j = n, m
    if best_i ≠ n || best_j ≠ m
        i, j = best_i, best_j
        for k in n - 1:-1:i
            push!(s1, AA_A)  # Dummy sequence
            push!(s2, AA_Gap)
            push!(features, Set{String}())
            push!(indices, -length(s1))
        end
        for k in m - 1:-1:j
            push!(s1, AA_Gap)
            push!(s2, seq[k])
            push!(features, seq.features[k])
            push!(indices, seq.indices[k])
        end
    end
    current_matrix = :V
    aa1 = AA_A
    aa2 = AA_A
    feature = Set{String}()
    index = -1
    while i > 1 || j > 1
        if i > 1 && j > 1 && current_matrix ≡ :V
            aa1, aa2 = AA_A, seq[j - 1]
            feature = seq.features[j - 1]
            index = seq.indices[j - 1]
            t = traces_V[i - 1, j - 1]
            if t == 2
                current_matrix = :G
            elseif t == 3
                current_matrix = :H
            end
            i -= 1
            j -= 1
        elseif i > 1 && current_matrix ≡ :H
            aa1, aa2 = AA_A, AA_Gap
            feature = Set{String}()
            index = -length(s1) - 1
            t = j > 1 ? traces_H[i - 1, j - 1] : 0
            if t == 1
                current_matrix = :V
            elseif t == 2
                current_matrix = :G
            end
            i -= 1
        elseif j > 1 && current_matrix ≡ :G
            aa1, aa2 = AA_Gap, seq[j - 1]
            feature = seq.features[j - 1]
            index = seq.indices[j - 1]
            t = i > 1 ? traces_G[i - 1, j - 1] : 0
            if t == 1
                current_matrix = :V
            elseif t == 3
                current_matrix = :H
            end
            j -= 1
        end
        push!(s1, aa1)
        push!(s2, aa2)
        push!(features, feature)
        push!(indices, index)
    end
    aligned1 = FeaturedSequence(LongAminoAcidSeq(reverse(s1)))
    aligned2 = FeaturedSequence(LongAminoAcidSeq(reverse(s2)), reverse(features), reverse(indices), seq.length)
    aligned1, aligned2
end

function remove_gaps(s1::FeaturedSequence, s2::FeaturedSequence)
    aligned = AminoAcid[]
    features = Set{String}[]
    indices = Int[]
    n = length(s1)
    for i in 1:n
        if s1[i] ≠ AA_Gap
            push!(aligned, s2[i])
            push!(features, s2.features[i])
            push!(indices, s2.indices[i])
        end
    end
    s = LongAminoAcidSeq(aligned)
    aligned1 = FeaturedSequence(s, features)
    aligned2 = FeaturedSequence(copy(s), features, indices, s2.length)
    aligned1, aligned2
end

function pairwise_align(seq::FeaturedSequence, model::DisorderMSAModel, first_gapped::Bool)
    sp = model.score_profile[]
    fp = model.feature_profile
    n1 = size(sp.data, 2)
    n2 = length(seq)
    # Three matrices: V, H, G
    V = zeros(n1 + 1, n2 + 1)
    H = zeros(n1 + 1, n2 + 1)
    G = zeros(n1 + 1, n2 + 1)
    # Traces
    traces_V = zeros(Int8, n1, n2)
    traces_H = zeros(Int8, n1, n2)
    traces_G = zeros(Int8, n1, n2)

    ending_penalty = model.ending_gaps_penalty
    gap_opening = model.gap_opening_penalty
    gap_extension = model.gap_extension_penalty
    for i = 2:n1 + 1
        V[i, 1] = G[i, 1] = -∞
        H[i, 1] = (i - 1) * ending_penalty
    end
    for i = 2:n2 + 1
        V[1, i] = H[1, i] = -∞
        G[1, i] = (i - 1) * ending_penalty
    end
    for i = 2:n1 + 1
        for j = 2:n2 + 1
            @inbounds aa = seq[j - 1]
            profile_score = sp[aa, i - 1]
            feature_score = 0.0
            for motif in seq.features[j - 1]
                feature_score += fp[motif][i - 1]
            end
            V[i, j], traces_V[i - 1, j - 1] = max_score(
                (V[i - 1, j - 1], 1),
                (G[i - 1, j - 1], 2),
                (H[i - 1, j - 1], 3)
            )
            V[i, j] += profile_score + feature_score
            H[i, j], traces_H[i - 1, j - 1] = max_score(
                (V[i - 1, j] + gap_opening, 1),
                (G[i - 1, j] + gap_opening, 2),
                (H[i - 1, j] + gap_extension, 3)
            )
            G[i, j], traces_G[i - 1, j - 1] = max_score(
                (V[i, j - 1] + gap_opening, 1),
                (G[i, j - 1] + gap_extension, 2),
                (H[i, j - 1] + gap_opening, 3)
            )
        end
    end
    s1, s2 = backtrace_alignment(
        seq, model,
        V, H, G,
        traces_V, traces_H, traces_G
    )
    if !first_gapped
        s1, s2 = remove_gaps(s1, s2)
    end
    s1, s2
end


