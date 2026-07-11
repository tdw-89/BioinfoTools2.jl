module Paralogs

using DataFrames

function rbh_ds(paralog_df::DataFrame)
    @assert typeof(paralog_df[1,1]) <: AbstractString
    @assert typeof(paralog_df[1,2]) <: AbstractString
    @assert typeof(paralog_df[1,3]) <: AbstractFloat

    unique_ids = unique(vcat(paralog_df[:,1], paralog_df[:,2]))
    ids_to_ind_dict = Dict(unique_ids[i] => i for i in eachindex(unique_ids))
    ind_to_ids_dict = Dict(i => unique_ids[i] for i in eachindex(unique_ids))

    # Create a matrix of zeros
    rbh_matrix = zeros(Float64, length(unique_ids), length(unique_ids))
    orig_mat = zeros(Float64, length(unique_ids), length(unique_ids))

    # Fill in the matrix, treating gene 'i' as the query and gene 'j' as the subject
    for row in eachrow(paralog_df)
        
        i = ids_to_ind_dict[row[1]]
        j = ids_to_ind_dict[row[2]]

        orig_mat[i,j] = row[3]
        orig_mat[j,i] = row[3]
        rbh_matrix[i,j] = row[3]
        rbh_matrix[j,i] = row[3]
    end

    rbh_gene, rbh_paralog = String[], String[]
    matched_inds = Int[]
    score_i_origs, score_j_origs, max_scores, mean_scores = Float64[], Float64[], Float64[], Float64[]
    for i in 1:size(rbh_matrix)[1]

        _, max_i = findmax(rbh_matrix[i,:])[1:2]
        _, max_j = findmax(rbh_matrix[:,max_i])[1:2]
        score_i_orig = orig_mat[i,max_i]
        score_j_orig = orig_mat[max_i,i]
        max_score = max(score_i_orig, score_j_orig)
        mean_score = (score_i_orig + score_j_orig) / 2
        
        if i == max_j && max_i ∉ matched_inds && max_j ∉ matched_inds
            
            id_i = ind_to_ids_dict[max_i]
            id_j = ind_to_ids_dict[max_j]
            
            push!(rbh_gene, id_i)
            push!(rbh_paralog, id_j)
            push!(matched_inds, max_i)
            push!(matched_inds, max_j)
            push!(score_i_origs, score_i_orig)
            push!(score_j_origs, score_j_orig)
            push!(max_scores, max_score)
            push!(mean_scores, mean_score)
        end
    end

    return DataFrame(
        "GeneID" => rbh_gene, 
        "ParalogID" => rbh_paralog, 
        "ds" => score_i_origs
    )
end

"""
    rbh(paralog_df; scoring="max")

Identify reciprocal best hits (RBH) between paralogs based on similarity scores.

# Arguments
- `paralog_df::DataFrame`: DataFrame with at least 4 columns:
  1. GeneID (String)
  2. ParalogID (String)
  3. Percent identity from gene to paralog (Float), or dS value if using `scoring="ds"`
  4. Percent identity from paralog to gene (Float), or ignored if using `scoring="ds"`
- `scoring::String="max"`: Scoring method for determining best hits
  - `"max"` or `"maximum"`: Use maximum of bidirectional scores
  - `"mean"`, `"avg"`, or `"average"`: Use mean of bidirectional scores
  - `"double_max"`: Use original bidirectional scores
  - `"ds"`: Use dS values (column 3 should contain dS values)

# Returns
- `DataFrame`: Contains reciprocal best hit pairs with columns:
  - `GeneID`: First gene in RBH pair
  - `ParalogID`: Second gene in RBH pair
  - `perc_1`: Original percent identity (gene → paralog)
  - `perc_2`: Original percent identity (paralog → gene)
  - `max_perc`: Maximum of the two scores
  - `mean_perc`: Mean of the two scores

# Examples
```julia
df = DataFrame(
    GeneID = ["A", "B", "C"],
    ParalogID = ["B", "A", "D"],
    Perc1 = [95.0, 94.0, 80.0],
    Perc2 = [94.0, 95.0, 85.0]
)
rbh_pairs = rbh(df; scoring="max")
```
"""
function rbh(paralog_df::DataFrame; scoring::String="max")

    scoring = lowercase(scoring)

    if !(scoring in ["max", "maximum", "double_max", "mean", "average", "avg", "ds"])

        error("Invalid scoring method. Must be 'ds', 'max', 'maximum', 'double_max', 'mean', 'avg', or 'average'.")
    end

    scoring = scoring in ["max", "maximum"] ? "max" : scoring in ["mean", "avg", "average"] ? "mean" : scoring

    if scoring == "ds"
        return rbh_ds(paralog_df)
    end

    @assert typeof(paralog_df[1,1]) <: AbstractString
    @assert typeof(paralog_df[1,2]) <: AbstractString
    @assert typeof(paralog_df[1,3]) <: AbstractFloat
    @assert typeof(paralog_df[1,4]) <: AbstractFloat

    unique_ids = unique(vcat(paralog_df[:,1], paralog_df[:,2]))
    ids_to_ind_dict = Dict(unique_ids[i] => i for i in eachindex(unique_ids))
    ind_to_ids_dict = Dict(i => unique_ids[i] for i in eachindex(unique_ids))

    # Create a matrix of zeros
    rbh_matrix = zeros(Float64, length(unique_ids), length(unique_ids))
    orig_mat = zeros(Float64, length(unique_ids), length(unique_ids))

    # Fill in the matrix, treating gene 'i' as the query and gene 'j' as the subject
    for row in eachrow(paralog_df)

        i = ids_to_ind_dict[row[1]]
        j = ids_to_ind_dict[row[2]]
        
        orig_mat[i,j] = row[3]
        orig_mat[j,i] = row[4]

        if scoring == "max"

            max_perc_temp = max(row[3], row[4])
            rbh_matrix[i,j] = max_perc_temp
            rbh_matrix[j,i] = max_perc_temp
        elseif scoring == "mean"

            mean_perc_temp = (row[3] + row[4]) / 2
            rbh_matrix[i,j] = mean_perc_temp
            rbh_matrix[j,i] = mean_perc_temp
        else

            rbh_matrix[i,j] = orig_mat[i,j]
            rbh_matrix[j,i] = orig_mat[j,i]
        end
    end

    rbh_gene, rbh_paralog = String[], String[]
    matched_inds = Int[]
    perc_i_origs, perc_j_origs, max_percs, mean_percs = Float64[], Float64[], Float64[], Float64[]
        
    for i in 1:size(rbh_matrix)[1]

        _, max_i = findmax(rbh_matrix[i,:])[1:2]
        _, max_j = findmax(rbh_matrix[:,max_i])[1:2]
        perc_i_orig = orig_mat[i,max_i]
        perc_j_orig = orig_mat[max_i,i]
        max_perc = max(perc_i_orig, perc_j_orig)
        mean_perc = (perc_i_orig + perc_j_orig) / 2
        
        if i == max_j && max_i ∉ matched_inds && max_j ∉ matched_inds
            
            id_i = ind_to_ids_dict[max_i]
            id_j = ind_to_ids_dict[max_j]
            
            push!(rbh_gene, id_i)
            push!(rbh_paralog, id_j)
            push!(matched_inds, max_i)
            push!(matched_inds, max_j)
            push!(perc_i_origs, perc_i_orig)
            push!(perc_j_origs, perc_j_orig)
            push!(max_percs, max_perc)
            push!(mean_percs, mean_perc)
        end
    end

    return DataFrame("GeneID" => rbh_gene, "ParalogID" => rbh_paralog, "perc_1" => perc_i_origs, "perc_2" => perc_j_origs, "max_perc" => max_percs, "mean_perc" => mean_percs)
end

export 
    rbh, 
    rbh_ds

end