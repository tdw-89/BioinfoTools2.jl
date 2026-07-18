module SOTerms
#========= Build the SO term reference =========#
using JSON

const ASSETS = joinpath(pkgdir(SOTerms), "assets")
const SO_JSON = joinpath(ASSETS, "SOFA.json")
include_dependency(SO_JSON)

"""Stack-allocated, flexible lookup for SO terms (SOFA subset) for GFF3 files (column 3)"""
struct SOTermLookup{T1<:Tuple,T2<:Tuple,T3<:Tuple}
    by_short::T1
    by_full::T2
    by_label::T3
end

function Base.getindex(lookup::SOTermLookup, short::UInt16)
    if short >= 1 && short <= length(lookup.by_short)
        return lookup.by_short[short]
    else
        return nothing
    end
end

function Base.getindex(lookup::SOTermLookup, full::Int)
    data = lookup.by_full
    left, right = 1, length(data)
    while left <= right
        mid = (left + right) >>> 1
        val = data[mid][1]
        if val == full
            return (data[mid][2], data[mid][3])
        elseif val < full
            left = mid + 1
        else
            right = mid - 1
        end
    end
    return nothing
end

function Base.getindex(lookup::SOTermLookup, label::Symbol)
    data = lookup.by_label
    left, right = 1, length(data)
    while left <= right
        mid = (left + right) >>> 1
        val = data[mid][1]
        if val == label
            return (data[mid][2], data[mid][3])
        elseif val < label
            left = mid + 1
        else
            right = mid - 1
        end
    end
    return nothing
end

Base.getindex(lookup::SOTermLookup, label::AbstractString) = lookup[Symbol(label)]

Base.show(io::IO, lookup::SOTermLookup) =
    print(io, "SOTermLookup(", length(lookup.by_short), " terms)")

function _build_lookups()
    raw_json = JSON.parse(read(SO_JSON, String))

    arr_by_short = Tuple{Int,Symbol}[]
    arr_by_full = Tuple{Int,UInt16,Symbol}[]
    arr_by_label = Tuple{Symbol,UInt16,Int}[]

    def_dict = Dict{Int,String}()

    for node in raw_json["graphs"][1]["nodes"]
        if haskey(node, "id")
            id_str = String(node["id"])
            m = match(r"SO_(\d+)$", id_str)

            if m !== nothing
                full_id = parse(Int, m.captures[1])
                label = Symbol(node["lbl"])
                short_id = UInt16(length(arr_by_short) + 1)

                push!(arr_by_short, (full_id, label))
                push!(arr_by_full, (full_id, short_id, label))
                push!(arr_by_label, (label, short_id, full_id))

                if haskey(node, "meta") &&
                   haskey(node["meta"], "definition") &&
                   haskey(node["meta"]["definition"], "val")
                    def_dict[full_id] = String(node["meta"]["definition"]["val"])
                end
            end
        end
    end

    sort!(arr_by_full, by = x -> x[1])
    sort!(arr_by_label, by = x -> x[1])

    lookup_obj = SOTermLookup(Tuple(arr_by_short), Tuple(arr_by_full), Tuple(arr_by_label))

    return lookup_obj, def_dict
end

const SO_TERMS, SO_DEFS = _build_lookups()
export SO_TERMS, SO_DEFS

end
