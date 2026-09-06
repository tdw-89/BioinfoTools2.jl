module SOTerms

using JSON

"""Directory holding the package's bundled data files."""
const ASSETS = joinpath(pkgdir(SOTerms), "assets")

"""Snapshot of the SOFA ontology, refreshed by `./update_assets.sh`."""
const SO_JSON = joinpath(ASSETS, "SOFA.json")
include_dependency(SO_JSON)

"""
Stack-allocated lookup for the SO terms (SOFA subset) a GFF3 file's third column
can hold.

Each field is a `Tuple` of entries sorted by its own key, so a lookup is a
binary search over immutable data:
- `by_short`: `(full_id, label)` indexed directly by short code.
- `by_full`: `(full_id, short_code, label)` sorted by `full_id`.
- `by_label`: `(label, short_code, full_id)` sorted by `label`.
"""
struct SOTermLookup{T1<:Tuple,T2<:Tuple,T3<:Tuple}
    by_short::T1
    by_full::T2
    by_label::T3
end

"""
Binary search `entries` — sorted by the first element of each entry — for
`key`, returning that entry's remaining elements, or `nothing` when absent.
"""
function _search(entries::Tuple, key)
    left, right = 1, length(entries)
    while left <= right
        middle = (left + right) >>> 1
        found = entries[middle][1]
        if found == key
            return entries[middle][2:end]
        elseif found < key
            left = middle + 1
        else
            right = middle - 1
        end
    end
    return nothing
end

"""Look a term up by its short code, returning `(full_id, label)`."""
Base.getindex(lookup::SOTermLookup, short::UInt16) =
    1 <= short <= length(lookup.by_short) ? lookup.by_short[short] : nothing

"""Look a term up by its full ontology ID, returning `(short_code, label)`."""
Base.getindex(lookup::SOTermLookup, full::Int) = _search(lookup.by_full, full)

"""Look a term up by its label, returning `(short_code, full_id)`."""
Base.getindex(lookup::SOTermLookup, label::Symbol) = _search(lookup.by_label, label)

Base.getindex(lookup::SOTermLookup, label::AbstractString) = lookup[Symbol(label)]

Base.show(io::IO, lookup::SOTermLookup) =
    print(io, "SOTermLookup(", length(lookup.by_short), " terms)")

"""
Parse [`SO_JSON`](@ref) into an [`SOTermLookup`](@ref) and a `full_id =>
definition` dictionary. Short codes are assigned in the order the ontology
lists its `SO_<digits>` nodes, so they are stable for a given snapshot.
"""
function _build_lookups()
    raw_json = JSON.parse(read(SO_JSON, String))

    by_short = Tuple{Int,Symbol}[]
    by_full = Tuple{Int,UInt16,Symbol}[]
    by_label = Tuple{Symbol,UInt16,Int}[]
    definitions = Dict{Int,String}()

    for node in raw_json["graphs"][1]["nodes"]
        haskey(node, "id") || continue
        matched = match(r"SO_(\d+)$", String(node["id"]))
        matched === nothing && continue

        full_id = parse(Int, matched.captures[1])
        label = Symbol(node["lbl"])
        short_code = UInt16(length(by_short) + 1)

        push!(by_short, (full_id, label))
        push!(by_full, (full_id, short_code, label))
        push!(by_label, (label, short_code, full_id))

        meta = get(node, "meta", nothing)
        definition = meta === nothing ? nothing : get(meta, "definition", nothing)
        if definition !== nothing && haskey(definition, "val")
            definitions[full_id] = String(definition["val"])
        end
    end

    sort!(by_full; by = first)
    sort!(by_label; by = first)

    return SOTermLookup(Tuple(by_short), Tuple(by_full), Tuple(by_label)), definitions
end

# The package-wide SO term lookup and its `full_id => definition` companion,
# built from `SO_JSON` at load time.
const SO_TERMS, SO_DEFS = _build_lookups()

export SO_TERMS, SO_DEFS

end
