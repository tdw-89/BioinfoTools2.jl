# Implementation Blueprint: In-Memory & Arrow-Backed Methylation Data Storage

## Context & Objectives
We are implementing an in-memory and disk-backed storage format for single-base methylation calls in this package. 

### Key Architectural Decisions
1. **Target Data:** Short-read bisulfite sequencing (or pre-aggregated single-base call) data. 
2. **Data Model:** **Pre-aggregated positional counts** (per-strand, per-genomic site). Read IDs are discarded during parsing to save memory and eliminate string-dictionary overhead.
3. **Memory Layout:** **Struct of Arrays (SoA)** via `StructArrays.jl` for optimal cache locality, fast binary searching (`searchsorted`), and direct compatibility with columnar formats.
4. **Disk / Virtual Memory Layout:** Apache Arrow (`Arrow.jl`) with Zstandard compression. Allows zero-copy memory mapping (`mmap`) of large chromosome files.
5. **Footprint Target:** Exactly **8 bytes per genomic site** (`pos::UInt32` + `payload::UInt32`).

---

## Data Structure & Bit-Packing Specification

Each genomic call is stored as an 8-byte entry composed of two 32-bit values:

```julia
struct AggregatedCall
    pos::UInt32      # 4 bytes: 1-indexed genomic coordinate (partitioned by chromosome)
    payload::UInt32  # 4 bytes: Bit-packed metadata and coverage counts
end
```

### `payload::UInt32` Bit Field Allocation

| Bit Range | Length | Field | Description | Range / Values |
| :--- | :--- | :--- | :--- | :--- |
| `0..11` | 12 bits | `meth_count` | Methylated read count | 0 <= N <= 4095 (clamped) |
| `12..23` | 12 bits | `unmeth_count` | Unmethylated read count | 0 <= N <= 4095 (clamped) |
| `24..25` | 2 bits | `context` | Cytosine context | `0x00` = CpG, `0x01` = CHG, `0x02` = CHH |
| `26` | 2 bits | `strand` | Strand orientation | `0x00` = Forward (`+`), `0x01` = Reverse (`-`), `0x02` = Both, `0x03` = neither/NA |
| `27..31` | 4 bits | *Reserved* | Unused flags | Reserved for future quality filters/SNPs |

---

## Instructions for Implementation LLM

> **IMPORTANT CODEBASE CONSISTENCY INSTRUCTION:**
> Before defining new custom bit-shifting routines or bitmask macros, **search the codebase for existing bit-packing modules, types, or helper macros**. Re-use the existing bitwise paradigms, masks, or utility functions present elsewhere in this package to maintain idiomatic consistency across the project.

### Step 1: Audit Existing Bit-Packing Patterns
* Search for where bitpacking is currently used in the package.
* Align the naming conventions, bit-masking macros, or bit-field utilities of `AggregatedCall` with those existing package standards.
* Note that the handling of strand is different (only 4 bits here). Ideally I would like this to be uniform across the library, so if you can change the *other* implementations of strand info (i.e., storing and parsing the info) in bits that would be great. Just make sure this one stays only 4 bits.

### Step 2: Define Core Types & Bitwise Accessors
* Define `AggregatedCall` and constant values for contexts (`CPG`, `CHG`, `CHH`).
* Implement constructor/encoding methods
  *(Ensure read counts exceeding 4095 are safely clamped to `4095` to prevent bit-overflow).*
* Implement inline accessor functions (or extend `Base.getproperty` if consistent with the rest of the package):
  * `get_meth(payload)`
  * `get_unmeth(payload)`
  * `get_depth(payload)`
  * `get_context(payload)`
  * `is_forward(payload)`

### Step 3: Implement StructArray & Search Interface
* Provide constructor functions for initializing a `StructArray{AggregatedCall}` given expected site counts.
* Implement range-query utilities using Julia's binary search on the contiguous `pos` vector:
  ```julia
  function find_calls_in_range(calls::StructArray{AggregatedCall}, start_pos::Integer, stop_pos::Integer)
      # Binary search over calls.pos
  end
  ```

### Step 4: Arrow Read/Write & Memory-Mapping Interface
* Implement export functions to write `StructArray{AggregatedCall}` to an Arrow file:
  ```julia
  write_methylation_arrow(filepath::String, calls::StructArray{AggregatedCall}; compress=:zstd)
  ```
* Implement import functions to read back memory-mapped Arrow tables and re-wrap them into zero-copy `StructArray` instances:
  ```julia
  read_methylation_arrow(filepath::String)::StructArray{AggregatedCall}
  ```

### Step 5: Unit Tests
Write comprehensive unit tests in `test/` covering:
1. Encoding and decoding round-trips for edge cases (0 counts, 4095 counts, >4095 overflow clamping, all contexts, both strands).
2. Range searches across millions of simulated positions.
3. Writing to Arrow with Zstd compression and reading back via `mmap`, verifying equality of data.