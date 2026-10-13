"""
    clamp_idx(i, n, Val(mode)) -> Int

Fold the index `i` into `1:n` for the boundary extension `mode` of a signal of
length `n`:

- `:symmetric`: whole-sample symmetric, the edge sample is not repeated
  (`x[2], x[1], x[2], x[3], …`);
- `:periodic`: `mod1(i, n)`;
- `:zero`: `i` unchanged; the caller treats indices outside `1:n` as zeros.

The mode is a `Val` so each caller is specialised, with no per-sample dispatch.

```jldoctest
julia> using MultiscaleFilterBanks: clamp_idx

julia> [clamp_idx(i, 4, Val(:symmetric)) for i in -1:6]
8-element Vector{Int64}:
 3
 2
 1
 2
 3
 4
 3
 2
```
"""
function clamp_idx end

@inline function clamp_idx(i::Int, n::Int, ::Val{:symmetric})::Int
    1 <= i <= n && return i
    n == 1 && return 1
    p = 2 * (n - 1)
    m = mod(i - 1, p)
    return n - abs(m - (n - 1))
end

@inline function clamp_idx(i::Int, n::Int, ::Val{:periodic})::Int
    return mod(i - 1, n) + 1
end

@inline function clamp_idx(i::Int, n::Int, ::Val{:zero})::Int
    return i
end
