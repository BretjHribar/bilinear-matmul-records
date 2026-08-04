# Field-aware tensor substitution: `<8,8,8>` in rank 329

## Result

The certificate proves

`R_GF(2)(<8,8,8>) <= 329`.

The FMM-Lille page audited on 2026-08-04 lists rank 336 and describes its
construction as the tensor product of rank-7 `<2,2,2>` with rank-48
`<4,4,4>`:

https://fmm.univ-lille.fr/8x8x8.html

## Combined technique

The change is a field-aware leaf substitution inside the catalogue's existing
tensor-product architecture:

1. retain the exact rank-7 Strassen `<2,2,2>` factor;
2. replace the general rank-48 `<4,4,4>` factor with an exact rank-47 GF(2)
   decomposition;
3. materialize the Kronecker product into standard row-major 8 by 8 masks;
4. exact-gate the complete result rather than relying on rank arithmetic.

For matrices `A` of shape `a*b` and `B` of shape `c*d`, the generator maps a
pair of factor coordinates `(i,j)` and `(k,l)` to
`(i*c+k, j*d+l)`. Applying that mapping independently to U, V and W turns
every pair of source terms into one rank-one term of the product tensor. There
are exactly `47 * 7 = 329` terms.

The two source certificates are retained under `sources/`; `tensor_product.py`
recreates the published bytes deterministically.

## Verification

The Python verifier accumulates W masks in independent `(U-coordinate,
V-coordinate)` residual slices. The C++ verifier toggles individual flattened
Brent coefficients in a hash set. Both check rank, factor bounds, duplicates,
all target coefficients, and require zero residual.

This is a characteristic-two result. It numerically beats the Lille entry over
GF(2), but does not assert rank 329 over arbitrary rings or characteristic zero.

Public source lineage:

- https://github.com/tungsten-lang/tungsten/blob/main/benchmarks/matmul/metaflip/matmul_4x4_rank47_d450_gf2.txt
- https://github.com/tungsten-lang/tungsten/blob/main/benchmarks/matmul/metaflip/matmul_2x2_rank7_strassen_gf2.txt
