# Exact GF(2) bilinear matrix-multiplication records

This repository contains complete decimal-mask certificates and independent
Brent-equation verifiers for two exact non-commutative tensor decompositions:

| Tensor | Rank | Certificate |
|---|---:|---|
| `<19,19,19>` | 3,987 | [`records/19x19x19`](records/19x19x19/) |
| `<15,15,17>` | 2,256 | [`records/15x15x17`](records/15x15x17/) |

Each term is `U V W`, where the three nonzero decimal integers encode row-major
binary coefficient masks. The sum of `U_t tensor V_t tensor W_t` is evaluated
over GF(2).

## Verify

```sh
cd records/19x19x19
python3 verify_19x19x19_rank3987_gf2.py matmul_19x19x19_rank3987_gf2.txt
```

Both implementations expand every sparse rank-one term, toggle every induced
Brent coefficient, toggle all target matrix-multiplication coefficients, and
require an empty residual. They also reject zero/out-of-range factors,
rank-header mismatches, and duplicate triples. Each record directory includes
the C++ command, retained logs, hashes, metadata, and derivation notes.

The record comparison in the frozen research notes was field- and model-aware:
commutative algorithms are not comparators for these ordinary non-commutative
`(U,V,W)` decompositions. A public-record statement cannot exclude unpublished
or uncatalogued decompositions.
