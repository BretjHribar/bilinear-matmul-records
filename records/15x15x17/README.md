# Exact GF(2) rank-2256 scheme for `<15,15,17>`

This directory contains a complete non-commutative bilinear decomposition of

`(15×15) · (15×17) → (15×17)`

over `GF(2)` using 2,256 scalar multiplications.

## Files

- `matmul_15x15x17_rank2256_gf2.txt` — the complete decimal-mask certificate.
- `verify_15x15x17_rank2256_gf2.py` — independent exhaustive Python verifier.
- `verify_15x15x17_rank2256_gf2.cpp` — independent exhaustive C++ verifier.
- `verification_manifest.json` — dimensions, format, expected invariants, and provenance.
- `DERIVATION.md` — technique and reconstruction provenance.
- `verification_python.log` and `verification_cpp.log` — successful runs.
- `SHA256SUMS` — SHA-256 digests for the retained files in this directory.

## Certificate format

After comment lines, the first non-comment line is the rank. Each subsequent
line contains three nonzero decimal integers `U V W`. Their binary expansions
encode one rank-one term `U_t ⊗ V_t ⊗ W_t` using row-major coordinates:

- `U`: 225 coordinates for the `15×15` left input;
- `V`: 255 coordinates for the `15×17` right input;
- `W`: 255 coordinates for the `15×17` output.

The certificate header records the source allocation as
`N=(3,4,4,4)`, `M=(4,4,4,3)`, and `P=(4,5,4,4)`.

## Reproduce

Python:

```sh
python3 verify_15x15x17_rank2256_gf2.py \
  matmul_15x15x17_rank2256_gf2.txt
```

C++:

```sh
g++ -O3 -std=c++20 -o verify_15x15x17_rank2256_gf2 \
  verify_15x15x17_rank2256_gf2.cpp
./verify_15x15x17_rank2256_gf2 matmul_15x15x17_rank2256_gf2.txt
```

Each verifier checks the declared and parsed rank, nonzero/in-range factors,
duplicate triples, all `225 × 255 × 255 = 14,630,625` Brent coefficients, and
zero residual against the matrix-multiplication tensor.
