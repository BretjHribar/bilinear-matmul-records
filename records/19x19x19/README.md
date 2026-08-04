# Exact GF(2) rank-3,987 scheme for `<19,19,19>`

This directory contains a complete non-commutative bilinear decomposition of
19 by 19 matrix multiplication over GF(2) using 3,987 scalar multiplications.

## Files

- `matmul_19x19x19_rank3987_gf2.txt` — complete decimal-mask certificate.
- `verify_19x19x19_rank3987_gf2.py` and `.cpp` — exhaustive verifiers.
- `DERIVATION.md` — reconstruction and rank-drop technique.
- `verification_manifest.json` — machine-readable invariants.
- `verification_*.log` — retained successful runs.
- `SHA256SUMS` — file digests.

## Reproduce

```sh
python3 verify_19x19x19_rank3987_gf2.py matmul_19x19x19_rank3987_gf2.txt
g++ -O3 -std=c++20 verify_19x19x19_rank3987_gf2.cpp -o verify19
./verify19 matmul_19x19x19_rank3987_gf2.txt
```

Both checks expand all sparse rank-one contributions and require the complete
Brent residual to be zero.
