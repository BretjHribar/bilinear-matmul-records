# Exact GF(2) rank-329 scheme for `<8,8,8>`

This directory contains a complete ordinary non-commutative bilinear
decomposition of 8 by 8 matrix multiplication over GF(2), using 329 scalar
multiplications.

The live FMM-Lille entry is rank 336 and is explicitly `7 * 48`. This record
uses the same tensor-product architecture with a characteristic-two rank-47
`<4,4,4>` factor, hence `7 * 47 = 329`.

## Reproduce

```sh
python3 tensor_product.py
python3 verify_8x8x8_rank329_gf2.py matmul_8x8x8_rank329_gf2.txt
g++ -O3 -std=c++20 verify_8x8x8_rank329_gf2.cpp -o verify8
./verify8 matmul_8x8x8_rank329_gf2.txt
```

See `DERIVATION.md` for the construction and comparison scope. `sources/`
retains the two exact input certificates needed by the deterministic generator.
