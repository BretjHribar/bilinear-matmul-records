# Exact GF(2) rank-705 scheme for `<8,12,12>`

This requested supplemental record contains a complete GF(2) decomposition of
`(8x12) * (12x12) -> (8x12)` using 705 scalar multiplications. It is not used
as the square-matrix research result.

FMM-Lille lists rank 720 as `15 * 48`. Field-aware substitution of the exact
rank-47 GF(2) `<4,4,4>` factor gives `15 * 47 = 705`.

```sh
python3 tensor_product.py
python3 verify_8x12x12_rank705_gf2.py matmul_8x12x12_rank705_gf2.txt
g++ -O3 -std=c++20 verify_8x12x12_rank705_gf2.cpp -o verify812
./verify812 matmul_8x12x12_rank705_gf2.txt
```
