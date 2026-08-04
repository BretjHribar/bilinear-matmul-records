# Supplemental `<8,12,12>` rank-705 construction

FMM-Lille lists rank 720 for `<8,12,12>` and explicitly constructs it as the
tensor product of rank-15 `<2,3,3>` with rank-48 `<4,4,4>`:

https://fmm.univ-lille.fr/8x12x12.html

The retained generator substitutes the exact GF(2) rank-47 `<4,4,4>` factor,
then materializes the tensor product in standard row-major coordinates. This
gives `15 * 47 = 705`. Both independent verifiers return zero residual.

This is a characteristic-two comparison, not a claimed arbitrary-ring bound.
