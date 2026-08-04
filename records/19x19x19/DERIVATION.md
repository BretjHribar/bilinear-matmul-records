# Derivation of the `<19,19,19>` rank-3,987 certificate

## Starting construction

The construction is support-aware block composition over GF(2). It uses the
public exact rank-47 `<4,4,4>` outer scheme and the block allocations

- `N = (4,5,5,5)`
- `M = (5,5,5,4)`
- `P = (4,5,5,5)`

with the standard exact leaf bank through dimension five. With the unmodified
AlphaEvolve rank-93 `<5,5,5>` leaf, exact embedding, zero-term deletion, and
parity cancellation give rank 3,993.

## Six-term improvement

In the AlphaEvolve leaf, the output mask `27263834` is the unique output factor
of multiplicity three. As a 5 by 5 binary matrix it is rank one:

`27263834 = (10011)_2 (11010)_2^T`.

An exact matrix-multiplication isotropy was applied to the leaf. If `A`, `B`,
and `C` are invertible 5 by 5 binary matrices, each leaf term transforms as

- `U -> A U B^-1`
- `V -> B V C^-1`
- `W -> A^-T W C^T`.

Here `B` is the identity; `A` and `C` map the two vectors of the repeated
rank-one output factor to local coordinate 4 (zero based). This preserves the
leaf tensor exactly. That coordinate is the clipped row/column at the outer
construction's size-four block boundaries. Three copies vanish wherever the
old presentation lost two. Across the affected recursive slots this removes
six additional scalar products, producing rank 3,987.

As a reconstruction sanity check, all other 24 placements of the repeated
factor among the 25 local coordinate pairs were scanned: each produced rank
4,005 with the deterministic complementary bases used here. The doubly clipped
pair `(4,4)` uniquely produced 3,987.

## Exact gate

The retained Python and C++ runs independently report:

- 3,987 parsed nonzero, in-range, distinct triples;
- 4,659,259 sparse coefficient toggles;
- 6,859 target coefficients;
- zero residual coefficients.

The exhaustive sparse check is algebraically identical to checking the full
`361^3 = 47,045,881` Brent coefficient space: coefficients never toggled are
zero and need not be materialized.

Public parent sources used in reconstruction:

- Tungsten rank-3,993 block construction and block composer
- Tungsten AlphaEvolve rank-93 5 by 5 leaf
- Tungsten rank-47 4 by 4 outer and exact leaf bank

These are in `tungsten-lang/tungsten/benchmarks/matmul/metaflip`.
