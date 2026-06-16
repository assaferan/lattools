# Portable BV invariant implementations

Reference implementations of the BV lattice invariant and its hash in
Sage/Python, PARI/GP, and Magma. They mirror the fast C routine
(`fast_marked_HBV` in `../eqfminim.c`, compiled to `../eqfminim.so`), but are
written to be readable and to produce identical output across all three
systems.

## Files

| File                 | Purpose                                                       |
|----------------------|---------------------------------------------------------------|
| `bv.py`              | Sage/Python implementation (with doctests)                    |
| `bv.gp`              | PARI/GP implementation                                        |
| `bv.m`               | Magma implementation                                          |
| `compute_hashes.py`  | Generate test-data lines `[gram_flat]:D:HBV_poly:HBV_xor`     |
| `bench.py`           | Per-phase benchmark (qfminim, graph, square, cols, hashes)    |
| `TIMINGS.md`         | Benchmark results across C, Magma, Sage, and GP               |

Each of `bv.py`, `bv.gp`, and `bv.m` exposes the same functions:
`graph`/`LatticeGraph`, `BV`, `HBV_poly`, `HBV_xor`, and `HBV`.

## The BV invariant

`BV(gram, d)` takes a positive-definite Gram matrix and a short-vector bound
`d`:

1. enumerate short vectors `v` with `0 < v^T gram v <= d`, one per `{v, -v}`
   pair (PARI `qfminim` convention);
2. build the mod-2 adjacency matrix `G = S gram S^T mod 2` of those vectors;
3. square it (`G^2` over a prime `p > #short vectors`, so the reduction is a
   no-op);
4. return the canonical multiset of column signatures of `G^2`.

The result is invariant under unimodular change of basis (isometry), so it
distinguishes many non-isometric lattices.

## The two hashes

`BV` returns a nested tuple; the `HBV_*` functions fold it down to a single
integer for fast comparison and storage. Both run over the same integer stream
of the canonical form, differing only in the mixing step.

### `HBV_xor`: the hash used by the C code

    h = 13282407956253574712
    for each x:  h = ((h ^ x) * 1111111111111111111) mod 2^64

This is the hash function used by the C routine in `../eqfminim.c` (init
constant `13282407956253574712`, multiplier `1111111111111111111`, unsigned
64-bit). Reproducing it needs unsigned 64-bit arithmetic, hence the explicit
`mod 2^64` in the GP and Magma versions.

Note: `HBV_xor` reuses the C hash's mixing primitive and constants over the
portable `BV` canonical form. The C routine `fast_marked_HBV` also folds in
per-vector marks and aggregates per-row hashes, so it is not a bit-for-bit
reproduction of the C output.

### `HBV_poly`: the portable default

    h = 0
    for each x:  h = (h * 1000003 + x) mod (2^61 - 1)

All intermediate values stay below the Mersenne prime `2^61 - 1`, so this
reproduces with plain signed 64-bit arithmetic in any language, with no
unsigned wraparound. `HBV(gram, d) = HBV_poly(BV(gram, d))` is the default.

## Usage

Generate test data from a file of Gram matrices (one flat `[a,b,...]`
row-major matrix per line):

    sage -python alt/compute_hashes.py input_file [D]

Benchmark per phase:

    sage -python alt/bench.py [test_data_file] [D]

To run the cross-implementation tests, see `../tests/README.md`.
