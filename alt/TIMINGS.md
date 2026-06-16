# BV invariant timings

399 Gram matrices (37 of size 10x10 D=4, 328 of size 12x12 D=4, 34 of size 16x16 D=3).

## Overall

| Implementation       | BV     | HBV_poly | HBV_xor | Total   |
|----------------------|--------|----------|---------|---------|
| C (eqfminim.so)      | —      | —        | —       | 2.0s    |
| Magma                | 5.9s   | 20ms     | 50ms    | 5.9s    |
| Sage                 | 34.3s  | 31ms     | 25ms    | 35.7s   |
| GP (portable)        | 45.2s  | 43ms     | 46ms    | 45.3s   |

The hash functions (HBV_poly, HBV_xor) are negligible.
The C implementation (now compiled with -march=native -O3) is ~3x faster than Magma
and ~18–23x faster than the portable Sage/GP implementations.
