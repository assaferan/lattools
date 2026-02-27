"""
Benchmark BV computation with per-phase timing.

Usage: sage -python alt/bench.py [test_data_file] [D]

Reads test data in the format used by tests/test_bv_data.txt and times
each phase of the BV computation separately:

  qfminim  -- short vector enumeration (PARI)
  graph    -- S * gram * S^T mod 2
  square   -- G^2 mod p (matrix squaring)
  cols     -- column signature computation
  HBV_poly -- polynomial hash
  HBV_xor  -- XOR-multiply hash
"""
import sys
import time
from collections import Counter

from sage.all import GF, Matrix, ZZ, next_prime, pari


def bench_one(gram, d):
    """Return per-phase wall times (seconds) for a single Gram matrix."""
    t = {}

    # Phase 1: qfminim
    t0 = time.perf_counter()
    S = pari(gram).qfminim(d)[2]
    t["qfminim"] = time.perf_counter() - t0

    # Phase 2: graph (S * gram * S^T mod 2)
    t0 = time.perf_counter()
    S = Matrix(ZZ, S).T
    m = S.nrows()
    if m == 0:
        t["graph"] = time.perf_counter() - t0
        t["square"] = 0.0
        t["cols"] = 0.0
        t["HBV_poly"] = 0.0
        t["HBV_xor"] = 0.0
        return t
    G = (S * gram * S.T).change_ring(GF(2))
    t["graph"] = time.perf_counter() - t0

    # Phase 3: square (G^2 mod p)
    t0 = time.perf_counter()
    p = next_prime(m)
    Gp = G.change_ring(GF(p))
    G2 = Gp * Gp
    t["square"] = time.perf_counter() - t0

    # Phase 4: column signatures
    t0 = time.perf_counter()
    cols = []
    for j in range(m):
        col = tuple(int(G2[i, j]) for i in range(m))
        sig = tuple(sorted(Counter(col).items()))
        cols.append(sig)
    bv = tuple(sorted(Counter(cols).items()))
    t["cols"] = time.perf_counter() - t0

    # Phase 5: HBV_poly
    t0 = time.perf_counter()
    M = (1 << 61) - 1
    h = 0
    for sig, count in bv:
        for value, multiplicity in sig:
            h = (h * 1000003 + int(value)) % M
            h = (h * 1000003 + int(multiplicity)) % M
        h = (h * 1000003 + int(count)) % M
    t["HBV_poly"] = time.perf_counter() - t0

    # Phase 6: HBV_xor
    t0 = time.perf_counter()
    MASK = (1 << 64) - 1
    MULT = 1111111111111111111
    h = 13282407956253574712
    for sig, count in bv:
        for value, multiplicity in sig:
            h = ((h ^ int(value)) * MULT) & MASK
            h = ((h ^ int(multiplicity)) * MULT) & MASK
        h = ((h ^ int(count)) * MULT) & MASK
    t["HBV_xor"] = time.perf_counter() - t0

    return t


def fmt(seconds):
    """Format seconds as a human-readable string."""
    if seconds < 0.001:
        return f"{seconds*1e6:.0f}us"
    if seconds < 1.0:
        return f"{seconds*1e3:.0f}ms"
    return f"{seconds:.1f}s"


def main():
    data_file = sys.argv[1] if len(sys.argv) > 1 else "tests/test_bv_data.txt"
    D = int(sys.argv[2]) if len(sys.argv) > 2 else 4

    with open(data_file) as f:
        lines = [l.strip() for l in f if l.strip() and l.strip().startswith('[')]

    phases = ["qfminim", "graph", "square", "cols", "HBV_poly", "HBV_xor"]
    totals = {k: 0.0 for k in phases}

    print(f"Benchmarking {len(lines)} Gram matrices, D = {D} ...")
    wall_start = time.perf_counter()

    for line in lines:
        gram_str, *_ = line.split(":")
        entries = list(map(int, gram_str.strip().strip("[]").split(",")))
        n = int(len(entries) ** 0.5)
        gram = Matrix(ZZ, n, n, entries)
        t = bench_one(gram, D)
        for k in phases:
            totals[k] += t[k]

    wall_total = time.perf_counter() - wall_start

    # Print results
    print()
    header = "| Phase    | Time    | % of BV |"
    sep    = "|----------|---------|---------|"
    print(header)
    print(sep)
    bv_total = totals["qfminim"] + totals["graph"] + totals["square"] + totals["cols"]
    for k in phases[:4]:
        pct = 100 * totals[k] / bv_total if bv_total > 0 else 0
        print(f"| {k:<8s} | {fmt(totals[k]):>7s} | {pct:5.1f}%  |")
    print(sep)
    print(f"| {'BV':8s} | {fmt(bv_total):>7s} | 100.0%  |")
    print(f"| {'HBV_poly':8s} | {fmt(totals['HBV_poly']):>7s} |         |")
    print(f"| {'HBV_xor':8s} | {fmt(totals['HBV_xor']):>7s} |         |")
    print(sep)
    print(f"| {'Total':8s} | {fmt(wall_total):>7s} |         |")
    print()


if __name__ == "__main__":
    main()
