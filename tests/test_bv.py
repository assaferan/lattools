"""
Test BV invariant -- Sage
Run: sage -python tests/test_bv.py [test_data_file]   (from repo root)

Uses PARI's qfminim convention (v^T * M * v <= d) to match test_bv.gp and test_bv.m.
Defaults to tests/test_bv_data.txt if no file is specified.
"""
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'alt'))
from bv import BV, HBV_poly, HBV_xor, graph
from sage.all import Matrix, ZZ

def load_test_data(path):
    """Parse test_bv_data.txt -> list of (gram, D, expected_poly, expected_xor)."""
    data = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            gram_str, d_str, poly_str, xor_str = line.split(':')
            entries = list(map(int, gram_str.strip().strip('[]').split(',')))
            n = int(len(entries)**0.5)
            assert n * n == len(entries), f"entry count {len(entries)} is not a perfect square"
            gram = Matrix(ZZ, n, n, entries)
            data.append((gram, int(d_str.strip()), int(poly_str.strip()), int(xor_str.strip())))
    return data

if len(sys.argv) > 1:
    data_path = sys.argv[1]
else:
    data_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_bv_data.txt')
test_data = load_test_data(data_path)

print("=" * 60)
print("BV invariant test -- Sage")
print("=" * 60)

from math import ceil, log2

t0 = time.time()
total_bv = 0
bv_by_bucket = {}
ok = True
for i, (gram, D, exp_poly, exp_xor) in enumerate(test_data):
    assert gram == gram.T, f"Matrix {i+1} is not symmetric"
    m = graph(gram, D).ncols()
    t1 = time.time()
    bv = BV(gram, D)
    t2 = time.time()
    hp = HBV_poly(bv)
    hx = HBV_xor(bv)
    total_bv += t2 - t1
    b = ceil(log2(max(m, 1)))
    bv_by_bucket.setdefault(b, [0.0, 0])
    bv_by_bucket[b][0] += t2 - t1
    bv_by_bucket[b][1] += 1
    if hp != exp_poly or hx != exp_xor:
        print(f"FAIL: Matrix {i+1} ({gram.nrows()}x{gram.ncols()}, m={m}): poly = {hp}  xor = {hx}  (BV {t2-t1:.4f}s)")
        ok = False
    else:
        print(f"  Matrix {i+1} ({gram.nrows()}x{gram.ncols()}, m={m}): BV {t2-t1:.4f}s")
print(f"  Total: {time.time() - t0:.3f}s  (BV {total_bv:.3f}s)")
for b in sorted(bv_by_bucket):
    t, cnt = bv_by_bucket[b]
    print(f"  m <= 2^{b}: {cnt} matrices, BV {t:.3f}s")

if ok:
    print("PASS: all hashes match expected values (cross-implementation verified)")
