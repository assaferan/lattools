/*
 * BV lattice invariant -- Magma implementation
 * Portable output format matching bv.py (Sage) and bv.gp (PARI/GP).
 *
 * Uses PARI's qfminim convention: v^T * gram * v <= d.
 * Magma's LatticeWithGram + ShortVectors uses the same convention.
 *
 * Attach: Attach("bv.m");   // package file: HBV_poly, HBV_xor, HBV, BV are intrinsics
 */

// Build the adjacency matrix from short vectors of a lattice.
// Returns an integer 0/1 matrix: R = C * gram * C^T mod 2,
// where rows of C are coordinates of short vectors (one per {v,-v} pair).
// Uses GF(2) matrix multiply (bit-packed) instead of per-pair inner products.
function LatticeGraph(gram, d)
    L := LatticeWithGram(gram);
    S := ShortVectors(L, d);
    if #S eq 0 then
        return Matrix(Integers(), 0, 0, []);
    end if;
    m := #S;
    n := Nrows(gram);
    C := Matrix(GF(2), m, n, &cat[Eltseq(S[i][1]) : i in [1..m]]);
    R := C * ChangeRing(gram, GF(2)) * Transpose(C);
    return ChangeRing(R, Integers());
end function;

// Column signature: sorted [value, multiplicity] pairs from a multiset
function ColSig(col)
    return Sort([<v, m> : v -> m in {* x : x in col *}]);
end function;

// BV invariant -- returns [ <sig, count>, ... ]
// Squares over Z instead of GF(p): since G is 0/1 and p > m,
// entries of G^2 are at most m < p, so mod p is a no-op.
// Uses row iteration (S is symmetric) for cache-friendly access.
intrinsic BV(gram::Mtrx, d::RngIntElt) -> SeqEnum
{BV invariant: sorted sequence of <signature, count> tuples for the Gram matrix}
    G := LatticeGraph(gram, d);
    m := Nrows(G);
    if m eq 0 then return []; end if;
    S := G^2;
    cols := [ColSig(Eltseq(r)) : r in Rows(S)];
    return Sort([<sig, cnt> : sig -> cnt in {* c : c in cols *}]);
end intrinsic;

// Portable polynomial hash matching bv.py HBV_poly
intrinsic HBV_poly(bv::SeqEnum[Tup]) -> RngIntElt
{Poly based hash of a sequence of integers}
    M := 2^61 - 1;
    h := 0;
    for entry in bv do
        for vc in entry[1] do
            h := (h * 1000003 + vc[1]) mod M;
            h := (h * 1000003 + vc[2]) mod M;
        end for;
        h := (h * 1000003 + entry[2]) mod M;
    end for;
    return h;
end intrinsic;

// Portable XOR-multiply hash matching bv.py HBV_xor
intrinsic HBV_xor(bv::SeqEnum[Tup]) -> RngIntElt
{Xor based hash of a sequence of integers}
    M := 2^64;
    MULT := 1111111111111111111;
    h := 13282407956253574712;
    for entry in bv do
        for vc in entry[1] do
            h := (BitwiseXor(h, vc[1]) * MULT) mod M;
            h := (BitwiseXor(h, vc[2]) * MULT) mod M;
        end for;
        h := (BitwiseXor(h, entry[2]) * MULT) mod M;
    end for;
    return h;
end intrinsic;

intrinsic HBV(gram::Mtrx, d::RngIntElt) -> RngIntElt
{The main BV-based hash for a Gram matrix}
    return HBV_poly(BV(gram, d));
end intrinsic;
