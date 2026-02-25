/*
 * BV lattice invariant -- Magma implementation
 * Portable output format matching bv.py (Sage) and bv.gp (PARI/GP).
 *
 * Uses PARI's qfminim convention: v^T * gram * v <= d.
 * Magma's LatticeWithGram + ShortVectors uses the same convention.
 *
 * Load: load "bv.m";
 */

function LatticeGraph(gram, d)
    L := LatticeWithGram(gram);
    S := ShortVectors(L, d);
    if #S eq 0 then
        return Matrix(Integers(), 0, 0, []);
    end if;
    m := #S;
    R := Matrix(GF(2), m, m,
         [InnerProduct(S[i][1], S[j][1]) : i in [1..m], j in [1..m]]);
    return ChangeRing(R, Integers());
end function;

// Column signature: sorted [value, multiplicity] pairs from a multiset
function ColSig(col)
    return Sort([<v, m> : v -> m in {* x : x in col *}]);
end function;

// BV invariant -- returns [ <sig, count>, ... ]
function BV(gram, d)
    G := LatticeGraph(gram, d);
    m := Nrows(G);
    if m eq 0 then return []; end if;
    p := NextPrime(m);
    Gp := ChangeRing(G, GF(p));
    S := Gp^2;
    cols := [ColSig([Integers()!S[i,j] : i in [1..m]]) : j in [1..m]];
    return Sort([<sig, cnt> : sig -> cnt in {* c : c in cols *}]);
end function;

// Portable polynomial hash matching bv.py HBV_poly
function HBV_poly(bv)
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
end function;

// Portable XOR-multiply hash matching bv.py HBV_xor
function HBV_xor(bv)
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
end function;

function HBV(gram, d)
    return HBV_poly(BV(gram, d));
end function;
