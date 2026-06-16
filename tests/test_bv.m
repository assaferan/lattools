/*
 * Test BV invariant -- Magma
 * Run: magma -b tests/test_bv.m [test_data_file]   (from repo root)
 *
 * Defaults to tests/test_bv_data.txt if no file is specified.
 */

SetColumns(0);

load "alt/bv.m";

function LoadTestData(path)
    data := [* *];
    F := Open(path, "r");
    while true do
        line := Gets(F);
        if IsEof(line) then break; end if;
        line := StripWhiteSpace(line);
        if #line eq 0 or line[1] eq "#" then continue; end if;
        parts := Split(line, ":");
        gram_str := StripWhiteSpace(parts[1]);
        // strip [ and ]
        gram_str := gram_str[2..#gram_str-1];
        entries := [StringToInteger(StripWhiteSpace(x)) : x in Split(gram_str, ",")];
        n := Isqrt(#entries);
        assert n^2 eq #entries;
        gram := Matrix(Integers(), n, n, entries);
        d := StringToInteger(StripWhiteSpace(parts[2]));
        poly_h := StringToInteger(StripWhiteSpace(parts[3]));
        xor_h := StringToInteger(StripWhiteSpace(parts[4]));
        Append(~data, <gram, d, poly_h, xor_h>);
    end while;
    return data;
end function;

// ---- Tests ----

if assigned args and #args ge 1 then
    data_file := args[1];
else
    data_file := "tests/test_bv_data.txt";
end if;
data := LoadTestData(data_file);

printf "============================================================\n";
printf "BV invariant test -- Magma\n";
printf "============================================================\n";

t0 := Cputime();
tot_bv := 0.0;
bv_by_bucket := AssociativeArray(); cnt_by_bucket := AssociativeArray();
ok := true;
for i -> entry in data do
    gram := entry[1];
    D := entry[2];
    exp_poly := entry[3];
    exp_xor := entry[4];
    n := Nrows(gram);
    assert gram eq Transpose(gram);
    t1 := Cputime();
    bv := BV(gram, D);
    t2 := Cputime();
    hp := HBV_poly(bv);
    hx := HBV_xor(bv);
    tot_bv +:= t2 - t1;
    m := #ShortVectors(LatticeWithGram(gram), D);
    b := Ceiling(Log(2, Max(m, 1)));
    if not IsDefined(bv_by_bucket, b) then
        bv_by_bucket[b] := 0.0; cnt_by_bucket[b] := 0;
    end if;
    bv_by_bucket[b] +:= t2 - t1; cnt_by_bucket[b] +:= 1;
    if hp ne exp_poly or hx ne exp_xor then
        printf "FAIL: Matrix %o (%ox%o, m=%o): poly = %o  xor = %o  (BV %os)\n",
               i, n, n, m, hp, hx, t2-t1;
        ok := false;
    else
        printf "  Matrix %o (%ox%o, m=%o): BV %os\n", i, n, n, m, t2-t1;
    end if;
end for;
printf "  Total: %os  (BV %os)\n", Cputime(t0), tot_bv;
for b in Sort(SetToSequence(Keys(bv_by_bucket))) do
    printf "  m <= 2^%o: %o matrices, BV %os\n",
           b, cnt_by_bucket[b], bv_by_bucket[b];
end for;

if ok then
    printf "PASS: all hashes match expected values (cross-implementation verified)\n";
end if;

quit;
