\\
\\ Test BV invariant -- PARI/GP
\\ Run: sage -gp -q < tests/test_bv.gp   (from repo root)
\\
\\ To use a custom test data file, set test_data_file before loading:
\\   test_data_file = "path/to/data.txt"; \r tests/test_bv.gp
\\ Defaults to tests/test_bv_data.txt.
\\

default(parisize, 1024000000); \\ 1 GB stack
\r alt/bv.gp

\\ ---- Load test data ----

if(type(test_data_file) != "t_STR", test_data_file = "tests/test_bv_data.txt");

{
my(lines, parts, entries, n, gram, d, ep, ex, td = List());
lines = externstr(Str("grep -v '^#' ", test_data_file, " | grep -v '^[[:space:]]*$'"));
for(k = 1, #lines,
  parts = strsplit(lines[k], ":");
  entries = eval(parts[1]);
  n = sqrtint(#entries);
  if(n*n != #entries, error("entry count not a perfect square: ", #entries));
  gram = matrix(n, n, i, j, entries[n*(i-1)+j]);
  d = eval(parts[2]);
  ep = eval(parts[3]);
  ex = eval(parts[4]);
  listput(td, [gram, d, ep, ex]);
);
test_data = Vec(td);
}

\\ ---- Tests ----

print("============================================================");
print("BV invariant test -- GP");
print("============================================================");

{
my(ok = 1, gram, D, ep, ex, bv, hp, hx, n, m,
   t0, t1, t2, tot_bv = 0,
   buckets = Map(), bk, bt, bc);
t0 = getabstime();
for(i = 1, #test_data,
  gram = test_data[i][1];
  D = test_data[i][2];
  ep = test_data[i][3];
  ex = test_data[i][4];
  n = #gram;
  if(gram != gram~, error("Matrix ", i, " not symmetric"));
  m = #qfminim(gram, D)[3];
  t1 = getabstime();
  bv = BV(gram, D);
  t2 = getabstime();
  hp = HBV_poly(bv);
  hx = HBV_xor(bv);
  tot_bv += t2 - t1;
  bk = ceil(log(max(m, 1))/log(2));
  if(!mapisdefined(buckets, bk),
    mapput(buckets, bk, [0, 0]));
  bt = mapget(buckets, bk);
  mapput(buckets, bk, [bt[1] + t2 - t1, bt[2] + 1]);
  if(hp != ep || hx != ex,
    print("FAIL: Matrix ", i, " (", n, "x", n, ", m=", m, "): poly = ", hp, "  xor = ", hx,
          "  (BV ", t2-t1, "ms)");
    ok = 0;
  ,
    print("  Matrix ", i, " (", n, "x", n, ", m=", m, "): BV ", t2-t1, "ms");
  );
);
print("  Total: ", getabstime() - t0, "ms  (BV ", tot_bv, "ms)");
my(keys = vecsort(Vec(mattranspose(Mat(buckets))[1,])));
for(j = 1, #keys,
  bt = mapget(buckets, keys[j]);
  print("  m <= 2^", keys[j], ": ", bt[2], " matrices, BV ", bt[1], "ms");
);
if(ok, print("PASS: all hashes match expected values (cross-implementation verified)"));
}

\\ ---- C-level fast_marked_HBV comparison ----
\\ Load the C shared library (run from repo root where eqfminim.so lives)
\r tools.gp
\r rs.gp
\r eqfminim.gp

print("");
print("============================================================");
print("C-level fast_marked_HBV (eqfminim.so) -- timing comparison");
print("============================================================");

{
my(t0 = getabstime(), gram, D, h, m, t1, n, tot = 0,
   buckets = Map(), bk, bt, keys);
for(i = 1, #test_data,
  gram = test_data[i][1];
  D = test_data[i][2];
  n = #gram;
  m = #qfminim(gram, D)[3];
  t1 = getabstime();
  h = fast_marked_HBV(gram, [], D);
  my(dt = getabstime() - t1);
  tot += dt;
  bk = ceil(log(max(m, 1))/log(2));
  if(!mapisdefined(buckets, bk),
    mapput(buckets, bk, [0, 0]));
  bt = mapget(buckets, bk);
  mapput(buckets, bk, [bt[1] + dt, bt[2] + 1]);
  print("  Matrix ", i, " (", n, "x", n, ", m=", m, "): ", dt, "ms");
);
print("  Total: ", getabstime() - t0, "ms");
keys = vecsort(Vec(mattranspose(Mat(buckets))[1,]));
for(j = 1, #keys,
  bt = mapget(buckets, keys[j]);
  print("  m <= 2^", keys[j], ": ", bt[2], " matrices, ", bt[1], "ms");
);
}
