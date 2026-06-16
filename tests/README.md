# Running BV invariant tests

All commands run from the repository root.

## Sage

    sage -python tests/test_bv.py

## PARI/GP

    sage -gp -q < tests/test_bv.gp

## Magma

    magma -b tests/test_bv.m

## Custom test data

Each test accepts an optional path to a different data file:

    sage -python tests/test_bv.py tests/test_bv_data_hard.txt
    magma -b tests/test_bv.m tests/test_bv_data_hard.txt

For GP, set the variable before loading:

    echo 'test_data_file = "tests/test_bv_data_hard.txt"; \r tests/test_bv.gp' | sage -gp -q
