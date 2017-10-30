# Implementation of the BRW polynomial based hash functions BRW128 and BRW256.

This corresponds to the paper: "Evaluating Bernstein-Rabin-Winograd Polynomials" authored by:

    Sebati Ghosh,                                        Indian Statistical Institute;

    Palash Sarkar,                                       Indian Statistical Institute;

This package contains the implementations for BRW128 and BRW256 for t = 2, 3, 4 and 5, where t is a parameter to the algorithm EvalBRW, which is the main part of the hash functions and is a non-recursive algorithm computing BRW polynomials on any number of field elements.

In each package, 

    the file mult128.h or mult256.h contains routines for basic field multiplications and routines to calculate 
    unreducedBRW on 1 to 2^t-1 field elements.

    the file brw.c contains the routine to evaluate the hash function and the main function.


For compiling any of the packages, we need to use the corresponding flags:

    -mavx2 (to compile the SSE and AVX instructions)

    -mpclmul (for pclmulqdq)

    -O3 (optimisation flag)
    
    -lm (for computing log)

So, for example we use the command:

    gcc brw.c -mavx2 -mpclmul -O3 -lm

    and the executable a.out is produced.
    
For measuring performance of each package we have used the strategy used in http://dx.doi.org/10.1007/978-3-642-21702-9_18
and the corresponding file is timedefs.h.
