
#ifndef _typedefs_h
#define _typedefs_h 1

struct Param{
    int   runNum;
    int   distnSteps;
    int   gen;
    int   loci;
    int   popsize;
    float mutation;
    float recombination;
    float s;        // strength of selection
    float minFail;  // min failure rate per component
    float failExp;  // exponent to scale failure rate
};

using Param = struct Param;

using schar = signed char;
using uchar = unsigned char;
// using uint  = unsigned int;
using ulong  = unsigned long;

using Allele    = float;
using AllelePtr = Allele*;
using FLOAT     = double;            // use double here to change to double precision
using FLOATPTR  = FLOAT*;
using FLTMATRIX = FLOAT**;
using GSL       = double;           // GSL uses double, so for all stats use this type
using GSLPTR    = GSL*;
using GSLMATRIX = GSL**;

using IntMatrix = int**;

#endif
