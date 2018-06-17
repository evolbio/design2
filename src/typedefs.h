
#ifndef _typedefs_h
#define _typedefs_h 1

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

struct Param{
    int   runNum;
    int   distnSteps;
    int   gen;
    int   loci;
    int   popsize;
    float mutation;
    float recombination;
    Allele maxAllele; // max allelic value
    FLOAT fitVar;   // width of fitness gradient
};

using Param = struct Param;


#endif
