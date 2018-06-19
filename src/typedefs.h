
#ifndef _typedefs_h
#define _typedefs_h 1

using schar = signed char;
using uchar = unsigned char;
// using uint  = unsigned int;
using ulong  = unsigned long;

using Allele    = float;
using AllelePtr = Allele*;

struct Param{
    int   runNum;
    int   distnSteps;
    int   gen;
    int   loci;
    int   popsize;
    float mutation;
    float recombination;
    Allele maxAllele; // max allelic value
    double fitVar;   // width of fitness gradient
};

using Param = struct Param;


#endif
