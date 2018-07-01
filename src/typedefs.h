
#ifndef _typedefs_h
#define _typedefs_h 1

#include APPL_H

using schar = signed char;
using uchar = unsigned char;
// using uint  = unsigned int;
using ulong  = unsigned long;

using Allele    = float;
using AllelePtr = Allele*;

enum class Loop {open, close, dclose};

struct Param {
    int   runNum;
    int   distnSteps;
    int   gen;
    int   loci;
    int   popsize;
    int   mutLocus;
    unsigned long rndSeed;  // seed is 64bit, must type cast if rndType is 32bit
    double mutation;
    double recombination;
    double stochWt;        // weighting of stochastic fluctuations
    bool   stoch;          // (stochWt == 0) ? false : true
    Allele mutStep;     // max allelic value
    double aSD;
    double fitVar;      // width of fitness gradient
    double gamma;
    Loop loop;
    std::string rec;
};

#endif
