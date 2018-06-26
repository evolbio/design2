
#ifndef _typedefs_h
#define _typedefs_h 1

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
    double mutation;
    double recombination;
    Allele mutStep;     // max allelic value
    double aSD;
    double fitVar;      // width of fitness gradient
    double gamma;
    Loop loop;
    std::string rec;
};

#endif
