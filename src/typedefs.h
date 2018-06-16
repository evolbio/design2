
#ifndef _typedefs_h
#define _typedefs_h 1

struct Param{
    int runNum;
    int time;               // end time for run
    int tinit;              // initial time before collecting stats
    char bndry;             // p => periodic; r => reflective
    double step;            // length of time step
    int oincr;              // number of steps between data collection
    int sinit;              // initial # of S molecules
    int xinit;              // initial # of X molecules
    int rinit;              // initial # of R molecules
    double diff;            // base diffusion rate, added starting w/ExpF
    double ds,dx,da,dr;     // diffusion coefficients, each multiplied by diff
    double rsxf,rsxb;       // forward and backward rates S+X->SX
    double rapr;            // A prod rate, SX -> SX + A
    double radc;            // A decay rate, A -> 0
    double gmsx;            // geminate recombination prob for SX back rx
    char *pfile;            // smoldyn param file in /tmp/ directory
    char *ofile;            // smoldyn output file in /tmp/ directory
};

typedef struct Param Param;

typedef signed char schar;
typedef unsigned char uchar;
typedef unsigned int uint;

#endif
