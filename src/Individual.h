#ifndef _Individual_h
#define _Individual_h 1

#include "sensitivity.h"
#include "typedefs.h"
#include "Individual.h"

// Use float for fitness with about 6-7 places of precision because pop sizes will be less than 10^6 and so any fitness differences less than the pop size will be effectively neutral, if pop size gets to be on order of 10^6, should go to double precision.

// Use array of floats for genotype, each locus has a value on [0,1] which is probability of failure for that locus, phenotype is probability of failure across all loci, ie, product of probabilities across all loci; note optimum is zero.  If value at a locus is 0<=p<=1, then mutation causes new value that is worse than original, ie, p < p' <= 1, so must initialize pop with mostly 0 to provide best genotypes, and a small percentage of variable values to initialize with some variability, and let mutation cause only deterioration.

// Each individual is haploid hermaphrodite, so no distinct sexes
// Form diploid zygote and then make a gamete to produce haploid baby, ie, haploid dominant life cycle
// Thus, no dominance, dominance is favorable to maintenance of variability, so assumptions unfavorable for variability and therefore isolates effect of interlocus interactions over robustness

// Make global parameters as static variables for class, so can access them without always passing Param
// Static variables must be declared in Individual.cc, and values initialized by main program


class Individual
{
    friend void SetBabyGenotype(Individual& Parent1, Individual &Parent2, Individual& baby);
public:
    ~Individual();
    void			initialize(Param& param);
    void            setParam(Param& param);     // set static variables for class
    void			mutate();
    FLOAT			calcFitness();
    FLOAT           getFitness(){return fitness;};
    AllelePtr	    getGenotype(){return genotype;};
private:
    static FLOAT	mut;            // per genome mutation rate, param.mutation is per locus mutation rate
    static int		totalLoci;
    static FLOAT    s;              // selection coefficient
    static FLOAT    rec;            // recombination probability
    static FLOAT    minFail;        // minimum failure rate
    static FLOAT    failExp;        // failure exponent scaling
    AllelePtr 		genotype;       // array of failure probabilities, haploid so each allele is a single failure prob
    FLOAT           fitness;
};

#endif
