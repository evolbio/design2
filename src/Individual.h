#ifndef _Individual_h
#define _Individual_h 1

#include "sensitivity.h"
#include "typedefs.h"
#include "Individual.h"

// Test performance for float vs double for fitness. Could use float for fitness with about 6-7 places of precision because pop sizes will be less than 10^6 and so any fitness differences less than the pop size will be effectively neutral, if pop size gets to be on order of 10^6, should go to double precision. However, performance for double may be same, and not much space taken up.

// Use array of floats for genotype, each locus has a value on [-MAX,MAX]. For initial test problem, set phenotype, P, as sum of allelic values across all loci, and fitness as exp(-P^2/2*S^2). Then can test standard mut-selection pattern for one locus.

// Each individual is haploid hermaphrodite, so no distinct sexes
// Form diploid zygote and then make a gamete to produce haploid baby, ie, haploid dominant life cycle
// Thus, no dominance, dominance is favorable to maintenance of variability, so assumptions unfavorable for variability and therefore isolates effect of interlocus interactions over robustness

// Make global parameters as static variables for class, so can access them without always passing Param
// Static variables must be declared in Individual.cc, and values initialized by main program

class Individual;

// must declare in general scope to use pointer to function later
// Log version when recombination is given as -log2 = 1,2,..., ie, as 1/2, 1/4, 1/8, ...
// No recombination version, just copy parent genotype to baby
// No recombination, have Unused parameter so all functions have same args

void SetBabyGenotype(Individual&, Individual&, Individual&);
void SetBabyGenotypeLogRec(Individual&, Individual&, Individual&);
void SetBabyGenotypeNoRec(Individual&, Individual& Unused, Individual&);

class Individual
{
    friend void SetBabyGenotype(Individual& Parent1, Individual& Parent2, Individual& baby);
    friend void SetBabyGenotypeLogRec(Individual& Parent1, Individual& Parent2, Individual& baby);
    friend void SetBabyGenotypeNoRec(Individual& Parent, Individual& Unused, Individual& baby);
public:
    Individual(){};
    Individual(const Individual& other);                // copy constructor
    Individual&     operator=(const Individual& other); // assignment constructor
    void            setParam(Param& param);     // set static variables for class
    void            setNegLog2Rec(ulong r) {negLog2Rec = r;};
    void			initialize();
    void			mutate();
    auto            getRecombination(){return rec;}
    void            setRecombination(double r){rec = r;}
    double          calcJ();
    double			calcFitness();
    double          getFitness(){return fitness;};
    auto&   	    getGenotype(){return genotype;};
    Allele          mutateStep(Allele a);
private:
    static double	mut;            // per genome mutation rate, param.mutation is per locus mutation rate
    static int		totalLoci;
    static double   rec;            // recombination probability
    static ulong    negLog2Rec;     // -log2 recombination, used when rec = 1, 1/2, 1/4, ...
    static Allele   mutStep;        // size of mutational step
    static double   fitVar;         // variance of fitness scaling
    static double   gamma;          // weighting of performance components
    static Loop     loop;           // control loop type
    std::unique_ptr<Allele[]> genotype;
    double          fitness;
};

#endif
