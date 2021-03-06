#ifndef _Individual_h
#define _Individual_h 1

#include APPL_H
#include "typedefs.h"
#include "Individual.h"

// Use array of floats for genotype.

// Each individual is haploid hermaphrodite, so no distinct sexes
// Form diploid zygote and then make a gamete to produce haploid baby, ie, haploid dominant life cycle
// Thus, no dominance, dominance is favorable to maintenance of variability, so assumptions unfavorable for variability and therefore isolates effect of interlocus interactions over robustness

// Make global parameters as static variables for class, so can access them without always passing Param
// Static variables must be declared in Individual.cc, and values initialized by main program

// stochast is array of phenotypic stochasticity Alleles, with one-to-one map of stochasticity to genotype alleles. For recombination, stochast alleles linked to genotype alleles, ie, no recombination between each genotype and its associated stochasticity allele. Value of stochast is standard deviation of Gaussian fluctuations, each fluctuation weighted by param stochWt. If stochWt == 0, then param stoch = false and ignore.

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
    static void     setParam(Param& param);             // set static variables for class
    void            setNegLog2Rec(ulong r) {negLog2Rec = r;};
    void			initialize();
    void			mutate();
    void            mutateG(std::unique_ptr<Allele []>&, bool);
    auto            getRecombination(){return rec;}
    void            setRecombination(double r){rec = r;}
    double          calcJ();
    double			calcFitness();
    double          getFitness(){return fitness;};
    auto&           getGenotype(){return genotype;};
    auto&           getStochast(){return stochast;};
    Allele          mutateStep(Allele a);
private:
    static double	mut;            // per genome mutation rate, param.mutation is per locus mutation rate
    static int		totalLoci;
    static double   rec;            // recombination probability
    static ulong    negLog2Rec;     // -log2 recombination, used when rec = 1, 1/2, 1/4, ...
    static Allele   mutStep;        // size of mutational step
    static double   aSD;            // variability of plant parameter a
    static double   fitVar;         // variance of fitness scaling
    static double   gamma;          // weighting of performance components
    static Loop     loop;           // control loop type
    static int      mutLocus;       // if >= 0, then mutate only this locus
    static double   stochWt;        // weighting of stochastic fluctuations
    static bool     stoch;          // (stochWt == 0) ? false : true
    std::unique_ptr<Allele[]> genotype;
    std::unique_ptr<Allele[]> stochast;  // phenotypic stochasticity
    double          fitness;
};

#endif
