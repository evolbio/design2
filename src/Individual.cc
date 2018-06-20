#include <cmath>
#include "Individual.h"

// declare static private member variables

double Individual::mut;
double Individual::rec;
int Individual::totalLoci;
Allele Individual::maxAllele;
double Individual::fitVar;

// Algorithm for fast Poisson for lambda < 30
// from https://www.johndcook.com/blog/2010/06/14/generating-poisson-random-values/
// Test measurement suggests about twice as fast as rnd.poisson()

int MyRandomPoisson(double mean)
{
    int k = 0;
    double p = 1.0;
    double L = exp(-mean);
    while (p > L){
        k++;
        p *= rnd.rU01();
    }
    return k-1;
}

// copy constructor
Individual::Individual(const Individual& other)
{
    genotype = std::unique_ptr<Allele[]> {new Allele[totalLoci]};
    for (int i = 0; i < totalLoci; ++i){
        genotype[i] = other.genotype[i];
    }
    fitness = calcFitness();
}

// assignment constructor
Individual& Individual::operator=(const Individual& other)
{
    genotype = std::unique_ptr<Allele[]> {new Allele[totalLoci]};
    for (int i = 0; i < totalLoci; ++i){
        genotype[i] = other.genotype[i];
    }
    fitness = calcFitness();
    return *this;
}

void Individual::initialize()
{
    genotype = std::unique_ptr<Allele[]> {new Allele[totalLoci]};
    
    for (int i = 0; i < totalLoci; ++i){
        genotype[i] = static_cast<Allele>(rnd.rUniform(-maxAllele,maxAllele));
    }
    fitness = calcFitness();
}

// set static variables used by class

void Individual::setParam(Param& param)
{
    mut = param.mutation;
    rec = param.recombination;
    totalLoci = param.loci;
    maxAllele = param.maxAllele;
    fitVar = param.fitVar;
}

// mut is per genotype mutation rate

void Individual::mutate()
{
    // ulong hits = rnd.poisson(mut*totalLoci);  // about twice as slow, note change in type for hits
    int hits = MyRandomPoisson(mut*totalLoci);
    for (int i = 0; i < hits; ++i){
        ulong locus = rnd.rtop(totalLoci);
        genotype[locus] = static_cast<Allele>(rnd.rUniform(-maxAllele,maxAllele));
    }
}

// No recombination, choose just one parent and copy genotype to baby
// Maybe memcpy would be faster??

void SetBabyGenotype(Individual& Parent, Individual& baby)
{
    for (int i = 0; i < Parent.totalLoci; ++i){
        baby.genotype[i] = Parent.genotype[i];
    }
    baby.fitness = baby.calcFitness();
}

// This routine applies when -log2(rec) is integer 0,1,2,...
// rec = 1 => -log2 rec = 0 is OK here, each successive locus chosen from alternate parent
// assumes that random integer has random bits

void SetBabyGenotypeLog(Individual& Parent1, Individual &Parent2, Individual& baby)
{
    auto& g1 = Parent1.genotype;
    auto& g2 = Parent2.genotype;
    auto& gb = baby.genotype;
    ulong rawint = rnd.rawint();
    ulong recShift = 2;                     // -log 2 recombination, w/ rec = (1/2, 1/4, 1/8, ...)
    ulong mask = (1 << recShift) - 1;       // e.g., recShift = 2 => mask = 00...0011, ie, low two bits
    ulong chrFlag = rawint & 1;             // determines initial parent w/prob = 1/2, ie, random bit
    auto rbits = rnd.bitSize() - recShift;  // remaining bits available
    
    for (int i = 0; i < Parent1.totalLoci; ++i){
        gb[i] = (chrFlag) ? g1[i] : g2[i];
        rawint >>= recShift;                            // move used bits out
        if ((rawint & mask) == mask) chrFlag ^= 1;      // flip flag if recombination
        if ((rbits -= recShift) == 0){                  // reload random bits if all used up
            rawint = rnd.rawint();                      // new random int
            rbits = rnd.bitSize();                      // reset remaining bits left to use
        }
    }
    baby.fitness = baby.calcFitness();
}

void SetBabyGenotype(Individual& Parent1, Individual &Parent2, Individual& baby)
{
    auto& g1 = Parent1.genotype;
    auto& g2 = Parent2.genotype;
    auto& gb = baby.genotype;
    ulong chrFlag = rnd.rbit();        // determines which parent is used for copying

    for (int i = 0; i < Parent1.totalLoci; ++i){
        gb[i] = (chrFlag) ? g1[i] : g2[i];
        if (rnd.rU01() < baby.rec) chrFlag ^= 1;        // flip flag if recombination at rate 0.5
    }
    baby.fitness = baby.calcFitness();
}

double Individual::calcFitness()
{
    double sumAllele = 0;
    for (int i = 0; i < totalLoci; ++i){
        sumAllele += genotype[i];
    }
    return fitness = exp(-sumAllele*sumAllele/(2*fitVar));
}

