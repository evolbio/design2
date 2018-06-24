#include <cmath>
#include "Individual.h"
#include "Performance.h"

// declare static private member variables

double Individual::mut;
double Individual::rec;
int Individual::totalLoci;
Allele Individual::maxAllele;
double Individual::fitVar;
ulong Individual::negLog2Rec;

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
        genotype[i] = mutUniform();
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
    negLog2Rec = 1;         // set elsewhere when needed, here is just default value
}

Allele Individual::mutUniform()
{
    return static_cast<Allele>(rnd.rUniform(-maxAllele,maxAllele));
}

// could use bit cache for random bits to speed up

Allele Individual::mutStep(Allele a)
{
    const Allele d = (Allele)1.3;
    Allele c = (rnd.rbit()) ? -d : d;
    return a + c;
}

// mut is per genotype mutation rate

void Individual::mutate()
{
    // ulong hits = rnd.poisson(mut*totalLoci);  // about twice as slow, note change in type for hits
    int hits = MyRandomPoisson(mut*totalLoci);
    for (int i = 0; i < hits; ++i){
        ulong locus = rnd.rtop(totalLoci);
        genotype[locus] = mutStep(genotype[locus]);
        //genotype[locus] = mutUniform();
    }
}

// Possible modification: For example, if rec = 0.4, then use 10/25 to test, ie, use random numbers on [0..24] and success if random number is <= 9. Not sure how much performance boost would be obtained. Would only need five bits to generate random number, so could use bit operations on 64 bit random number.

void SetBabyGenotype(Individual& Parent1, Individual& Parent2, Individual& baby)
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

// This routine applies when -log2(rec) is integer 0,1,2,...
// rec = 1 => -log2 rec = 0 is OK here, each successive locus chosen from alternate parent
// assumes that random integer has random bits

void SetBabyGenotypeLogRec(Individual& Parent1, Individual& Parent2, Individual& baby)
{
    auto& g1 = Parent1.genotype;
    auto& g2 = Parent2.genotype;
    auto& gb = baby.genotype;
    ulong rawint = rnd.rawint();
    ulong recShift = Individual::negLog2Rec; // -log 2 rec, w/rec = (1/2, 1/4, 1/8, ...), set in Popul
    ulong mask = (1 << recShift) - 1;        // e.g., recShift = 2 => mask = 00...0011, ie, low two bits
    ulong chrFlag = rawint & 1;              // determines initial parent w/prob = 1/2, ie, random bit
    auto rbits = rnd.bitSize() - recShift;   // remaining bits available
    
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

// No recombination, choose just one parent and copy genotype to baby. Maybe memcpy would be faster?? Note how to turn off warning for unused parameter

void SetBabyGenotypeNoRec(Individual& Parent, Individual& Unused __attribute__((unused)), Individual& baby)
{
    for (int i = 0; i < Parent.totalLoci; ++i){
        baby.genotype[i] = Parent.genotype[i];
    }
    baby.fitness = baby.calcFitness();
}

// Calculation of num and den take from openVclose.h in pagmo optimization code; assumes dentilde = den, ie, not studying role of variable plant w/regard to stability margin.

double Individual::calcFitness()
{
    double tmax = 20.0;
    double gamma = 2.0;
    double a = sqrt(1+gamma);
    double p0 = genotype[0];
    double p1 = 1.0;
    double p2 = genotype[1];
    double q0 = genotype[2];
    double q1 = genotype[3];
    double q2 = genotype[4];
    std::vector<double> num {q2,q1,q0};
    std::vector<double> den {p2, p1+a*p2, p0+ a*p1 + p2, a*p0+p1, p0};

    double perf = performance(num, den, den, gamma, tmax, signalType::output);
    // update to center at optimum performance
    return fitness = exp(-perf*perf/(2*fitVar));
}

