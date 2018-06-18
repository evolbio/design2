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

void SetBabyGenotype(Individual& Parent1, Individual &Parent2, Individual& baby)
{
	auto& g1 = Parent1.genotype;
	auto& g2 = Parent2.genotype;
	auto& gb = baby.genotype;
	ulong chrFlag = rnd.rbit();		// determines which parent is used for copying
	
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

