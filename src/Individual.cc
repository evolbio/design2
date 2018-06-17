#include <cmath>
#include "Individual.h"

// declare static private member variables

FLOAT Individual::mut;
FLOAT Individual::rec;
int Individual::totalLoci;
Allele Individual::maxAllele;
FLOAT Individual::fitVar;


Individual::~Individual()
{
    delete [] genotype;
}

void Individual::initialize()
{
    genotype = new Allele[totalLoci];

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
    ulong hits = rnd.poisson(mut);
    for (unsigned i = 0; i < hits; ++i){
        ulong locus = rnd.rtop(totalLoci);
        genotype[locus] = static_cast<Allele>(rnd.rUniform(-maxAllele,maxAllele));
    }
}

void SetBabyGenotype(Individual& Parent1, Individual &Parent2, Individual& baby)
{
	AllelePtr g1 = Parent1.genotype;
	AllelePtr g2 = Parent2.genotype;
	AllelePtr gb = baby.genotype;
	ulong chrFlag = rnd.rbit();		// determines which parent is used for copying
    FLOAT probFailure = 1.0;        // fitness is prob of not failing, which is 1 - product of failure prob
                                    // at each locus, calc here for efficiency
	
	for (int i = 0; i < Parent1.totalLoci; i++, g1++, g2++, gb++){
		probFailure *= *gb = (chrFlag) ? *g1 : *g2;
		if (rnd.rU01() < baby.rec) chrFlag ^= 1;		// flip flag if recombination at rate 0.5
	}
    baby.fitness = baby.calcFitness();
}

FLOAT Individual::calcFitness()
{
    FLOAT sumAllele = 0;
    for (int i = 0; i < totalLoci; ++i){
        sumAllele += genotype[i];
    }
    return fitness = exp(-sumAllele*sumAllele/(2*fitVar));
}

