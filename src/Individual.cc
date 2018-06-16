#include <cmath>
#include "Individual.h"

// declare static private member variables

FLOAT Individual::mut;
FLOAT Individual::rec;
FLOAT Individual::s;
FLOAT Individual::minFail;
FLOAT Individual::failExp;
int Individual::totalLoci;

Individual::~Individual()
{
    delete [] genotype;
}

void Individual::initialize(Param& param)
{
    totalLoci = param.loci;
    genotype = new Allele[totalLoci];
    mut = param.mutation * totalLoci;

    // initialize values
    
    // randomly, set 95% of alleles to min failure probability, and 5% randomly on [minFail,1]
    for (int i = 0; i < totalLoci; ++i){
        genotype[i] = static_cast<Allele>((rnd.rU01() < 0.95) ? param.minFail : param.minFail + (1-param.minFail)*rnd.rU01());
    }
    fitness = calcFitness();
}

// set static variables used by class

void Individual::setParam(Param& param)
{
    mut = param.mutation;
    rec = param.recombination;
    s = param.s;
    totalLoci = param.loci;
    minFail = param.minFail;
    failExp = param.failExp;
}

// mut is per genotype mutation rate

void Individual::mutate()
{
    ulong hits = rnd.poisson(mut);
    for (unsigned i = 0; i < hits; ++i){
        ulong locus = rnd.rtop(totalLoci);
        genotype[locus] = static_cast<Allele>(minFail + (1.0-minFail)*rnd.rU01());  // uniform on [minFail, 1]
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

// s is selective coefficient for failure
// failure at each locus is m + (1-m)*((g-m)/(1-m))^failExp
// where m is minFail and g is genotypic value on range [m,1]
// This gives diminishing reduction in failure rate as genotypic value declines if failExp > 1, and 
// accelerating reduction in failure rate as genotypic value declines if failExp < 1.
// Diminishing return tends to push all loci to intermediate value for failure rate, since further declines
// would not give great benefits, but further increases cause large losses, like internal equilibrium
// under diminishing returns

FLOAT Individual::calcFitness()
{
    FLOAT probFailure = 1.0;
    AllelePtr g = genotype;
    for (int i = 0; i < totalLoci; ++i){
        probFailure *= minFail + (1-minFail)*pow((*g++-minFail)/(1-minFail),failExp);
    }
    return fitness = 1.0 - s*probFailure;
}

