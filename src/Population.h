#ifndef _Population_h
#define _Population_h 1

#include "sensitivity.h"
#include "typedefs.h"
#include "Individual.h"
#include "SumStat.h"

// Life cycle is make a baby, mutate the baby, calculate its fitness,
// analyze the population characteristics every so often, reproduce

class Population
{
public:
	Population(Param& param);
	~Population();
	int			getPopSize(){return popSize;}
	Individual&	getInd(int i){return ind[i];}
    Individual& chooseInd(){return ind[chooseMember(cumFit,popSize)];}
    void		setFitnessArray();
	void		reproduceMutateCalcFit(Population& oldPop);
    void		calcStats(Param& param, SumStat& stats);
private:
    int     	chooseMember(double *array, int n);
	int 		popSize;		// # females = # males = popSize
    int			totalLoci;      // must store for destructor
	Individual*	ind;			// array of individuals
    double*		cumFit;         // cumfitness of individuals
};

#endif
