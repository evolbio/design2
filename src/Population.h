#ifndef _Population_h
#define _Population_h 1

#include APPL_H
#include "typedefs.h"
#include "Individual.h"
#include "SumStat.h"

// Life cycle is make a baby, mutate the baby, calculate its fitness,
// analyze the population characteristics every so often, reproduce

class Population
{
public:
	Population(Param& param);
	int			getPopSize(){return popSize;}
	Individual&	getInd(int i){return ind[i];}
    Individual& chooseInd(){return ind[getRandIndex()];}
    void        partialSortInd(int sortToIndex);  // sort first percent of individuals by fitness
    void		setFitnessArray();
	void		reproduceMutateCalcFit(Population& oldPop);
    void        reproduceNoMutRec(Population& oldPop);
    void		calcStats(Param& param, SumStat& stats);
    void        createAliasTable();
private:
    int     	chooseMember(double *array, int n);
	int 		popSize;
    std::vector<Individual>	ind;			// vector of individuals
    std::vector<double>		indFitness;     // fitness of individuals
    std::vector<uint64_t>   hvec;
    std::vector<uint32_t>   avec;
    uint32_t getRandIndex();
    void (*SetBaby)(Individual&, Individual&, Individual&);
};

#endif
