// Collect data from simulations here for later analysis and printing

#ifndef _SumStat_h
#define _SumStat_h 1

#include "sensitivity.h"
#include "typedefs.h"

class SumStat
{
public:
    ~SumStat();
    void		initialize(Param& param);
    auto&       getGMean(){return gMean;}
    auto&      getGSD(){return gSD;}
    DBLMATRIX	getGDistn(){return gDistn;}
    DBLMATRIX	getGCorr(){return gCorr;}
    auto&       getFitnessDistn(){return fitnessDistn;}
    double      getAveFitness(){return aveFitness;}
    double      getSDFitness(){return sdFitness;}
    void        setAveFitness(double x){aveFitness = x;}
    void        setSDFitness(double x){sdFitness = x;}
private:
    std::vector<double> gMean;	// mean values of component failure rates
    std::vector<double> gSD;	// sd values of component failure rates
    DBLMATRIX	gDistn;			// percentiles [0..100]  for component failure rates
    DBLMATRIX	gCorr;			// phenotypic correlation matrix for component failure rates
    int			loci;
    int         distnSteps;     // steps in percentile distns
    double      aveFitness;
    double      sdFitness;
    std::vector<double> fitnessDistn;
};

#endif
