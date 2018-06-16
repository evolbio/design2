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
    GSLPTR      getGMean(){return gMean;}
    GSLPTR      getGSD(){return gSD;}
    GSLMATRIX	getGDistn(){return gDistn;}
    GSLMATRIX	getGCorr(){return gCorr;}
    GSL         getAveDisease(){return aveDisease;}
    GSL         getSDDisease(){return sdDisease;}
    GSLPTR      getDiseaseDistn(){return diseaseDistn;}
    GSLPTR      getDiseaseDistnNormal(){return diseaseDistnNormal;}
    GSLPTR      getFitnessDistn(){return fitnessDistn;}
    GSL         getAveFitness(){return aveFitness;}
    GSL         getSDFitness(){return sdFitness;}
    void        setAveFitness(GSL x){aveFitness = x;}
    void        setSDFitness(GSL x){sdFitness = x;}
    void        setAveDisease(GSL x){aveDisease = x;}
    void        setSDDisease(GSL x){sdDisease = x;}
private:
    GSLPTR      gMean;			// mean values of component failure rates
    GSLPTR      gSD;			// sd values of component failure rates
    GSLMATRIX	gDistn;			// percentiles [0..100]  for component failure rates
    GSLMATRIX	gCorr;			// phenotypic correlation matrix for component failure rates
    int			loci;
    int         distnSteps;     // steps in percentile distns
    GSL         aveDisease;       // freq of disease
    GSL         sdDisease;
    GSLPTR      diseaseDistn;     // percentiles [0..100]
    GSLPTR      diseaseDistnNormal;     // normalized, see calcStats() in Population.cc
    GSL         aveFitness;
    GSL         sdFitness;
    GSLPTR      fitnessDistn;
};

#endif
