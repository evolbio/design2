// Collect data from simulations here for later analysis and printing

#ifndef _SumStat_h
#define _SumStat_h 1

#include "sensitivity.h"
#include "typedefs.h"

class SumStat
{
public:
    void		initialize(Param& param);
    auto&       getGMean(){return gMean;}
    auto&       getGSD(){return gSD;}
    auto&       getGDistn(){return gDistn;}
    auto&       getGCorr(){return gCorr;}
    auto&       getSMean(){return sMean;}
    auto&       getSSD(){return sSD;}
    auto&       getSDistn(){return sDistn;}
    auto&       getSCorr(){return sCorr;}
    auto&       getSGCorr(){return sgCorr;}
    auto&       getFitnessDistn(){return fitnessDistn;}
    double      getAveFitness(){return aveFitness;}
    double      getSDFitness(){return sdFitness;}
    void        setAveFitness(double x){aveFitness = x;}
    void        setSDFitness(double x){sdFitness = x;}
    auto&       getPerfDistn(){return perfDistn;}
    double      getAvePerf(){return avePerf;}
    double      getSDPerf(){return sdPerf;}
    void        setAvePerf(double x){avePerf = x;}
    void        setSDPerf(double x){sdPerf = x;}
private:
    std::vector<double> gMean;                  // mean values of alleles
    std::vector<double> gSD;                    // sd values of alleles
    std::vector<std::vector<double>> gDistn;    // percentiles [0..100]  of alleles
    std::vector<std::vector<double>> gCorr;     // correlation matrix of alleles
    std::vector<double> sMean;                  // mean values of stochastic
    std::vector<double> sSD;                    // sd values of stochastic
    std::vector<std::vector<double>> sDistn;    // percentiles [0..100]  of stochastic
    std::vector<std::vector<double>> sCorr;     // correlation matrix of stochastic
    std::vector<std::vector<double>> sgCorr;    // cross correlation of genotype & stochastic
    int			loci;
    int         distnSteps;     // steps in percentile distns
    double      aveFitness;
    double      sdFitness;
    double      avePerf;
    double      sdPerf;
    std::vector<double> fitnessDistn;
    std::vector<double> perfDistn;
};

#endif
