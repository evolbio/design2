#include "Population.h"
#include "gsl/gsl_statistics_double.h"
#include "util.h"

Population::Population(Param& param)
{
	popSize = param.popsize;
    ind = std::vector<Individual>(popSize);
    cumFit = std::vector<double>(popSize);
    ind[0].setParam(param);
	for (int i = 0; i < popSize; i++){
		ind[i].initialize();
	}
}

void Population::setFitnessArray()
{
    double fit = 0;
    for (int i = 0; i < popSize; ++i){
        cumFit[i] = fit += ind[i].getFitness();
    }
}

// Do everything on population in one loop

void Population::reproduceMutateCalcFit(Population& oldPop)
{
    double fit = 0;
	for (int i = 0; i < popSize; ++i){
        SetBabyGenotype(oldPop.chooseInd(), oldPop.chooseInd(), ind[i]);
        ind[i].mutate();
        cumFit[i] = fit += ind[i].getFitness();
	}
}

// gsl takes double C arrays
void Population::calcStats(Param& param, SumStat& stats)
{
    int i, j;
    std::vector<double> fitness(param.popsize);
    std::vector<std::vector<double>> gMatrix(param.loci,std::vector<double>(param.popsize));
    
    for (i = 0; i < popSize; ++i){
        fitness[i] = cumFit[i] - ((i>0) ? cumFit[i-1] : 0.0);
        auto& genotype = ind[i].getGenotype();
        for (j = 0; j < param.loci; ++j){
            gMatrix[j][i] = genotype[j];
        }
    }
    auto& gMean = stats.getGMean();
    auto& gSD = stats.getGSD();
    auto& gCorr = stats.getGCorr();
    for (i = 0; i < param.loci; ++i){
        gMean[i] = vecMean<double>(gMatrix[i]);
        gSD[i] = vecSD<double>(gMatrix[i], gMean[i]);
    }
    for (i = 0; i < param.loci; ++i){
        for (j = i; j < param.loci; ++j){
            double cov = vecCov(gMatrix[i], gMatrix[j], gMean[i], gMean[j]);
            double prodSD = gSD[i]*gSD[j];
            gCorr[i][j] = gCorr[j][i] = (prodSD < 1e-10) ? 0.0 : cov/(prodSD);
        }
    }
    
    // distn of allelic values
    
    std::vector<unsigned> ptiles(param.distnSteps);
    std::iota(ptiles.begin(), ptiles.end(), 0);     // assign [0..n-1] for distnSteps = n, use n = 101
    auto& gDistn = stats.getGDistn();

    for (i = 0; i < param.loci; ++i){
        gDistn[i] = percentiles_interpol<std::vector<double>>(gMatrix[i], ptiles);
    }
    
    // fitness distn
    
    double mean = vecMean<double>(fitness);
    stats.setAveFitness(mean);
    stats.setSDFitness(vecSD<double>(fitness, mean));
    
    auto& fitnessDistn = stats.getFitnessDistn();
    fitnessDistn = percentiles_interpol<std::vector<double>>(fitness, ptiles);
}

// array is cumulative fitness or weighting, n is number of elements
// returns random index from array sampled by cumulative success weighting

int Population::chooseMember(double *array, int n)
{
    double ran = rnd.rU01();
    double cumFitness = array[n-1];
    if (cumFitness < 1e-8)
        return static_cast<int>(rnd.rtop(n));

    double target = ran * cumFitness;
    int k = static_cast<int>((ran * (n-1)));

    while (array[k++] < target);
    while ((--k != -1) && (array[k] >= target));
    return ++k;
}
