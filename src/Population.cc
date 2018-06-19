#include "Population.h"
#include "gsl/gsl_statistics_double.h"
#include "gsl/gsl_spline.h"
extern "C"{
#include "util.h"
};

int dcompare(const void *e1, const void *e2)
{
    return ((*(const double *)e1) > (*(const double *)e2)) ? 1 :
              (((*(const double *)e1) < (*(const double *)e2)) ? -1 : 0);
}

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

// use type double here, which is type double, because double stat routines use double

void Population::calcStats(Param& param, SumStat& stats)
{
    int i, j;
    double d = param.distnSteps - 1.0;
    std::vector<double> fitness(param.popsize);
    DBLMATRIX gMatrix = new DBLPTR[param.loci];
    for (i = 0; i < param.loci; ++i){
        gMatrix[i] = new double[param.popsize];
    }
    
    for (i = 0; i < popSize; ++i){
        fitness[i] = cumFit[i] - ((i>0) ? cumFit[i-1] : 0.0);
        auto& genotype = ind[i].getGenotype();
        for (j = 0; j < param.loci; ++j){
            gMatrix[j][i] = genotype[j];
        }
    }
    DBLPTR gMean = stats.getGMean();
    DBLPTR gSD = stats.getGSD();
    DBLMATRIX gDistn = stats.getGDistn();
    DBLMATRIX gCorr = stats.getGCorr();
    for (i = 0; i < param.loci; ++i){
        gMean[i] = gsl_stats_mean(gMatrix[i], 1, param.popsize);
        gSD[i] = gsl_stats_sd(gMatrix[i], 1, param.popsize);
    }
    for (i = 0; i < param.loci; ++i){
        for (j = i; j < param.loci; ++j){
            double cov = gsl_stats_covariance(gMatrix[i], 1, gMatrix[j], 1, param.popsize);
            double prodSD = gSD[i]*gSD[j];
            gCorr[i][j] = gCorr[j][i] = (prodSD < 1e-10) ? 0.0 : cov/(prodSD);
        }
    }
    for (i = 0; i < param.loci; ++i){
        qsort(gMatrix[i], param.popsize, sizeof(double), dcompare);
        for (j = 0; j < param.distnSteps; ++j){
            gDistn[i][j] = gsl_stats_quantile_from_sorted_data(gMatrix[i], 1, param.popsize, (double)j/d);
        }
    }
    
    // fitness distn
    
    std::sort(fitness.begin(),fitness.end());
    stats.setAveFitness(gsl_stats_mean(fitness.data(), 1, param.popsize));
    stats.setSDFitness(gsl_stats_sd(fitness.data(), 1, param.popsize));
    DBLPTR fitnessDistn = stats.getFitnessDistn();
    for (j = 0; j < param.distnSteps; ++j){
        fitnessDistn[j] = gsl_stats_quantile_from_sorted_data(fitness.data(),
                                                                  1, param.popsize, (double)j/d);
    }
        
    for (i = 0; i < param.loci; ++i){
        delete [] gMatrix[i];
    }
    delete [] gMatrix;
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
