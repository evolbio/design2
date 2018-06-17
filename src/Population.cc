#include "Population.h"
#include "gsl/gsl_statistics_double.h"
#include "gsl/gsl_spline.h"
extern "C"{
#include "util.h"
};

int dcompare(const void *e1, const void *e2)
{
    return ((*(const GSL *)e1) > (*(const GSL *)e2)) ? 1 :
              (((*(const GSL *)e1) < (*(const GSL *)e2)) ? -1 : 0);
}

Population::Population(Param& param)
{
	popSize = param.popsize;
    totalLoci = param.loci;
	ind = new Individual[popSize];
    cumFit = new FLOAT[popSize];
    ind[0].setParam(param);
	for (int i = 0; i < popSize; i++){
		ind[i].initialize();
	}
}

Population::~Population()
{
	delete [] ind;
    delete [] cumFit;
}

void Population::setFitnessArray()
{
    Individual *iPtr = ind;
    FLOAT *cumPtr = cumFit;
    FLOAT fit = 0;
    for (int i = 0; i < popSize; ++i, ++iPtr, ++cumPtr){
        *cumPtr = fit += iPtr->getFitness();
    }
}

// Do everything on population in one loop

void Population::reproduceMutateCalcFit(Population& oldPop)
{
    Individual *iPtr = ind;
    FLOAT *cumPtr = cumFit;
    FLOAT fit = 0;
	for (int i = 0; i < popSize; ++i, ++iPtr, ++cumPtr){
        SetBabyGenotype(oldPop.chooseInd(), oldPop.chooseInd(), *iPtr);
        iPtr->mutate();
        *cumPtr = fit += iPtr->getFitness();
	}
}

// use type GSL here, which is type double, because GSL stat routines use double

void Population::calcStats(Param& param, SumStat& stats)
{
    int i, j;
    GSL d = param.distnSteps - 1.0;
    Individual *indPtr = ind;
    GSLPTR fitness = new GSL[popSize];
    GSLMATRIX gMatrix = new GSLPTR[param.loci];
    for (i = 0; i < param.loci; ++i){
        gMatrix[i] = new GSL[param.popsize];
    }
    
    for (i = 0; i < popSize; ++i, ++indPtr){
        fitness[i] = cumFit[i] - ((i>0) ? cumFit[i-1] : 0.0);
        const AllelePtr genotype = indPtr->getGenotype();
        for (j = 0; j < param.loci; ++j){
            gMatrix[j][i] = genotype[j];
        }
    }
    GSLPTR gMean = stats.getGMean();
    GSLPTR gSD = stats.getGSD();
    GSLMATRIX gDistn = stats.getGDistn();
    GSLMATRIX gCorr = stats.getGCorr();
    for (i = 0; i < param.loci; ++i){
        gMean[i] = gsl_stats_mean(gMatrix[i], 1, param.popsize);
        gSD[i] = gsl_stats_sd(gMatrix[i], 1, param.popsize);
    }
    for (i = 0; i < param.loci; ++i){
        for (j = i; j < param.loci; ++j){
            GSL cov = gsl_stats_covariance(gMatrix[i], 1, gMatrix[j], 1, param.popsize);
            GSL prodSD = gSD[i]*gSD[j];
            gCorr[i][j] = gCorr[j][i] = (prodSD < 1e-10) ? 0.0 : cov/(prodSD);
        }
    }
    for (i = 0; i < param.loci; ++i){
        qsort(gMatrix[i], param.popsize, sizeof(GSL), dcompare);
        for (j = 0; j < param.distnSteps; ++j){
            gDistn[i][j] = gsl_stats_quantile_from_sorted_data(gMatrix[i], 1, param.popsize, (double)j/d);
        }
    }
    
    // fitness distn
    
    qsort(fitness, param.popsize, sizeof(GSL), dcompare);
    stats.setAveFitness(gsl_stats_mean(fitness, 1, param.popsize));
    stats.setSDFitness(gsl_stats_sd(fitness, 1, param.popsize));
    GSLPTR fitnessDistn = stats.getFitnessDistn();
    for (j = 0; j < param.distnSteps; ++j){
        fitnessDistn[j] = gsl_stats_quantile_from_sorted_data(fitness,
                                                                  1, param.popsize, (double)j/d);
    }
        
    for (i = 0; i < param.loci; ++i){
        delete [] gMatrix[i];
    }
    delete [] gMatrix;
    delete [] fitness;
}


// array is cumulative fitness or weighting, n is number of elements
// returns random index from array sampled by cumulative success weighting

int Population::chooseMember(FLOAT *array, int n)
{
    FLOAT ran = rnd.rU01();
    FLOAT cumFitness = array[n-1];
    if (cumFitness < 1e-8)
        return static_cast<int>(rnd.rtop(n));

    FLOAT target = ran * cumFitness;
    int k = static_cast<int>((ran * (n-1)));

    while (array[k++] < target);
    while ((--k != -1) && (array[k] >= target));
    return ++k;
}
