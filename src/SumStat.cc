#include "SumStat.h"


void SumStat::initialize(Param& param)
{
    int i;
    loci = param.loci;
    distnSteps = param.distnSteps;
    gMean = std::vector<double>(loci);
    gSD = std::vector<double>(loci);
    gDistn = new DBLPTR[loci];
    gCorr = new DBLPTR[loci];
    for (i = 0; i < loci; ++i){
        gDistn[i] = new double[param.distnSteps];		// steps for percentiles
        gCorr[i] = new double[loci];
    }
    fitnessDistn = std::vector<double>(param.distnSteps);
}

SumStat::~SumStat()
{
    int i;
    for (i = 0; i < loci; ++i){
        delete [] gDistn[i];
        delete [] gCorr[i];
    }
    delete [] gDistn;
    delete [] gCorr;
}
