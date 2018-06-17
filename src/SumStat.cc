#include "SumStat.h"


void SumStat::initialize(Param& param)
{
    int i;
    loci = param.loci;
    distnSteps = param.distnSteps;
    gMean = new GSL[loci];
    gSD = new GSL[loci];
    gDistn = new GSLPTR[loci];
    gCorr = new GSLPTR[loci];
    for (i = 0; i < loci; ++i){
        gDistn[i] = new GSL[param.distnSteps];		// steps for percentiles
        gCorr[i] = new GSL[loci];
    }
    fitnessDistn = new GSL[param.distnSteps];
}

SumStat::~SumStat()
{
    int i;
    delete [] gMean;
    delete [] gSD;
    for (i = 0; i < loci; ++i){
        delete [] gDistn[i];
        delete [] gCorr[i];
    }
    delete [] gDistn;
    delete [] gCorr;
    delete [] fitnessDistn;
}
