#include "SumStat.h"


void SumStat::initialize(Param& param)
{
    loci = param.loci;
    distnSteps = param.distnSteps;
    gMean = std::vector<double>(loci);
    gSD = std::vector<double>(loci);
    gDistn = std::vector<std::vector<double>>(loci);
    gCorr = std::vector<std::vector<double>>(loci, std::vector<double>(loci));
    if (param.stoch){
        sMean = std::vector<double>(loci);
        sSD = std::vector<double>(loci);
        sDistn = std::vector<std::vector<double>>(loci);
        sCorr = std::vector<std::vector<double>>(loci, std::vector<double>(loci));
        sgCorr = std::vector<std::vector<double>>(loci, std::vector<double>(loci));
    }
}
