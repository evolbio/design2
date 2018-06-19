#include "SumStat.h"


void SumStat::initialize(Param& param)
{
    loci = param.loci;
    distnSteps = param.distnSteps;
    gMean = std::vector<double>(loci);
    gSD = std::vector<double>(loci);
    gDistn = std::vector<std::vector<double>>(loci);
    gCorr = std::vector<std::vector<double>>(loci, std::vector<double>(loci));
}
