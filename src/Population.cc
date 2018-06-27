#include <cmath>

#include "Population.h"
#include "util.h"

Population::Population(Param& param)
{
    static bool flag = true;
    std::string showRec;
    popSize = param.popsize;
    ind = std::vector<Individual>(popSize);
    indFitness = std::vector<double>(popSize);
    hvec = std::vector<uint64_t>(popSize);
    avec = std::vector<uint32_t>(popSize);

    Individual::setParam(param);            // static function to set static variables
	for (int i = 0; i < popSize; i++){
		ind[i].initialize();
        indFitness[i] = ind[i].getFitness();
	}
    auto rec = param.recombination;
    double logRec = -log2(rec);
    ulong logRecRound = round<ulong>(logRec);
    if (rec < 1e-7){
        SetBaby = SetBabyGenotypeNoRec;
        showRec = "Rec: no recombination\n";
        param.rec = "None";
    }
    else if (abs(logRec - logRecRound) < 1e-2){
        SetBaby = SetBabyGenotypeLogRec;
        ind[0].setNegLog2Rec(logRecRound);
        showRec = fmt::format("Rec: using Log = {} -> {}\n", logRec,logRecRound);
        param.rec = fmt::format("Log {}", logRecRound);
    }
    else{
        SetBaby = SetBabyGenotype;
        showRec = fmt::format("Rec: using Uniform, Log = {}\n", logRec);
        param.rec = "Uniform";
    }
    if (flag && showProgress){
        flag=false;
        std::cout << showRec;
    }
}

void Population::setFitnessArray()
{
    for (int i = 0; i < popSize; ++i){
        indFitness[i] = ind[i].getFitness();
    }
}

// Do everything on population in one loop

void Population::reproduceMutateCalcFit(Population& oldPop)
{
    oldPop.createAliasTable();
    for (int i = 0; i < popSize; ++i){
        SetBaby(oldPop.chooseInd(), oldPop.chooseInd(), ind[i]);
        ind[i].mutate();
        indFitness[i] = ind[i].getFitness();
    }
}

// Round of reproduction without mutation or recombination, useful for testing models in which final population for stats is a population formed after selection but before recombination or mutation

void Population::reproduceNoMutRec(Population& oldPop)
{
    auto r = ind[0].getRecombination();
    oldPop.createAliasTable();
    for (int i = 0; i < popSize; ++i){
        SetBaby(oldPop.chooseInd(), oldPop.chooseInd(), ind[i]);
        indFitness[i] = ind[i].getFitness();
    }
    ind[0].setRecombination(r);
}

// If using stochastic loci for phenotypic variability, then simply double number of loci for allocation of vectors and matrix, and use 1..L for genotype and L+1,...,2L for stochastic alleles

void Population::calcStats(Param& param, SumStat& stats)
{
    int i, j;
    int loci = param.loci;
    int lociTot = (param.stoch) ? 2*loci : loci;
    std::vector<std::vector<double>> gMatrix(lociTot,std::vector<double>(param.popsize));
    
    for (i = 0; i < popSize; ++i){
        auto& genotype = ind[i].getGenotype();
        auto& stochast = ind[i].getStochast();
        for (j = 0; j < loci; ++j){
            gMatrix[j][i] = genotype[j];
            if (param.stoch) gMatrix[loci+j][i] = stochast[j];
        }
    }
    auto& gMean = stats.getGMean();
    auto& gSD = stats.getGSD();
    auto& gCorr = stats.getGCorr();
    auto& sMean = stats.getSMean();
    auto& sSD = stats.getSSD();
    auto& sCorr = stats.getSCorr();
    for (i = 0; i < loci; ++i){
        gMean[i] = vecMean<double>(gMatrix[i]);
        gSD[i] = vecSD<double>(gMatrix[i], gMean[i]);
        if (param.stoch){
            sMean[i] = vecMean<double>(gMatrix[loci+i]);
            sSD[i] = vecSD<double>(gMatrix[loci+i], sMean[i]);
        }
    }
    for (i = 0; i < loci; ++i){
        for (j = i; j < loci; ++j){
            // corr of g loci
            double cov = vecCov(gMatrix[i], gMatrix[j], gMean[i], gMean[j]);
            double prodSD = gSD[i]*gSD[j];
            gCorr[i][j] = gCorr[j][i] = (prodSD < 1e-10) ? 0.0 : cov/(prodSD);
            if (param.stoch){
                // corr of s loci
                double covS = vecCov(gMatrix[loci+i], gMatrix[loci+j], sMean[i], sMean[j]);
                double prodSDS = sSD[i]*sSD[j];
                sCorr[i][j] = sCorr[j][i] = (prodSDS < 1e-10) ? 0.0 : covS/(prodSDS);
                // cross corr of g[i] and s[j]
                double covSG = vecCov(gMatrix[i], gMatrix[loci+j], gMean[i], sMean[j]);
                double prodSDSG = gSD[i]*sSD[j];
                sCorr[i][j] = sCorr[j][i] = (prodSDSG < 1e-10) ? 0.0 : covSG/(prodSDSG);
            }
        }
    }
    
    // distn of allelic values
    
    std::vector<unsigned> ptiles(param.distnSteps);
    std::iota(ptiles.begin(), ptiles.end(), 0);     // assign [0..n-1] for distnSteps = n, use n = 101
    auto& gDistn = stats.getGDistn();
    auto& sDistn = stats.getSDistn();

    for (i = 0; i < loci; ++i){
        gDistn[i] = percentiles_interpol<std::vector<double>>(gMatrix[i], ptiles);
        if (param.stoch) sDistn[i] = percentiles_interpol<std::vector<double>>(gMatrix[loci+i], ptiles);
    }
    
    // fitness distn
    
    double mean = vecMean<double>(indFitness);
    stats.setAveFitness(mean);
    stats.setSDFitness(vecSD<double>(indFitness, mean));
    
    auto& fitnessDistn = stats.getFitnessDistn();
    fitnessDistn = percentiles_interpol<std::vector<double>>(indFitness, ptiles);
    
    // perf distn
    
    std::vector<double> indPerf(popSize);
    for (i = 0; i < popSize; ++i){
        indPerf[i] = ind[i].calcJ();
    }
    double pmean = vecMean<double>(indPerf);
    stats.setAvePerf(pmean);
    stats.setSDPerf(vecSD<double>(indPerf, pmean));
    
    auto& perfDistn = stats.getPerfDistn();
    perfDistn = percentiles_interpol<std::vector<double>>(indPerf, ptiles);
}

// Alias method for sampling from discrete distribution, see https://pandasthumb.org/archives/2012/08/lab-notes-the-a.html and https://en.wikipedia.org/wiki/Alias_method and http://www.keithschwarz.com/darts-dice-coins/
// Example code for testing in ~/sim/02_SmallTests/aliasMethod.cc

uint32_t Population::getRandIndex(){
    uint64_t u = rnd.rawint();
    uint32_t x = u % popSize;
    return (u < hvec[x]) ? x : avec[x];
}

void Population::createAliasTable() {
    const double *pp = indFitness.data();
    uint64_t *h = hvec.data();
    std::uint32_t *a = avec.data();
    uint32_t n = popSize;
    
    // normalize pp and copy into buffer
    double f = 0.0;
    std::vector<double> p(n);
    for (uint32_t i = 0; i < n; ++i)
        f += pp[i];
    f = static_cast<double>(n) / f;
    for (uint32_t i = 0; i < n; ++i)
        p[i] = pp[i] * f;
    
    // find starting positions, g => less than target, m greater than target
    uint32_t g, m, mm;
    for (g = 0; g < n && p[g] < 1.0; ++g)
    /*noop*/;
    for (m = 0; m < n && p[m] >= 1.0; ++m)
    /*noop*/;
    mm = m + 1;

    // build alias table until we run out of large or small bars
    while (g < n && m < n) {
        // convert double to 64-bit integer, control for precision
        // 9007... is max integer in double format w/53 bits, then shift by 11
        h[m] = (static_cast<uint64_t>(ceil(p[m] * 9007199254740992.0)) << 11);
        a[m] = g;
        p[g] = (p[g] + p[m]) - 1.0;
        if (p[g] >= 1.0 || mm <= g) {
            for (m = mm; m < n && p[m] >= 1.0; ++m)
            /*noop*/;
            mm = m + 1;
        } else
            m = g;
        for (; g < n && p[g] < 1.0; ++g)
        /*noop*/;
    }

    // any bars that remain have no alias
    for (; g < n; ++g) {
        if (p[g] < 1.0)
            continue;
        h[g] = std::numeric_limits<uint64_t>::max();
        a[g] = g;
    }
    if (m < n) {
        h[m] = std::numeric_limits<uint64_t>::max();
        a[m] = m;
        for (m = mm; m < n; ++m) {
            if (p[m] > 1.0)
                continue;
            h[m] = std::numeric_limits<uint64_t>::max();
            a[m] = m;
        }
    }
}

