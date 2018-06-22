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

    ind[0].setParam(param);
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
    double fit = 0;
    for (int i = 0; i < popSize; ++i){
        indFitness[i] = fit = ind[i].getFitness();
    }
}

// Do everything on population in one loop

void Population::reproduceMutateCalcFit(Population& oldPop)
{
    oldPop.createAliasTable(indFitness.data(), hvec.data(), avec.data(), popSize);
    for (int i = 0; i < popSize; ++i){
        SetBaby(oldPop.chooseInd(), oldPop.chooseInd(), ind[i]);
        ind[i].mutate();
        indFitness[i] = ind[i].getFitness();
	}
}

void Population::calcStats(Param& param, SumStat& stats)
{
    int i, j;
    std::vector<std::vector<double>> gMatrix(param.loci,std::vector<double>(param.popsize));
    
    for (i = 0; i < popSize; ++i){
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
    
    double mean = vecMean<double>(indFitness);
    stats.setAveFitness(mean);
    stats.setSDFitness(vecSD<double>(indFitness, mean));
    
    auto& fitnessDistn = stats.getFitnessDistn();
    fitnessDistn = percentiles_interpol<std::vector<double>>(indFitness, ptiles);
}

// Alias method for sampling from discrete distribution, see https://pandasthumb.org/archives/2012/08/lab-notes-the-a.html and https://en.wikipedia.org/wiki/Alias_method and http://www.keithschwarz.com/darts-dice-coins/
// Example code for testing in ~/sim/02_SmallTests/aliasMethod.cc

uint32_t Population::getRandIndex(){
    uint64_t u = rnd.rawint();
    uint32_t x = u % popSize;
    return (u < hvec[x]) ? x : avec[x];
}

void Population::createAliasTable(const double *pp, uint64_t *h, std::uint32_t *a, uint32_t n) {
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

// array is cumulative fitness or weighting, n is number of elements, returns random index from array sampled by cumulative success weighting, this is faster than using std::lower_bound as binary search, see below

int Population::chooseMember(double *array, int n)
{
    double cumFitness = array[n-1];
    if (cumFitness < 1e-8)
        return static_cast<int>(rnd.rtop(n));
    
    double ran = rnd.rU01();
    double target = ran * cumFitness;
    int k = static_cast<int>((ran * (n-1)));
    
    while (array[k++] < target);
    while ((--k != -1) && (array[k] >= target));
    return ++k;
}

// lower_bound does binary search to find iterator to least index that is >= target, distance gives the index value as an int, and then call ind[] to get individual associated with index, slower than my chooseMember code, not used but shown here to document alternative algorithms. Tested and gives same output as chooseMember. Test done with code as inlined in Population.h, but still much slower than chooseMember.

//Individual& Population::chooseInd()
//{
//    auto it = std::lower_bound(indFitness.begin(), indFitness.end(), rnd.rU01() * indFitness.back());
//    return ind[std::distance(indFitness.begin(), it)];
//}



