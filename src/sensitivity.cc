#include <cmath>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <climits>

#include APPL_H
#include "fmt/format.h"
#include "SAFrand_pcg.h"
#include "SAFtimer.h"

#include "Population.h"
#include "Performance.h"

const int 	linesPerRun = 3;
SAFrand_pcg<pcgT> rnd;

// start with result and fix all other strings and files

/*************************** prototypes **************************/

void 		GetParam(Param& p, std::istringstream& parmBuf);
void 		LifeCycle(Param& param, std::ostringstream& resultss);
std::string PrintParam(Param& p);
void        PrintSummary(Param& param, std::ostringstream& resultss, SumStat& stats);

/*****************************************************************/


// paramBuf holds parameters for runs.  First three entries are first and last run, and random seed. The remaining data are the same as in the standard parm file and can be parsed accordingly.

std::string Control(std::istringstream& parmBuf)
{
    Param param;		
    int i, first, last;
    std::ostringstream resultss;
    
    std::string temp{parmBuf.str()};
    parmBuf >> first >> last >> param.rndSeed;
    if (parmBuf.bad())
        ThrowError(__FILE__, __LINE__, "Failed reading from parameter string stream.");
    // seed may be 64bit, but rndType may be 32 bit, if so, truncate seed
    if (sizeof(rndType) != sizeof(unsigned long)){
        rndType seed = static_cast<rndType>(param.rndSeed);
        param.rndSeed = seed;
        rnd.setRandSeed(seed);
    }
    else{
        rnd.setRandSeed(param.rndSeed);
    }
    setGSLErrorHandle(1);       // 1 => turn on my error handler, 0 => turn off handler
	for (i = first; i <= last; i++){
		GetParam(param, parmBuf);
		LifeCycle(param, resultss);
	}
	return resultss.str();
}

void LifeCycle(Param& param, std::ostringstream& resultss)
{
    if (showProgress){
        std::cout << fmt::format("\nrunNum = {:>3}\n\n", param.runNum);
        std::cout.flush();
    }
    Population p1(param);
    Population p2(param);
    Population *op, *np, *swap;     // oldpop and newpop
    op = &p1;
    np = &p2;
    SumStat stats;
    stats.initialize(param);
    int gen = param.gen;
    int i;

    op->setFitnessArray();
    for (i = 0; i < gen; ++i){
        if (showProgress && ((i % 100) == 0))
            std::cout << fmt::format("Rep {:8} of {:8}\n", i, param.gen);
        np->reproduceMutateCalcFit(*op);
        swap = op;
        op = np;
        np = swap;
    }
    // run round of selection without mutation or recombination before collecting stats
    np->reproduceNoMutRec(*op);
    np->calcStats(param, stats);
    PrintSummary(param, resultss, stats);
}

void GetParam(Param& p, std::istringstream& parmBuf)
{
    double tmpLoop, tmpGen, tmpPop, tmpMutLoc;
    
    parmBuf >> p.runNum >> tmpLoop >> tmpGen >> tmpPop >> p.mutation >> p.recombination >> p.mutStep >>
        p.aSD >>p.fitVar >> p.gamma >> p.stochWt >> tmpMutLoc;
    if (parmBuf.bad())
        ThrowError(__FILE__, __LINE__, "Failed reading from parameter string stream.");

    p.loop = static_cast<Loop>(round<int>(tmpLoop));
    p.loci = (p.loop == Loop::dclose) ? 7 : 5;
    p.gen = round<int>(tmpGen);
    p.popsize = round<int>(tmpPop);
    p.stoch = (abs(p.stochWt) < 1e-6) ? false : true;
    p.mutLocus = round<int>(tmpMutLoc);
    if (p.mutLocus >= p.loci)
        ThrowError(__FILE__, __LINE__, "Mutated locus number greater than number loci.");

    int newseed;
    unsigned long seed;
    parmBuf >> p.distnSteps >> seed >> newseed;
    if (parmBuf.bad())
        ThrowError(__FILE__, __LINE__, "Failed reading from parameter string stream.");
    
    std::cout << p.rndSeed << std::endl;
    // seed may be 64bit, but rndType may be 32 bit, if so, truncate seed
    if (!newseed){
        if ((sizeof(rndType) != sizeof(unsigned long)) && seed > UINT32_MAX)
            ThrowError(__FILE__, __LINE__, "64 bit seed in design file > 32 bit generator type");
        else{
            p.rndSeed = seed;
            rnd.setRandSeed(p.rndSeed);	// initial random numbers
        }
    }
    std::cout << seed << " " << p.rndSeed << std::endl;

}

std::string PrintParam(Param& p)
{
    std::string outString;
    std::string format = "{:<10} = {:>9}\n";
    std::string formatf = "{:<10} = {:>9.3e}\n";
    outString += fmt::format(format,  "Run", p.runNum);
    outString += fmt::format(format,  "disStp", p.distnSteps);
    outString += fmt::format("{:<10} = {}\n",  "seed", p.rndSeed);
    outString += fmt::format(format,  "loop", static_cast<int>(p.loop));
    outString += fmt::format(format,  "gen", p.gen);
    outString += fmt::format(format,  "loci", p.loci);
    outString += fmt::format(format,  "popSz", p.popsize);
    outString += fmt::format(formatf, "mut", p.mutation);
    outString += fmt::format(formatf, "rec", p.recombination);
    outString += fmt::format(format,  "recT", p.rec);
    outString += fmt::format(formatf, "mutStep", p.mutStep);
    outString += fmt::format(formatf, "aSD", p.aSD);
    outString += fmt::format(formatf, "fitVar", p.fitVar);
    outString += fmt::format(formatf, "gamma", p.gamma);
    outString += fmt::format(formatf, "stochWt", p.stochWt);
    outString += fmt::format(format,  "mutLocus", p.mutLocus);
    outString += "\n";
    return outString;
}

void PrintGDistn(Param& param, std::ostringstream& resultss,
                 std::vector<double> mean, std::vector<double> sd, std::vector<std::vector<double>> distn)
{
    resultss <<  "     ";
    for (int i = 0; i < param.loci; ++i){
        resultss << fmt::format((i < 10) ? "{:>7}{:1}" : "{:>6}{:2}", "g", i);
    }
    resultss << "\n";
    
    resultss <<  " Mean";
    for (int i = 0; i < param.loci; ++i){
        resultss << fmt::format("{:8.3f}", mean[i]);
    }
    resultss << "\n";
    
    resultss <<  "   SD";
    for (int i = 0; i < param.loci; ++i){
        resultss << fmt::format("{:8.3f}", sd[i]);
    }
    resultss << "\n";
    
    for (int j = 0; j < param.distnSteps; ++j){
        resultss << fmt::format("{:5.1f}", 100.0*(double)j/(double)(param.distnSteps-1));
        for (int i = 0; i < param.loci; ++i){
            resultss << fmt::format("{:8.3f}", distn[i][j]);
        }
        resultss << "\n";
    }
    
    resultss << "\n\n";
}

void PrintGCorr(Param& param, std::ostringstream& resultss, std::vector<std::vector<double>> corr)
{
    resultss <<  "     ";
    for (int i = 0; i < param.loci; ++i){
        resultss << fmt::format((i < 10) ? "{:>6}{:1}" : "{:>5}{:2}", "g", i);
    }
    resultss << "\n";
    
    for (int i = 0; i < param.loci; ++i){
        resultss << fmt::format((i < 10) ? "{:>4}{:1}" : "{:>3}{:2}", "g", i);
        for (int j = 0; j < param.loci; ++j){
            resultss << fmt::format("{:7.3f}", corr[i][j]);
        }
        resultss << "\n";
    }
    
    resultss << "\n\n";
}

void PrintSummary(Param& param, std::ostringstream& resultss, SumStat& stats)
{
    resultss << PrintParam(param);

    // print fitness distn
    
    resultss << fmt::format("{}", "Fitness distribution\n\n");
    resultss << fmt::format("{:5}{:11.3e}\n", " Mean", stats.getAveFitness());
    resultss << fmt::format("{:5}{:11.3e}\n", "   SD", stats.getSDFitness());
    
    const auto& fitness = stats.getFitnessDistn();
    
    for (int j = 0; j < param.distnSteps; ++j){
        resultss << fmt::format("{:5.1f}", 100.0*(double)j/(double)(param.distnSteps-1));
        resultss << fmt::format("{:11.3e}\n", fitness[j]);
    }
    resultss << "\n";
    
    // print performance distn
    
    resultss << fmt::format("{}", "Performance distribution\n\n");
    resultss << fmt::format("{:5}{:11.3e}\n", " Mean", stats.getAvePerf());
    resultss << fmt::format("{:5}{:11.3e}\n", "   SD", stats.getSDPerf());
    
    const auto& perf = stats.getPerfDistn();
    
    for (int j = 0; j < param.distnSteps; ++j){
        resultss << fmt::format("{:5.1f}", 100.0*(double)j/(double)(param.distnSteps-1));
        resultss << fmt::format("{:11.3e}\n", perf[j]);
    }
    resultss << "\n";
    
    // print genotype distn statistics for each locus
    
    const auto& gMean = stats.getGMean();
    const auto& gSD = stats.getGSD();
    const auto& gDistn = stats.getGDistn();
    
    resultss << "Distribution of genotype values, rows are percentiles\n\n";
    PrintGDistn(param, resultss, gMean, gSD, gDistn);
    
    // print G corr matrix
    
    const auto& gCorr = stats.getGCorr();
    
    resultss << "Correlation of G values\n\n";
    PrintGCorr(param, resultss, gCorr);
    
    if (param.stoch){
        // print stoch distn statistics for each locus
        
        const auto& mean = stats.getSMean();
        const auto& sd = stats.getSSD();
        const auto& distn = stats.getSDistn();
        
        resultss << "Distribution of stochastic values, rows are percentiles\n\n";
        PrintGDistn(param, resultss, mean, sd, distn);
        
        // print S corr matrix
        
        const auto& corr = stats.getSCorr();
        
        resultss << "Correlation of stoch values\n\n";
        PrintGCorr(param, resultss, corr);
        
        // print S X G corr matrix
        
        const auto& corrsg = stats.getSGCorr();
        
        resultss << "Correlation of stoch x genotype values\n\n";
        PrintGCorr(param, resultss, corrsg);
    }
}

