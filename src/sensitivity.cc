#include <cmath>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <climits>

#include "fmt/format.h"
#include "sensitivity.h"
#include "SAFrand_pcg.h"
#include "SAFtimer.h"

#include "Population.h"

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
    rndType seed;
    std::ostringstream resultss;
    
    std::string temp{parmBuf.str()};
    parmBuf >> first >> last >> seed;
    if (parmBuf.bad())
        ThrowError(__FILE__, __LINE__, "Failed reading from parameter string stream.");
	rnd.setRandSeed(seed);
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
        if (showProgress && ((i % 1000) == 0))
            std::cout << fmt::format("Rep {:8} of {:8}\n", i, param.gen);
        np->reproduceMutateCalcFit(*op);
        swap = op;
        op = np;
        np = swap;
    }
    op->calcStats(param, stats);
    PrintSummary(param, resultss, stats);
}

void GetParam(Param& p, std::istringstream& parmBuf)
{
    double tmpGen, tmpLoci, tmpPop;
    
    parmBuf >> p.runNum >> tmpGen >> tmpLoci >> tmpPop >> p.mutation >> p.recombination >> p.maxAllele >> p.fitVar;
    if (parmBuf.bad())
        ThrowError(__FILE__, __LINE__, "Failed reading from parameter string stream.");

    p.gen = round<int>(tmpGen);
    p.loci = round<int>(tmpLoci);
    p.popsize = round<int>(tmpPop);
    
    rndType newseed;
    rndType parmSeed;
    parmBuf >> p.distnSteps >> parmSeed >> newseed;
    if (parmBuf.bad())
        ThrowError(__FILE__, __LINE__, "Failed reading from parameter string stream.");
	if (!newseed)
		rnd.setRandSeed(parmSeed);	// initial random numbers
}

std::string PrintParam(Param& p)
{
    std::string outString;
    std::string format = "{:<10} = {:>9}\n";
    std::string formatf = "{:<10} = {:>9.3e}\n";
    outString += fmt::format(format, "Run", p.runNum);
    outString += fmt::format(format, "disStp", p.distnSteps);
    outString += fmt::format(format, "gen", p.gen);
    outString += fmt::format(format, "loci", p.loci);
    outString += fmt::format(format, "popSz", p.popsize);
    outString += fmt::format(formatf, "mut", p.mutation);
    outString += fmt::format(formatf, "rec", p.recombination);
    outString += fmt::format(format, "recT", p.rec);
    outString += fmt::format(formatf, "maxAllele", p.maxAllele);
    outString += fmt::format(formatf, "fitVar", p.fitVar);
    outString += "\n";
    return outString;
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

    // print genotype distn statistics for each locus
    
    const auto& gMean = stats.getGMean();
    const auto& gSD = stats.getGSD();
    const auto& gDistn = stats.getGDistn();
    
    resultss << "Distribution of genotype values, rows are percentiles\n\n";
    
    resultss <<  "     ";
    for (int i = 0; i < param.loci; ++i){
        resultss << fmt::format((i < 10) ? "{:>7}{:1}" : "{:>6}{:2}", "g", i);
    }
    resultss << "\n";
    
    resultss <<  " Mean";
    for (int i = 0; i < param.loci; ++i){
        resultss << fmt::format("{:8.3f}", gMean[i]);
    }
    resultss << "\n";

    resultss <<  "   SD";
    for (int i = 0; i < param.loci; ++i){
        resultss << fmt::format("{:8.3f}", gSD[i]);
    }
    resultss << "\n";

    for (int j = 0; j < param.distnSteps; ++j){
        resultss << fmt::format("{:5.1f}", 100.0*(double)j/(double)(param.distnSteps-1));
        for (int i = 0; i < param.loci; ++i){
            resultss << fmt::format("{:8.3f}", gDistn[i][j]);
        }
        resultss << "\n";
    }
    
    resultss << "\n\n";

    // print G corr matrix
    
    const auto& gCorr = stats.getGCorr();
    
    resultss << "Correlation of G values\n\n";
    
    resultss <<  "     ";
    for (int i = 0; i < param.loci; ++i){
        resultss << fmt::format((i < 10) ? "{:>6}{:1}" : "{:>5}{:2}", "g", i);
    }
    resultss << "\n";

    for (int i = 0; i < param.loci; ++i){
        resultss << fmt::format((i < 10) ? "{:>4}{:1}" : "{:>3}{:2}", "g", i);
        for (int j = 0; j < param.loci; ++j){
            resultss << fmt::format("{:7.3f}", gCorr[i][j]);
        }
        resultss << "\n";
    }
    
    resultss << "\n\n";

}

