#include <cmath>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <random>
#include <climits>

#include "fmt/format.h"
#include "sensitivity.h"
#include "SAFrand_pcg.h"
#include "SAFtimer.h"

const int 	linesPerRun = 3;
SAFrand_pcg<pcgT> rnd;

// start with result and fix all other strings and files

/*************************** prototypes **************************/

void 		GetParam(Param& p, std::istringstream& parmBuf);
void 		LifeCycle(Param& param, std::ostringstream& resultss);
std::string PrintParam(Param& p);
void        PrintSummary(Param& param, std::ostringstream& resultss);

/*****************************************************************/


// paramBuf holds parameters for runs.  First three entries are first and last run, and random seed. The remaining data are the same as in the standard parm file and can be parsed accordingly.

std::string Control(std::istringstream& parmBuf)
{
    Param param;		
    int i, first, last;
    unsigned seed;
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
    PrintSummary(param, resultss);
}

void GetParam(Param& p, std::istringstream& parmBuf)
{
    double tmpGen, tmpLoci, tmpPop;
    
    parmBuf >> p.runNum >> tmpGen >> tmpLoci >> tmpPop >> p.mutation >> p.recombination >> p.s >> p.minFail >> p.failExp;
    if (parmBuf.bad())
        ThrowError(__FILE__, __LINE__, "Failed reading from parameter string stream.");

    p.gen = round<int>(tmpGen);
    p.loci = round<int>(tmpLoci);
    p.popsize = round<int>(tmpPop);
    
    int newseed;
    unsigned parmSeed;
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
    outString += fmt::format(formatf, "s", p.s);
    outString += fmt::format(formatf, "minFail", p.minFail);
    outString += fmt::format(formatf, "failExp", p.failExp);
    outString += "\n";
    return outString;
}

void PrintSummary(Param& param, std::ostringstream& resultss)
{
    resultss << PrintParam(param);
    resultss << "-------------------------\n\n";
    
    constexpr int reps = 10'000'000;
    SAFTimer mytimer;
    
    double x = 0;
    auto g = rnd.generator();
    mytimer.start();
    for (int i = 0; i <= reps; i++){
        //x += std::uniform_real_distribution<double>{0.0, 1.0}(pcg);
        x += g();
    }
    mytimer.stop();
    resultss << fmt::format("pcg32  = {:8.3e}, elapsed time = {} micro s\n", x, mytimer.time<microseconds>());
    
    resultss << "\n--------------------------\n\n";
        
   // rtop, rbit, rbitB

    x = 0;
    mytimer.start();
    for (int i = 0; i <= reps; i++){
        x += rnd.rbit();
    }
    mytimer.stop();
    resultss << fmt::format("rbit  = {:8.3e}, elapsed time = {} micro s\n", x, mytimer.time<microseconds>());
    
    x = 0;
    mytimer.start();
    for (int i = 0; i <= reps; i++){
        x += g(2);
    }
    mytimer.stop();
    resultss << fmt::format("rtop  = {:8.3e}, elapsed time = {} micro s\n", x, mytimer.time<microseconds>());
    
    x = 0;
    mytimer.start();
    for (int i = 0; i <= reps; i++){
        x += rnd.rbit30();
    }
    mytimer.stop();
    resultss << fmt::format("mask  = {:8.3e}, elapsed time = {} micro s\n", x, mytimer.time<microseconds>());
    
}

