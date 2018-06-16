#include <cmath>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <random>
#include <climits>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include "format.h"
#include "sensitivity.h"
#include "SAFrand_pcg.h"
#include "SAFrand_mt.h"
#include "SAFtimer.h"
#include "pcg_random.hpp"

const int 	linesPerRun = 4;
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
    double tmpOincr, tmpSinit, tmpXinit, tmpRinit;
    
    parmBuf >> p.runNum >> p.step >> tmpOincr >> tmpSinit >> tmpXinit >> tmpRinit >> p.diff >> p.ds >> p.dx >> p.da >> p.dr >> p.rsxf >> p.rsxb >> p.rapr >> p.radc >> p.gmsx;
    if (parmBuf.bad())
        ThrowError(__FILE__, __LINE__, "Failed reading from parameter string stream.");

    p.oincr = round<int>(tmpOincr);
    p.sinit = round<int>(tmpSinit);
    p.xinit = round<int>(tmpXinit);
    p.rinit = round<int>(tmpRinit);
    p.ds *= p.diff;
    p.dx *= p.diff;
    p.da *= p.diff;
    p.dr *= p.diff;
    
    int tmpBndry, newseed;
    unsigned parmSeed;
    parmBuf >> p.time >> p.tinit >> tmpBndry >> parmSeed >> newseed;
    if (parmBuf.bad())
        ThrowError(__FILE__, __LINE__, "Failed reading from parameter string stream.");
    p.bndry = (tmpBndry == 0) ? 'p' : 'r';
	if (!newseed)
		rnd.setRandSeed(parmSeed);	// initial random numbers
}

std::string PrintParam(Param& p)
{
    std::string outString;
    std::string format = "{:<10} = {:>9}\n";
    std::string formatf = "{:<10} = {:>9.3e}\n";
    outString += fmt::format(format, "Run", p.runNum);
    outString += fmt::format(format, "time", p.time);
    outString += fmt::format(format, "tinit", p.tinit);
    outString += fmt::format(formatf, "step", p.step);
    outString += fmt::format(format, "bndry", p.bndry);
    outString += fmt::format(format, "oincr", p.oincr);
    outString += fmt::format(format, "sinit", p.sinit);
    outString += fmt::format(format, "xinit", p.xinit);
    outString += fmt::format(format, "rinit", p.rinit);
    outString += fmt::format(formatf, "ds", p.ds);
    outString += fmt::format(formatf, "dx", p.dx);
    outString += fmt::format(formatf, "da", p.da);
    outString += fmt::format(formatf, "dr", p.dr);
    outString += fmt::format(formatf, "rsxf", p.rsxf);
    outString += fmt::format(formatf, "rsxb", p.rsxb);
    outString += fmt::format(formatf, "rapr", p.rapr);
    outString += fmt::format(formatf, "radc", p.radc);
    outString += fmt::format(formatf, "gmsx", p.gmsx);
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

