
#ifndef _sensitivity_h
#define _sensitivity_h 1

#include "SAFrand_pcg.h"

// pcg32 or pcg64 for pcgT; using 32bit not tested, use care with seeds and test
using pcgT = pcg64;
using rndType = pcgT::result_type;
extern SAFrand_pcg<pcgT> rnd;
extern const int linesPerRun;
extern bool showProgress;

std::string Control(std::istringstream& parmBuf);

#include "typedefs.h"
#include "util.h"       // includes percentiles, rounding of floats, ThrowError()

#endif
