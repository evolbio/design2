
#ifndef _sensitivity_h
#define _sensitivity_h 1

#include "SAFrand_pcg.h"

using pcgT = pcg64;     // pcg32 or pcg64
using rndType = pcgT::result_type;
extern SAFrand_pcg<pcgT> rnd;
extern const int linesPerRun;
extern bool showProgress;

std::string Control(std::istringstream& parmBuf);

#include "typedefs.h"
#include "util.h"       // includes percentiles, rounding of floats, ThrowError()

#endif
