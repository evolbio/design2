
#ifndef _sensitivity_h
#define _sensitivity_h 1

#include "SAFrand_pcg.h"

using pcgT = pcg64;     // pcg32 or pcg64
extern SAFrand_pcg<pcgT> rnd;
extern const int linesPerRun;
extern bool showProgress;

std::string Control(std::istringstream& parmBuf);

#include "typedefs.h"
#include "util.h"       // includes percentiles, rounding of floats, ThrowError()

#endif
