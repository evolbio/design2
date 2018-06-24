#ifndef Performance_h
#define Performance_h 1
#include <vector>

// if p0 is less than chop, then set to zero, associated with calling routine chop
const double chop = 1e-3;
const double p1bound = 1e-5;

void debugPerformanceOn(unsigned);

enum class signalType {output, controlOpen, controlClosed};

double performance(const std::vector<double>& num, const std::vector<double>& den, 
					const std::vector<double>& dentilde,
					double gamma, double tmax, signalType s);


#endif
