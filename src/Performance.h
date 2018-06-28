#ifndef Performance_h
#define Performance_h 1
#include <vector>

#include <gsl/gsl_errno.h>

// if p0 is less than chop, then set to zero, associated with calling routine chop
const double chop = 1e-3;
const double p1bound = 1e-5;

void debugPerformanceOn(unsigned);

enum class signalType {output, controlOpen, controlClosed};

double performance(const std::vector<double>& num, const std::vector<double>& den, 
					double gamma, double tmax, signalType s);

// GSL error handling: call setGSLErrorHandle(s), s = 0 turns off error handler, 1 sets my handler; should check return status of all significant GSL calls and take appropriate action within code, for example return high performance value and thus zero fitness if cannot evaluate performance for parameter combination
inline void my_gsl_handler (const char *reason, const char *file, int line, int gsl_errno __attribute__((unused)))
{std::cout << fmt::format("GSL error: {}:{}, {}\n", file, line, reason);}
inline void setGSLErrorHandle(int s)
    {(s) ? gsl_set_error_handler(&my_gsl_handler) : gsl_set_error_handler_off();}

#endif
