// Method for converting transfer function to system of first order differential equations for evaluation numerically is in Ogata, 4th ed, page 78.
// Use Mathematica StateSpaceModel instead of Ogata as basis for conversion of TF to State Space.

// Step response achieved by setting u(t) = 1, in derivative term, see MMA State Space model in Plasticity.nb, subsection on Utilities for preparing C++ code.

// Chop p0 in calling routine
// Bound p1 >= 1e-7 in calling routine

#include <iostream>
#include <vector>
#include <cassert>
#include <string>
#include <map>

#include <gsl/gsl_poly.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_spline.h>

#include <boost/math/tools/polynomial.hpp>

#include "fmt/format.h"
#include "Performance.h"

// steps for interpolation, 5000 comes very close to Mathematica numerical results, check timing
const int steps = 5000;

unsigned debugPerformance = 0;
void debugPerformanceOn(unsigned d) {debugPerformance=d;}

double 	MaxRootRealPart(const std::vector<double>& coeff);
double 	H2sq(const std::vector<double>& num, const std::vector<double>& den);
double 	integrandH2(double w, void *p);
int 	deriv (double t, const double x[], double f[], void *p);
double 	stepPerformance(const std::vector<double>& num, const std::vector<double>& den, double tmax, signalType s);
double 	integrandStep(double y, void *p);

struct my_params {const std::vector<double>& num; const std::vector<double>& den;};
struct stepParam {gsl_spline *spline; gsl_interp_accel *acc;};

// Do chopping to set param to zero if close and check bounds before call

double performance(const std::vector<double>& num, const std::vector<double>& den,
					const std::vector<double>& dentilde, double gamma, double tmax, signalType s)
{
	if (MaxRootRealPart(den) > -1e-6) return 1e20;
    // next line checks for stability when Plant parameter a is variant, if system is unstable, then performance failure, can turn off if not checking for stability margin
	if (MaxRootRealPart(dentilde) > 1e-6) return 1e20;
	if (debugPerformance){
		double sf = stepPerformance(num, den, tmax, s);
		double hf = H2sq(num,den);
		std::map<signalType,std::string> signal = 
			{{signalType::output, "output"}, {signalType::controlOpen,"cntrlO"}, {signalType::controlClosed,"cntrlC"}};
		std::cout << fmt::format("step = {}, h2 = {} for {}\n", sf, hf, signal[s]);
		return sf + gamma*hf;
	}
	else
		return stepPerformance(num, den, tmax, s) + gamma*H2sq(num,den);
}

// coeff of polynomial from low order to high order terms
double MaxRootRealPart(const std::vector<double>& coeff) 
{
	auto n = coeff.size();
	// std::cout << "n = " << n << std::endl;
    std::vector<double> z(2*(n-1));
	gsl_poly_complex_workspace * w
	  = gsl_poly_complex_workspace_alloc (n);
	if (GSL_SUCCESS != gsl_poly_complex_solve(coeff.data(), n, w, z.data())){
		gsl_poly_complex_workspace_free (w);
		return 1.0;		// positive value is marked as unstable
	}
	gsl_poly_complex_workspace_free (w);

	double max = -1e20;
	for (unsigned i = 0; i < n-1; i++) {
	  // std::cout << fmt::format("z{} = {:+.18f} {:+.18f}\n", i, z[2*i], z[2*i+1]);
	  if (z[2*i] > max) max = z[2*i];
	}
	return max;
}

// In the following, check that integrand gets small at large frequency, required for H2 integral to converge. If not then if order of numerator = order of denominator, then integration will be infinite, because at high frequency, high order terms in numerator and denominator dominate, and so given same order of those terms, the system does not go to zero at high frequency. Correct for this by redefining the num/den transfer function as (num/den - high order num coeff/high order den coeff). This subtraction removes the Dirac delta impulse component of the dynamics, which corresponds in frequency space to a uniform addition to the absolute value of the tf at all frequencies equal to the amount subtracted off. Does not change dynamics, except to remove impulse at time zero acting over infinitesimal duration.

// returns square of H2 value
double H2sq(const std::vector<double>& num, const std::vector<double>& den)
{
	bool sizeFix = false;
	boost::math::tools::polynomial<double> poly_newnum;
	boost::math::tools::polynomial<double> poly_newden;
	if (num.size() == den.size()){			// deleting dirac impulse for equal-sized num & den
		double numHiOrder = num.back();
		double denHiOrder = den.back();
		boost::math::tools::polynomial<double> polyn(num.begin(), num.end());
		boost::math::tools::polynomial<double> polyd(den.begin(), den.end());
		poly_newnum = denHiOrder*polyn - numHiOrder*polyd;
		poly_newden = denHiOrder*polyd;
		// high order term in new numerator should be zero and automatically dropped, reducing size
		if (poly_newnum.data().size() != polyn.data().size() - 1){
			if (poly_newnum.data().back() < 1e-5)
				poly_newnum.data().pop_back();
			else
				assert(poly_newnum.data().size() == polyn.data().size() - 1);
		}
		sizeFix = true;
// 		std::cout << "Num = " << polyn << std::endl;
// 		std::cout << "Den = " << polyd << std::endl;
// 		std::cout << "NN  = " << poly_newnum << std::endl;
// 		std::cout << "ND  = " << poly_newden << std::endl;
	}
	
	gsl_function F;
	F.function = &integrandH2;
    // if size fixed above, use newnum and newden
    const std::vector<double>& n = (sizeFix) ? poly_newnum.data() : num;
    const std::vector<double>& d = (sizeFix) ? poly_newden.data() : den;
    my_params params {n, d};
    F.params = &params;
	if (integrandH2(1e10, F.params) > 1e-3) return 1e20; 	// should not happen, because den.size > num.size
	double result, error;
	gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);
	// std::cout << "h2 int start" << std::endl;
	if (GSL_SUCCESS != gsl_integration_qagi(&F, 0, 1e-7, 1000, w, &result, &error))
		result = 1e30;
	// std::cout << "h2 int end" << std::endl;
	gsl_integration_workspace_free(w);
	return result / (2.0*M_PI);
}

double integrandH2(double w, void *p)
{
	gsl_complex s;
	my_params *params = static_cast<my_params *>(p);
	const std::vector<double>& num = params->num;
	const std::vector<double>& den = params->den;
	GSL_SET_COMPLEX(&s, 0, w);
    int numSize = static_cast<int>(num.size());
    int denSize = static_cast<int>(den.size());
	return gsl_complex_abs2(gsl_complex_div(gsl_poly_complex_eval(num.data(), numSize, s),
								gsl_poly_complex_eval(den.data(), denSize, s)));
}

// deriv works for both output signal and control signal, see MMA file

int deriv(double t, const double x[], double f[], void *p)
{
	(void)(t); /* avoid unused parameter warning */
	my_params *params = static_cast<my_params *>(p);
	const std::vector<double>& den = params->den;
	auto lastrow = den.size()-2;
	f[lastrow] = 0;
	for (unsigned i = 0; i < lastrow; ++i){
		f[i] = x[i+1];
		f[lastrow] -= den[i]*x[i];
	}
	f[lastrow] -= den[lastrow]*x[lastrow];
	f[lastrow] /= den.back();
	// next line for input, u(t) = 1 for step input
	f[lastrow] += 1;
	return GSL_SUCCESS;
}

// must make different ycoeff for output and control signals
// for output, same coeff for open and closed loops
// for control signals, diff coeff for open and closed loops
// see MMA file

double stepPerformance(const std::vector<double>& num, const std::vector<double>& den, double tmax, signalType s)
{
	double time[steps+1];
	double y[steps+1];
	my_params params = {num, den};
	auto dim = den.size()-1;	// dimensions of state space model for dynamics
	gsl_odeiv2_system sys = {deriv, NULL, dim, &params};
	gsl_odeiv2_driver *d =
		gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);
	double t = 0.0;
	double* x = new double[den.size()-1](); // () at end causes zero initialization
	double denBack = den.back();
	
	// coefficients to get output, initialize with values for each case, max dim is 4, so use that
	double ycoeff[4];
	unsigned long ydim;		// number of output coefficients to get output y, varies by problem, set explicitly
	double yinputCoeff = 0;	// must add yinputCoeff * input to output; input=1 for the step response 
	
	if (s == signalType::output){ // case of output signal, same for open and closed loops
		ydim = 3; 
		ycoeff[0] = num[0]/denBack; ycoeff[1] = num[1]/denBack; ycoeff[2] = num[2]/denBack;
	}
	else {
		// cases of control signal, # coeff is always dimension of problem, dim
		ydim = dim;
		double numBack = num.back();
		yinputCoeff = numBack / denBack;
		if (s == signalType::controlOpen) { // control signal, open loop
			for (unsigned i = 0; i < ydim; ++i)
				ycoeff[i] = (num[i] - den[i]*yinputCoeff) / denBack;
		}
		else if (s == signalType::controlClosed) { // control signal, closed loop
			for (unsigned i = 0; i < ydim; ++i){
				ycoeff[i] = (num[i] - den[i]*yinputCoeff) / denBack;
				if (debugPerformance >= 2) std::cout << ycoeff[i] << " ";
			}
			if (debugPerformance >= 2) std::cout << yinputCoeff << std::endl;
			if (debugPerformance >= 2){
				for (auto& n : den) std::cout << -n/denBack << " ";
				std::cout << std::endl;
				for (auto& n : num) std::cout << n << " ";
				std::cout << std::endl;
				for (auto& n : den) std::cout << n << " ";
				std::cout << std::endl;
			}
		}
		else {
			assert(false);	// should not be here, signal must be one of above types
		}
	}

	time[0] = 0.0;
	// at time zero with step, dominated by infinite freq, so if order of den > num, then at time zero,
	// initial value is zero, if order den=num, then ratio of highest order terms,
	// order num>den should not happen
	y[0] = (den.size() > num.size()) ? 0.0 : num.back()/den.back(); 
	for (int i = 1; i <= steps; i++){
		// add one to tmax, so that interpolation extends past boundary for integration
		// otherwise, can get an error when interpolating close to the upper boundary
	  	double ti = i * (tmax+1.0) / static_cast<double>(steps);
	  	int status = gsl_odeiv2_driver_apply(d, &t, ti, x);
		assert(status == GSL_SUCCESS);
		time[i] = t;
		y[i] = 0;
		for (unsigned j = 0; j < ydim; ++j)
			y[i] += ycoeff[j] * x[j];
		y[i] += yinputCoeff;
		// if (i % 100 == 0) std::cout << fmt::format("{:7.3f} {:8.6f}\n", time[i], y[i]);
		if (debugPerformance >= 3 && i % 100 == 0 && s == signalType::controlClosed) 
			std::cout << fmt::format("{:7.3f} {:8.6f}\n", time[i], y[i]);
	}
	gsl_odeiv2_driver_free (d);
	delete [] x;
	
	gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, steps);
    if (GSL_SUCCESS != gsl_spline_init (spline, time, y, steps)){
    	gsl_spline_free (spline);
    	gsl_interp_accel_free (acc);
    	return 1e20;
    }
    
    // std::cout << gsl_spline_eval(spline, 19.9689, acc) << std::endl;
    
	stepParam integrParam = {spline, acc};
    gsl_function F;
	F.function = &integrandStep;
	F.params = &integrParam;
	double result, error;
	gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);
	// std::cout << "step int start" << std::endl;
	// using 1e-6 for abs and rel error
	double errtol = 1e-6;
	if (GSL_SUCCESS != gsl_integration_qag(&F, 0.0, tmax, errtol, errtol, 1000, 6, w, &result, &error)){
		gsl_integration_cquad_workspace *ctable = gsl_integration_cquad_workspace_alloc(200);
		if (GSL_SUCCESS != gsl_integration_cquad(&F, 0, tmax, errtol, errtol, ctable, &result, &error, NULL)){
     		boost::math::tools::polynomial<double> poly_newnum(num.begin(), num.end());
			boost::math::tools::polynomial<double> poly_newden(den.begin(), den.end());
			if (true){
				std::cout << "Step integration error, with num = " << poly_newnum << std::endl;
				std::cout << "                         and den = " << poly_newden << std::endl;
				std::cout << "                      and result = " << result << std::endl;
			}
    		result = 1e20;
    	}
    	gsl_integration_cquad_workspace_free(ctable);
	}
	// std::cout << "step int end" << std::endl;
	gsl_integration_workspace_free(w);

    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);

	return result;
}

double integrandStep(double t, void *p)
{
	stepParam *params = static_cast<stepParam *>(p);
	double z = 1.0 - gsl_spline_eval(params->spline, t, params->acc);
	return z*z;
}
