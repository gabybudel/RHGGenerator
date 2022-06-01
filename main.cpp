#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include <cmath>
#include <math.h>
#include <cstdlib>
#include <random>
#include <chrono>
#include <ctime>   
#include <thread>
#include <algorithm>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_monte_miser.h>
#include "./lambert.h"

/* Constants */
const double pi = 3.14159265358979323846;

/* Integration parameters */
struct int_params { int nodes; int d; double radius; double tau; double aCoeff; double constantId; };

class InputParser {
	/**
		Parse Cmd line input.
	*/
public:
	InputParser(int& argc, char** argv) {
		for (int i = 1; i < argc; ++i)
			this->tokens.push_back(std::string(argv[i]));
	}

	const std::string& getCmdOption(const std::string& option) const {
		/* Get cmd option value. */
		std::vector<std::string>::const_iterator itr;
		itr = find(this->tokens.begin(), this->tokens.end(), option);
		if (itr != this->tokens.end() && ++itr != this->tokens.end()) {
			return *itr;
		}
		static const std::string empty_string("");
		return empty_string;
	}

	bool cmdOptionExists(const std::string& option) const {
		/* Check if cmd option exists. */
		return std::find(this->tokens.begin(), this->tokens.end(), option)
			!= this->tokens.end();
	}
private:
	std::vector<std::string> tokens;
};

template <typename T> int sgn(T val) {
	/* Sign function. */

	return (T(0) < val) - (val < T(0));
}

template <typename T> int heaviside(T val) {
	/* Heaviside step function. */

	return (val > T(0));
}

double lambert_lookup(double ratio) {
	/* Find corresponding value Lambert function for ratio n/nu. */

	int row = 0;
	double res = -1.0;
	while (lambert_lookup_table[row++][0] < ratio)
	{
		// Find between which two values the given ratio is located
		res = lambert_lookup_table[row][1];
	}

	if(res > 0 && (unsigned) row < sizeof(lambert_lookup_table)/sizeof(lambert_lookup_table[0]) - 1) {
		// Find Lambert value of ratio via linear interpolation
		double x1 = lambert_lookup_table[row-1][0];
		double x2 = lambert_lookup_table[row][0];
		res = res + ((ratio - x1)/(x2 - x1))*(lambert_lookup_table[row][1] - res);
	} else {
		res = -1.0;
	}
	return res;
}

/* Constants pseudorandom generator */
#define IA 16807

#define IM 2147483647

#define AM (1.0/IM)

#define IQ 127773

#define IR 2836

#define NTAB 32

#define NDIV (1+(IM-1)/NTAB)

#define EPS 1.2e-7

#define RNMX (1.0-EPS)

float ran1(long* idum) {
	/**
		"Minimal" random number generator of Park and Miller with Bays-Durham shuffle and added
		safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint
		values). Call with idum a negative integer to initialize; thereafter, do not alter idum between
		successive deviates in a sequence. RNMX should approximate the largest floating value that is
		less than 1. 
	 */

	// Variable declaration
	int j;
	long k;
	static long iy = 0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0 || !iy) {
		// Initialize

		if (-(*idum) < 1)
			*idum = 1; /*Be sure to prevent idum = 0.*/
		else
			*idum = -(*idum);


		for (j = NTAB + 7; j >= 0; j--)
		{	// Load the shuffle table (after 8 warm-ups)
			k = (*idum) / IQ;
			*idum = IA * (*idum - k * IQ) - IR * k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}

		iy = iv[0];
	}

	k = (*idum) / IQ;	// Start here when not initializing
	*idum = IA * (*idum - k * IQ) - IR * k; // Compute idum=(IA*idum) % IM without overlows by Schrage's method
	if (*idum < 0)
		*idum += IM;

	j = iy / NDIV;	// Will be in the range 0..NTAB-1
	iy = iv[j];	// Output previously stored value and refill the shuffle table
	iv[j] = *idum;
	if ((temp = AM * iy) > RNMX)
		return RNMX;	// Because users don't expect endpoint values
	else
		return temp;
}

double avg_k(double* x, size_t dim, void* passed_parameters) {
	/* Integrand for average degree <k> integral. */
	(void)(dim); // avoid unused parameter warnings

	struct int_params* params = (struct int_params*)passed_parameters;

	double factor = (params->nodes - 1.0) * (pow(params->aCoeff, 2.0) / params->constantId);
	return factor * std::exp(params->aCoeff * (x[0] - params->radius)) * std::exp(params->aCoeff * (x[1] - params->radius)) * \
		(pow(sin(x[2]), params->d - 1.0) / (1.0 + std::exp((x[0] + x[1] - params->radius) / params->tau) * \
			pow(sin(x[2] / 2.0), params->d / params->tau)));
}

double avg_k_tzero(double* x, size_t dim, void* passed_parameters) {
	/* Integrand for average degree <k> integral at tau = 0. */
	(void)(dim); // avoid unused parameter warnings 

	struct int_params* params = (struct int_params*)passed_parameters;

	double factor = (params->nodes - 1.0) * (pow(params->aCoeff, 2.0) / params->constantId);
	return factor * std::exp(params->aCoeff * (x[0] - params->radius)) * std::exp(params->aCoeff * (x[1] - params->radius)) * \
		pow(sin(x[2]), params->d - 1.0) * heaviside(params->radius - x[0] - x[1] - (params->d * std::log(sin(x[2]/2.0))));
}

double integralAvgDegMISER(const double& radius, const int& nodes, const int& dim, const double& rescaledTemp, \
	const double& aCoeff, const double& constantId, const double& cutoff1, gsl_rng* r, const bool& debug) {
	/* Integral computation for desired average degree <k> with the MISER algorithm Monte Carlo integration. */

	if (radius == 0) {   
		return 0.0;
	}

	// Integration limits
	double xl[3] = { 0, 0, 0 };
	double xu[3] = { radius, radius, pi };

	// GSL integration objects
	gsl_monte_function integrand;
	if (rescaledTemp < cutoff1) {
		// For tau = 0 we integrate over the step function
		integrand = { &avg_k_tzero, 3, 0 };
	}
	else {
		// For tau > 0 we integrate over the regular connection probability function
		integrand = { &avg_k, 3, 0 };
	}
	struct int_params params = { nodes, dim, radius, rescaledTemp, aCoeff, constantId };
	integrand.params = &params;

	// GSL settings
	size_t calls = 5e5;

	// Monte Carlo integration MISER
	double res, err;
	gsl_monte_miser_state* s = gsl_monte_miser_alloc(3);
	gsl_monte_miser_integrate(&integrand, xl, xu, 3, calls, r, s,
		&res, &err);

	gsl_monte_miser_free(s);

	return res;

}

double integralAvgDeg(const double& radius, const int& nodes, const int& dim, const double& rescaledTemp, \
	const double& aCoeff, const double& constantId, const double& cutoff1, gsl_rng* r, const bool& debug) {
	/* Integral computation for desired average degree <k> with importance sampling Monte Carlo integration. */

	if (radius == 0) {
		return 0.0;
	}

	// Integration limits
	double xl[3] = { 0, 0, 0 };
	double xu[3] = { radius, radius, pi };

	// GSL integration objects
	gsl_monte_function integrand;
	if (rescaledTemp < cutoff1) {
		// For tau = 0 we integrate over the step function
		integrand = { &avg_k_tzero, 3, 0 };
	}
	else {
		// For tau > 0 we integrate over the regular connection probability function
		integrand = { &avg_k, 3, 0 };
	}
	struct int_params params = { nodes, dim, radius, rescaledTemp, aCoeff, constantId };
	integrand.params = &params;

	// GSL settings
	size_t warmup = 1e4, calls = 5e5;

	// Monte Carlo integration Vegas
	double res, err;
	gsl_monte_vegas_state* s = gsl_monte_vegas_alloc(3);
	gsl_monte_vegas_integrate(&integrand, xl, xu, 3, warmup, r, s, &res, &err);

	do {
		gsl_monte_vegas_integrate(&integrand, xl, xu, 3, calls, r, s, &res, &err);
	} while (fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.5);

	gsl_monte_vegas_free(s);

	return res;

}

double targetAvgDegMISER(const double& radius, const double& avgDeg, const int& nodes, const int& dim, \
	const double& rescaledTemp, const double& aCoeff, const double& constantId, const double& cutoff1, gsl_rng* r, const bool& debug) {
	/* Target function of integral computation for desired average degree <k> with the MISER algorithm Monte Carlo integration. */

	return integralAvgDegMISER(radius, nodes, dim, rescaledTemp, aCoeff, constantId, cutoff1, r, debug) - avgDeg;
}

double targetAvgDeg(const double& radius, const double& avgDeg, const int& nodes, const int& dim, \
	const double& rescaledTemp, const double& aCoeff, const double& constantId, const double& cutoff1, gsl_rng* r, const bool& debug) {
	/* Target function computation integral for desired average degree <k> with importance sampling Monte Carlo integration. */

	return integralAvgDeg(radius, nodes, dim, rescaledTemp, aCoeff, constantId, cutoff1, r, debug) - avgDeg;
}


double bisection(const double& avgDeg, const double& leftInit, const double& rightInit, const int& nodes, \
	const int& dim, const double& rescaledTemp, const double& aCoeff, const double& constantId, const double& cutoff1, \
	gsl_rng* r, const double& TOL, const bool& debug) {
	/* Compute rescaled radius numerically with bisection method for desired average degree <k>. */

	double left = leftInit, right = rightInit;
	double mid = (left + right) / 2.0;

	
	// Compute initial values with the MISER routine
	double fLeft = targetAvgDegMISER(left, avgDeg, nodes, dim, rescaledTemp, aCoeff, constantId, cutoff1, r, debug);
	double fRight = targetAvgDegMISER(right, avgDeg, nodes, dim, rescaledTemp, aCoeff, constantId, cutoff1, r, debug);
	double fMid = targetAvgDegMISER(mid, avgDeg, nodes, dim, rescaledTemp, aCoeff, constantId, cutoff1, r, debug);

	// Find a valid left starting point with different sign from right function value
	while (sgn(fRight) == sgn(fLeft)) {
		left += 0.1;
		fLeft = targetAvgDegMISER(left, avgDeg, nodes, dim, rescaledTemp, aCoeff, constantId, cutoff1, r, false);
		if (left >= right) {
			return -1.0;
		}
	}

	// Bisection method
	bool converged = false;
	int count = 0, thresh = (int)1e4;
	while (!converged) {
		if (abs(fMid) < TOL || abs(left - mid) < TOL) {
			converged = true;
		}
		else {
			if (sgn(fMid) != sgn(fLeft)) {
				right = mid;
				fRight = fMid;
			}
			else {
				left = mid;
				fLeft = fMid;
			}
			mid = (left + right) / 2.0;

			// Use the more precise Vegas Monte Carlo routine
			fMid = targetAvgDeg(mid, avgDeg, nodes, dim, rescaledTemp, aCoeff, constantId, cutoff1, r, debug);
		}

		if (++count > thresh) {
			// No solution could be found
			return -1.0;
		}
	}

	return mid;
}

double constantIdk(const int &dim, const int &k) {
	/* Compute integral I_{d,k} from closed-form solution. */

	return (std::sqrt(pi) * (std::tgamma((dim - k + 1) / 2.0) / std::tgamma(1 + (dim - k) / 2.0)));
}

double getCosAngle(int& dim, std::vector<double>& sin_theta, std::vector<double>& sin_theta_prime, std::vector<double>& cos_theta, std::vector<double>& cos_theta_prime) {
	/* Compute cosine of angle between two points. Require pre-calculated sine/cosine of angle in each dimension. */

	double sineProduct = 1;
	double cosAngle = 0;
	for (int d = 0; d < dim; d++) {
		cosAngle += sineProduct * cos_theta[d] * cos_theta_prime[d];
		sineProduct *= sin_theta[d] * sin_theta_prime[d];
	}
	cosAngle += sineProduct;

	return cosAngle;
}

double getHypDist(double& zeta, double& radial, double& radial_prime, double& cosAngle) {
	/* Compute hyperbolic distance between two points. */

	return std::acosh((std::cosh(zeta * radial) * std::cosh(zeta * radial_prime)) - (std::sinh(zeta * radial) * std::sinh(zeta * radial_prime) * cosAngle));
}

double getConnProb(double& hypDist, double& zeta, double& temp, double& mu) {
	/* Compute connection probability in the RHG model for two points with given hyperbolic distance. */

	return 1.0 / (1.0 + std::exp((zeta / (2.0 * temp)) * (hypDist - mu)));
}

bool parseInput(InputParser& input, long& idum, bool& exportCoordinates, int& nodes, int& dim, int& mode, \
	double& avgDeg, double& rescaledTemp, double& aCoeff, double& gamma, double& nu, double& rescaledRadius, \
	double& cutoff1, double& cutoff2, double& cutoff3, std::string& fname, \
	std::string& networkName, std::string& exportFile, std::string& metaFile) {
	/* Parse user input and perform input checks. */
	
	// File names
	std::string extension = ".dat";
	bool correctFileName = true;
	if (input.cmdOptionExists("-f")) {
		fname = input.getCmdOption("-f");
		if (!fname.empty()) {
			size_t dotIndex = fname.find_last_of(".");
			size_t extIndex = fname.find(extension);

			// Check if the last four characters are '.dat'
			if (dotIndex == fname.length() - extension.length() && extIndex != std::string::npos) {
				networkName = fname.substr(0, dotIndex);
				exportFile = networkName + ".coord" + extension;
				metaFile = networkName + ".meta" + extension;
			}
			else if (dotIndex == std::string::npos) {
				networkName = fname;
				fname = networkName + extension;
				exportFile = networkName + ".coord" + extension;
				metaFile = networkName + ".meta" + extension;
			}
			else {
				// File is not of .dat extension
				correctFileName = false;
			}
		}
		else {
			// Filename is empty
			correctFileName = false;
		}
	}
	else {
		// Filename not provided
		correctFileName = false;
	}
	if (!correctFileName) {
		std::cout << "Incorrect filename provided. Use the -f option to provide a filename with " << extension \
			<< " extension." << std::endl;
		return true;
	}

	/* Parameters */
	// Network size n
	bool nSpecified = true;
	if (input.cmdOptionExists("-n")) {
		std::string nodesString = input.getCmdOption("-n");
		if (!nodesString.empty()) {
			if (nodesString.find(".") == std::string::npos) {
				nodes = std::stoi(nodesString);
				if (nodes <= 1) {
					std::cout << "The network size n must be an integer > 1. " << \
						"Use the -n option to provide the network size n." << std::endl;
					return true;
				}
			}
			else {
				std::cout << "The network size n must be an integer > 1. " << \
					"Use the -n option to provide the network size n." << std::endl;
				return true;
			}
		}
		else {
			nSpecified = false;
		}
	}
	else {
		nSpecified = false;
	}
	if (!nSpecified) {
		std::cout << "The network size n must be specified. Use the -n option to provide the network size n." << std::endl;
		return true;
	}

	// Dimensionality d
	bool dSpecified = true;
	if (input.cmdOptionExists("-d")) {
		std::string dimString = input.getCmdOption("-d");
		if (!dimString.empty()) {
			if (dimString.find(".") == std::string::npos) {
				dim = std::stoi(dimString);
				if (dim < 1) {
					std::cout << "The dimensionality d must be an integer > 0. " << \
						"Use the -d option to provide a valid dimensionality d" << std::endl;
					return true;
				}
			}
			else {
				std::cout << "The dimensionality d must be an integer > 0. " << \
					"Use the -d option to provide a valid dimensionality d." << std::endl;
				return true;
			}
		}
		else {
			dSpecified = false;
		}
	}
	else {
		dSpecified = false;
	}
	if (!dSpecified) {
		std::cout << "The dimensionality d must be specified. Use the -d option to provide the dimensionality d." << std::endl;
		return true;
	}

	// Mode
	int nModes = 0;
	if (input.cmdOptionExists("-u")) {
		// User-based
		mode = 1;
		nModes++;
	}
	if (input.cmdOptionExists("-h")) {
		// Hybrid
		mode = 2;
		nModes++;
	}
	if (input.cmdOptionExists("-m")) {
		// Model-based
		mode = 3;
		nModes++;
	}
	if(nModes == 0) {
		std::cout << "No valid mode specified. Use one of the options: {-u, -h, -m}." << std::endl;
		return true;
	}
	else if (nModes == 1) {
		std::cout << "Mode: ";
		switch (mode) {
		case 1: std::cout << "user "; break;
		case 2: std::cout << "hybrid "; break;
		case 3: std::cout << "model-based "; break;
		}
		std::cout << std::endl;
	}
	else {
		std::cout << "Multiple modes selected. Use only one of the options: {-u, -h, -m}." << std::endl;
		return true;
	}

	// Mode dependent parameters
	if (mode == 1) {
		// User-based mode
		std::cout << "User-based mode not yet available." << std::endl;
		return true;
	}
	else if (mode == 2) {
		// Hybrid mode
		if (!((input.cmdOptionExists("-g") || input.cmdOptionExists("-a")) && input.cmdOptionExists("-k") && input.cmdOptionExists("-t"))) {
			std::cout << "Incorrect parameters provided. In hybrid mode, one must provide negative power-law exponent " \
				<< "gamma(-g) or rescaled radial component a(-a), average degree <k>(-k), and scaled temperature tau(-t).";
			return true;
		}

		// Check validity parameters hybrid mode
		if (input.cmdOptionExists("-g")) {
			gamma = std::stof(input.getCmdOption("-g"));
			if (gamma < 2) {
				std::cout << "The negative power-law exponent gamma must be >= 2. " \
					<< "Use the -g option to provide a valid negative power-law exponent gamma." << std::endl;
				return true;
			}
		}
		else {
			aCoeff = std::stof(input.getCmdOption("-a"));
			if (aCoeff < 1) {
				std::cout << "The rescaled radial component a must be >= 1. " \
					<< "Use the -a option to provide a valid rescaled radial component a." << std::endl;
				return true;
			}
		}
		
		avgDeg = std::stof(input.getCmdOption("-k"));
		rescaledTemp = std::stof(input.getCmdOption("-t"));
		if (rescaledTemp < 0) {
			std::cout << "The rescaled temperature tau must be >= 0. " \
				<< "Use the -t option to provide a valid rescaled temperature tau." << std::endl;
			return true;
		}
		if (avgDeg <= 0) {
			std::cout << "The average degree <k> must be > 0. " \
				<< "Use the -k option to provide a valid average degree <k>." << std::endl;
			return true;
		}

		// Check if gamma and tau are compatible
		if (input.cmdOptionExists("-g")) {
			if (rescaledTemp < 1.0 + cutoff2) {
				// tau <= 1
				aCoeff = gamma - 1.0;
			}
			else {
				if (gamma < rescaledTemp + 1.0) {
					// Invalid gamma and tau combination
					std::cout << "Chosen parameter gamma not possible at tau = " << rescaledTemp << ". " \
						<< "Use the -g and -t options to provide a valid combination of the parameter gamma " \
						<< "and rescaled temperature tau." << std::endl;
					return true;
				}
				else {
					// tau > 1
					aCoeff = (gamma - 1.0) / rescaledTemp;
				}
			}
		}
		else {
			// Compute gamma
			if (rescaledTemp < 1.0 + cutoff2) {
				// tau <= 1
				gamma = aCoeff + 1.0;
			}
			else {
				// tau > 1
				gamma = aCoeff * rescaledTemp + 1.0;
			}
		}

	}
	else {
		// Model-based mode
		if (!((input.cmdOptionExists("-g") || input.cmdOptionExists("-a")) && \
			((input.cmdOptionExists("-nu") || input.cmdOptionExists("-radius"))) && \
			input.cmdOptionExists("-t"))) {
			std::cout << "Incorrect parameters provided. In model-based mode, one must provide negative power-law exponent negative power-law exponent " \
				<< "gamma(-g) or rescaled radial component a(-a), scaling parameter nu(-nu) or rescaled radius(-radius), and scaled temperature tau(-t)." << std::endl;
			return true;
		}

		// Check validity parameters model-based mode
		rescaledTemp = std::stof(input.getCmdOption("-t"));
		if (rescaledTemp < 0) {
			std::cout << "The rescaled temperature tau must be >= 0. " \
				<< "Use the -t option to provide a valid rescaled temperature tau." << std::endl;
			return true;
		}
		if (input.cmdOptionExists("-g")) {
			gamma = std::stof(input.getCmdOption("-g"));
			if (gamma < 2) {
				std::cout << "The negative power-law exponent gamma must be >= 2. " \
					<< "Use the -g option to provide a valid negative power-law exponent gamma." << std::endl;
				return true;
			}
			if (rescaledTemp <= 1.0 + cutoff2) {
				// tau <= 1
				aCoeff = gamma - 1.0;
			}
			else {
				if (gamma < rescaledTemp + 1.0) {
					// Invalid gamma and tau combination
					std::cout << "Chosen parameter gamma not possible at tau = " << rescaledTemp << ". " \
						<< "Use the -g and -t options to provide a valid combination of the parameter gamma " \
						<< "and rescaled temperature tau: tau <= gamma - 1." << std::endl;
					return true;
				}
				else {
					// tau > 1
					aCoeff = (gamma - 1.0) / rescaledTemp;
				}
			}
		}
		else {
			aCoeff = std::stof(input.getCmdOption("-a"));
			if (aCoeff < 1) {
				std::cout << "The rescaled radial component a must be >= 1. " \
					<< "Use the -a option to provide a valid rescaled radial component a." << std::endl;
				return true;
			}

			// Compute gamma
			if (rescaledTemp < 1.0 + cutoff2) {
				// tau <= 1
				gamma = aCoeff + 1.0;
			}
			else {
				// tau > 1
				gamma = aCoeff * rescaledTemp + 1.0;
			}
			
		}

		std::cout << "Regime: ";
		if (input.cmdOptionExists("-nu")) {
			nu = std::stof(input.getCmdOption("-nu"));
			if (nu <= 0) {
				std::cout << "The scaling parameter nu must be > 0. " \
					<< "Use the -nu option to provide a valid scaling parameter nu." << std::endl;
				return true;
			}

			if (rescaledTemp > 1.0 + cutoff2) {
				// tau > 1
				rescaledRadius = rescaledTemp * std::log(nodes / nu);
				std::cout << "tau > 1, a >= 1" << std::endl;
			}
			else if (std::abs(rescaledTemp - 1.0) < cutoff2) {
				// tau = 1
				if (aCoeff > 1.0 + cutoff3) {
					// a > 1
					std::cout << "tau = 1, a > 1" << std::endl;
					if(nodes / nu < std::exp(1)) {
						// Lambert W not defined
						std::cout << "Chosen parameter nu not possible for network size n = " << nodes << " at tau = 1. " \
						<< "Use the -n and -nu options to provide a valid combination of the network size n " \
						<< "and parameter nu: n/nu >= e." << std::endl;
						return true;
					} if(nodes / nu <= 1e9) {
						// Lambert W value contained within look-up table
						rescaledRadius = lambert_lookup(nodes/nu);
						if(rescaledRadius < 0) {
							std::cout << "Chosen parameter nu not possible for network size n = " << nodes << " at tau = 1. " \
							<< "Use the -n and -nu options to provide a valid combination of the network size n " \
							<< "and parameter nu: n/nu >= e." << std::endl;
							return true;
						}
					}
					else {
						// Approximate Lambert W for extremely large network size n 
						rescaledRadius = std::log(nodes / nu) + std::log(std::log(nodes / nu));
						std::cout << "WARNING: Approximating Lambert W function for extremely large network size n." << std::endl;
					}
				} 
				else {
					// a = 1
					rescaledRadius = std::log(nodes / nu);
					std::cout << "tau = 1, a = 1" << std::endl;
				}
			}
			else {
				// 0 <= tau < 1
				rescaledRadius = std::log(nodes / nu);
				if (aCoeff > 1.0 + cutoff3) {
					// a > 1
					std::cout << "tau < 1, a > 1" << std::endl;
				} 
				else {
					// a = 1
					std::cout << "tau < 1, a = 1" << std::endl;
				}
			}
		}
		else {
			rescaledRadius = std::stof(input.getCmdOption("-radius"));
			if (rescaledRadius <= 0) {
				std::cout << "The rescaled radius R_H a must be > 0. " << \
					"Use the -radius option to provide a valid rescaled radius R_H." << std::endl;
				return true;
			}

			if (rescaledTemp > 1.0 + cutoff2) {
				// tau > 1
				nu = nodes * std::exp(-rescaledRadius / rescaledTemp);
				std::cout << "tau > 1, a >= 1" << std::endl;
			}
			else if (std::abs(rescaledTemp - 1.0) < cutoff2) {
				// tau = 1
				if (aCoeff > 1.0 + cutoff3) {
					// a > 1
					nu = nodes * rescaledRadius * std::exp(-rescaledRadius);
					std::cout << "tau = 1, a > 1" << std::endl;
				} 
				else {
					// a = 1
					nu = nodes * std::exp(-rescaledRadius);
					std::cout << "tau = 1, a = 1" << std::endl;
				}
			}
			else {
				// 0 <= tau < 1
				nu = nodes * std::exp(-rescaledRadius);
				if (aCoeff > 1.0 + cutoff3) {
					std::cout << "tau < 1, a > 1" << std::endl;
				} 
				else {
					std::cout << "tau < 1, a = 1" << std::endl;
				}
			}
		}
		std::cout << std::endl;
	}
	

	// Export coordinates
	if (input.cmdOptionExists("-v")) {
		exportCoordinates = true;
	}
	else {
		exportCoordinates = false;
	}

	/* Initialize random generator */
	std::random_device r;
	bool seedSpecified = true;
	if (input.cmdOptionExists("-seed")) {
		std::string seedString = input.getCmdOption("-seed");
		if (!seedString.empty()) {
			if (seedString.find(".") == std::string::npos) {
				idum = std::stol(seedString);
				if (idum == 0) {
					std::cout << "The random generator seed must be non-zero. " << \
						"Use the -seed option to provide a seed." << std::endl;
					return true;
				}
				else if (idum > 0) {
					// We need a negative seed
					idum = idum * ((long)-1);
				}
			}
			else {
				std::cout << "The random generator seed must be an integer. " << \
					"Use the -seed option to provide a seed." << std::endl;
				return true;
			}
		}
		else {
			seedSpecified = false;
		}
	}
	else {
		seedSpecified = false;
	}

	if (!seedSpecified) {
		//std::cout << "No seed for the pseudorandom generator provided. Proceeding with an arbitrarily chosen seed." << std::endl;
		idum = r() * ((long)-1);
	}

	return false;
}

int main(int argc, char* argv[]) {

	/* Variable declaration */
	long idum;
	int nodes, dim, mode;
	double avgDeg, temp, rescaledTemp, alpha, aCoeff, gamma, mu, nu, zeta, \
		radiusHyp, rescaledRadius, rescaledMu, constantId;
	std::string fname, networkName, exportFile, metaFile;
	bool exportCoordinates, degenerateModel = false;

	/* Numerical settings */
	double TOL = 1e-3;
	double cutoff1 = 0.05, cutoff2 = 0.01, cutoff3 = 0.01;
	double numLimit = 600.0; 
	bool debug = false;

	/* Parse input */
	InputParser input(argc, argv);
	bool inputError = parseInput(input, idum, exportCoordinates, nodes, dim, mode, avgDeg, rescaledTemp, aCoeff, \
		gamma, nu, rescaledRadius, cutoff1, cutoff2, cutoff3, fname, networkName, exportFile, metaFile);
	if (inputError) {
		// Abort program
		return 0;
	}


	/* Timers */
	clock_t t_coords_start, t_connections_start, t_overall_start;
	t_overall_start = clock();
	

	/* RHG parameters */
	zeta = 1.0 , alpha = (zeta / 2.0) * aCoeff, temp = rescaledTemp / dim;
	constantId = std::sqrt(pi) * (tgamma(dim / 2.0) / tgamma((dim + 1) / 2.0));  // Constant I_{d,1}
	

	/* Numerical solving for radius hyperbolic ball */
	if (mode == 2) {
		// Hybrid mode: we need to find the radius numerically

		double leftInit = 0.0, rightInit; 
		if (nodes > 1e9) {
			// TODO: do estimate based on log(n)
			// Search for rescaled radius in (0, 50.0)
			rightInit = 50.0;
		}
		else {
			// Search for rescaled radius in (0, 35.0)
			rightInit = 35.0;
		}

		// GSL objects
		const gsl_rng_type* T;
		gsl_rng* r;
		gsl_rng_env_setup();
		T = gsl_rng_default;
		r = gsl_rng_alloc(T);

		std::cout << "Finding radius hyperbolic ball... ";
		std::cout.flush();
		rescaledRadius = bisection(avgDeg, leftInit, rightInit, nodes, dim, rescaledTemp, aCoeff, constantId, \
			cutoff1, r, TOL, debug);
		if (rescaledRadius < 0) {
			std::cout << "Chosen network cannot be generated. " << \
				"Please choose a different configuration of parameters." << std::endl;
			return 0;
		}
		else {
			std::cout << "Done." << std::endl;
		}
		radiusHyp = (2.0 / (dim * zeta)) * rescaledRadius;

		std::cout << "Regime: ";
		if (rescaledTemp > 1.0 + cutoff2) {
			// tau > 1
			nu = nodes * std::exp(-rescaledRadius / rescaledTemp);
			std::cout << "tau > 1, a >= 1" << std::endl;
		}
		else if (std::abs(rescaledTemp - 1.0) < cutoff2) {
			// tau = 1
			if (aCoeff > 1.0 + cutoff3) {
				// a > 1
				nu = nodes * rescaledRadius * std::exp(-rescaledRadius);
				std::cout << "tau = 1, a > 1" << std::endl;
			} 
			else {
				// a = 1
				nu = nodes * std::exp(-rescaledRadius);
				std::cout << "tau = 1, a = 1" << std::endl;
			}
		}
		else {
			// 0 <= tau < 1
			nu = nodes * std::exp(-rescaledRadius);
			if (aCoeff > 1.0 + cutoff3) {
				// a > 1
				std::cout << "tau < 1, a > 1" << std::endl;
			} 
			else {
				// a = 1
				std::cout << "tau < 1, a = 1" << std::endl;
			}
		}
		std::cout << std::endl;

		// Clear GSL objects
		gsl_rng_free(r);
	}

	// Get rescaled mu
	if (rescaledTemp > 1.0 + cutoff2) {
		// tau > 1
		rescaledMu = rescaledRadius;
	}
	else if (std::abs(rescaledTemp - 1.0) < cutoff2) {
		// tau = 1
		if (aCoeff > 1.0 + cutoff3) {
			// a > 1
			rescaledMu = rescaledRadius;
		} 
		else {
			// a = 1
			rescaledMu = rescaledRadius - 2.0 * std::log(rescaledRadius);
		}
	}
	else {
		// 0 <= tau < 1
		if (aCoeff > 1.0 + cutoff3) {
			// a > 1
			rescaledMu = rescaledRadius;
		} 
		else {
			// a = 1
			rescaledMu = rescaledRadius - std::log(rescaledRadius);
		}
	}

	// Calculate model parameters from rescaled 
	radiusHyp = (2.0 / (dim * zeta)) * rescaledRadius;
	mu = (2.0 / (dim * zeta)) * rescaledMu; 


	/* Print network information */
	std::cout << "------ RHG settings ------" << std::endl;
	std::cout << "File: " << fname << std::endl;
	std::cout << "n: " << nodes << std::endl << "d: " << dim << std::endl << "a: " << aCoeff << std::endl \
		<< "tau: " << rescaledTemp << std::endl << "nu: " << nu << std::endl << "radius: " << radiusHyp << std::endl \
		<< "rescaled radius: " << rescaledRadius << std::endl;

	if(alpha * dim * radiusHyp > numLimit) {
		degenerateModel = true;
		std::cout << std::endl << "WARNING: value of parameter a extremely large wrt dimension d and network size n, " \
		<< "resorting to degenerate RHG model where all radial coordinates r = R." << std::endl;
	}


	/* Output file streams */
	std::ofstream ofile, expfile, mfile, logfile;

	/* Generate n points in a hyperbolic ball of dimensionality d + 1 */
	std::vector<std::vector<double>> theta(nodes, std::vector<double>(dim, 0.0));  // Each node has angular coordinates \theta_1, ..., \theta_d
	std::vector<std::vector<double>> sin_theta(nodes, std::vector<double>(dim, 0.0));  // Pre-calculated sine of angular coordinates
	std::vector<std::vector<double>> cos_theta(nodes, std::vector<double>(dim, 0.0));  // Pre-calculated cosine of angular coordinates
	std::vector<double> radial(nodes, 0.0); // Each node has a radial coordinate r
	double candidate;
	bool accepted;
	std::cout << std::endl << "Generating coordinates... ";
	std::cout.flush();
	t_coords_start = clock();
	// Generate radial coordinates
	if(degenerateModel) {
		for (int i = 0; i < nodes; i++) {
			radial[i] = radiusHyp;
		}
	} else {
		for (int i = 0; i < nodes; i++) {
			radial[i] = (1.0 / (alpha * dim)) * std::log(ran1(&idum) * (std::exp(alpha * dim * radiusHyp) - 1.0) + 1.0);
		}
	}
	for (int i = 0; i < nodes; i++) {

		// Angular coordinates dimensions 1 to d-1
		for (int d = 0; d < dim - 1; d++) {

			// Draw angular coordinate from [sin(theta)]^{d - k} with acceptance-rejection method
			accepted = false;
			while (!accepted) {
				// Draw candidate 
				candidate = acos(1.0 - 2.0 * ran1(&idum));

				if (ran1(&idum) <= pow(sin(candidate), (double)(dim - d - 2.0))) {
					// Accept candidate when unif random nr <= pdf
					theta[i][d] = candidate;
					sin_theta[i][d] = sin(theta[i][d]);
					cos_theta[i][d] = cos(theta[i][d]);
					accepted = true;
				}
			}
		}


		// Draw angular coordinate final dimension
		theta[i][dim - 1] = ran1(&idum) * 2.0 * pi;
		sin_theta[i][dim - 1] = sin(theta[i][dim - 1]);
		cos_theta[i][dim - 1] = cos(theta[i][dim - 1]);

	}
	float time_coords = (float)(clock() - t_coords_start) / CLOCKS_PER_SEC;


	/* Write coordinates */
	if (exportCoordinates) {
		expfile.open(exportFile, std::ios_base::out | std::ios_base::trunc);
		expfile << "id radial";
		for (int d = 0; d < dim; d++) {
			expfile << " theta" << d + 1;
		}
		expfile << std::endl;

		// Write coordinates for all nodes
		for (int i = 0; i < nodes; i++) {
			expfile << i << ' ' << radial[i];

			// Write coordinates of all dimensions
			for (int d = 0; d < dim; d++) {
				expfile << ' ' << theta[i][d];
			}

			expfile << std::endl;
		}
		expfile.close();
	}
	std::cout << "Done." << std::endl;
	
	/* Write meta info */
	mfile.open(metaFile, std::ios_base::out | std::ios_base::trunc);
	mfile << "name: " << networkName << std::endl;
	mfile << "nodes: " << nodes << std::endl << "dim: " << dim << std::endl << "alpha: " << alpha << std::endl \
		<< "a: " << aCoeff << std::endl << "gamma: " << gamma << std::endl << "temp: " << temp << std::endl \
		<< "tau: " << rescaledTemp << std::endl;
	mfile << "nu: " << nu << std::endl << "radius: " << radiusHyp  << std::endl  << "scaled radius: " << rescaledRadius << std::endl;
	mfile.close();

	/* Simulate connections */
	double hypDist, connectionProb, cosAngle;
	std::cout << "Generating links... ";
	std::cout.flush();
	t_connections_start = clock();
	ofile.open(fname, std::ios_base::out | std::ios_base::trunc);
	if (rescaledTemp < cutoff1) {
		// At tau = 0 we have the step model
		for (int i = 0; i < nodes - 1; i++) {
			for (int j = i + 1; j < nodes; j++) {

				// Calculate angle between i and j
				cosAngle = getCosAngle(dim, sin_theta[i], sin_theta[j], cos_theta[i], cos_theta[j]);

				// Calculate hyperbolic distance between i and j
				hypDist = getHypDist(zeta, radial[i], radial[j], cosAngle);

				// Connect nodes in the step model if distance < radius
				if (hypDist < radiusHyp) {
					ofile << i << ' ' << j << std::endl;
				}
			}
		}
	}
	else {
		// Tau > 0
		for (int i = 0; i < nodes - 1; i++) {
			for (int j = i + 1; j < nodes; j++) {

				// Calculate angle between i and j
				cosAngle = getCosAngle(dim, sin_theta[i], sin_theta[j], cos_theta[i], cos_theta[j]);

				// Calculate hyperbolic distance between i and j
				hypDist = getHypDist(zeta, radial[i], radial[j], cosAngle);

				// Calculate connection probability nodes i and j
				connectionProb = getConnProb(hypDist, zeta, temp, mu);

				// Simulate connection probability and assign link accordingly
				if (ran1(&idum) < connectionProb) {
					ofile << i << ' ' << j << std::endl;
				}
			}
		}
	}
	ofile.close();
	float time_connections = (float)(clock() - t_connections_start) / CLOCKS_PER_SEC;
	float time_execution = (float)(clock() - t_overall_start) / CLOCKS_PER_SEC;
	std::cout << "Done." << std::endl << std::endl;

	/* Print execution time */
	std::cout << "------ Execution times ------" << std::endl;
	std::cout << "Finding radius: " << time_execution - time_coords - time_connections << std::endl;
	std::cout << "Generating coordinates: " << time_coords << std::endl;
	std::cout << "Generating links: " << time_connections << std::endl;
	std::cout << "Overall: " << time_execution << std::endl;

	return 0;
}