#include <iostream>
#include <fstream>
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
using namespace std;

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

	const string& getCmdOption(const string& option) const {
		/* Get cmd option value. */
		vector<string>::const_iterator itr;
		itr = find(this->tokens.begin(), this->tokens.end(), option);
		if (itr != this->tokens.end() && ++itr != this->tokens.end()) {
			return *itr;
		}
		static const string empty_string("");
		return empty_string;
	}

	bool cmdOptionExists(const string& option) const {
		/* Check if cmd option exists. */
		return std::find(this->tokens.begin(), this->tokens.end(), option)
			!= this->tokens.end();
	}
private:
	vector<string> tokens;
};

template <typename T> int sgn(T val) {
	/* Sign function. */

	return (T(0) < val) - (val < T(0));
}

template <typename T> int heaviside(T val) {
	/* Heaviside step function. */

	return (val > T(0));
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
	return factor * exp(params->aCoeff * (x[0] - params->radius)) * exp(params->aCoeff * (x[1] - params->radius)) * \
		(pow(sin(x[2]), params->d - 1.0) / (1.0 + exp((x[0] + x[1] - params->radius) / params->tau) * \
			pow(sin(x[2] / 2.0), params->d / params->tau)));
}

double avg_k_tzero(double* x, size_t dim, void* passed_parameters) {
	/* Integrand for average degree <k> integral at tau = 0. */
	(void)(dim); // avoid unused parameter warnings 

	struct int_params* params = (struct int_params*)passed_parameters;

	double factor = (params->nodes - 1.0) * (pow(params->aCoeff, 2.0) / params->constantId);
	return factor * exp(params->aCoeff * (x[0] - params->radius)) * exp(params->aCoeff * (x[1] - params->radius)) * \
		pow(sin(x[2]), params->d - 1.0) * heaviside(params->radius - x[0] - x[1] - (params->d * log(sin(x[2]/2.0))));
}

double integralAvgDegMISER(const double& radius, const int& nodes, const int& dim, const double& rescaledTemp, \
	const double& aCoeff, const double& constantId, gsl_rng* r, const bool& debug) {
	/* Integral computation for desired average degree <k> with the MISER algorithm Monte Carlo integration. */

	if (radius == 0) {   
		return 0.0;
	}

	// Integration limits
	double xl[3] = { 0, 0, 0 };
	double xu[3] = { radius, radius, pi };

	// GSL integration objects
	gsl_monte_function integrand;
	if (rescaledTemp > 0) {
		// For tau > 0 we integrate over the regular connection probability function
		integrand = { &avg_k, 3, 0 };
	}
	else {
		// For tau = 0 we integrate over the step function
		integrand = { &avg_k_tzero, 3, 0 };
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
	const double& aCoeff, const double& constantId, gsl_rng* r, const bool& debug) {
	/* Integral computation for desired average degree <k> with importance sampling Monte Carlo integration. */

	if (radius == 0) {
		return 0.0;
	}

	// Integration limits
	double xl[3] = { 0, 0, 0 };
	double xu[3] = { radius, radius, pi };

	// GSL integration objects
	gsl_monte_function integrand;
	if (rescaledTemp > 0) {
		// For tau > 0 we integrate over the regular connection probability function
		integrand = { &avg_k, 3, 0 };
	}
	else {
		// For tau = 0 we integrate over the step function
		integrand = { &avg_k_tzero, 3, 0 };
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
	const double& rescaledTemp, const double& aCoeff, const double& constantId, gsl_rng* r, const bool& debug) {
	/* Target function of integral computation for desired average degree <k> with the MISER algorithm Monte Carlo integration. */

	return integralAvgDegMISER(radius, nodes, dim, rescaledTemp, aCoeff, constantId, r, debug) - avgDeg;
}

double targetAvgDeg(const double& radius, const double& avgDeg, const int& nodes, const int& dim, \
	const double& rescaledTemp, const double& aCoeff, const double& constantId, gsl_rng* r, const bool& debug) {
	/* Target function computation integral for desired average degree <k> with importance sampling Monte Carlo integration. */

	return integralAvgDeg(radius, nodes, dim, rescaledTemp, aCoeff, constantId, r, debug) - avgDeg;
}


double bisection(const double& avgDeg, const double& leftInit, const double& rightInit, const int& nodes, \
	const int& dim, const double& rescaledTemp, const double& aCoeff, const double& constantId, gsl_rng* r, \
	const double& TOL, const bool& debug) {
	/* Compute rescaled radius numerically with bisection method for desired average degree <k>. */

	double left = leftInit, right = rightInit;
	double mid = (left + right) / 2.0;

	
	// Compute initial values with the MISER routine
	double fLeft = targetAvgDegMISER(left, avgDeg, nodes, dim, rescaledTemp, aCoeff, constantId, r, debug);
	double fRight = targetAvgDegMISER(right, avgDeg, nodes, dim, rescaledTemp, aCoeff, constantId, r, debug);
	double fMid = targetAvgDegMISER(mid, avgDeg, nodes, dim, rescaledTemp, aCoeff, constantId, r, debug);

	// Find a valid left starting point with different sign from right function value
	while (sgn(fRight) == sgn(fLeft)) {
		left += 0.1;
		fLeft = targetAvgDegMISER(left, avgDeg, nodes, dim, rescaledTemp, aCoeff, constantId, r, false);
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
			fMid = targetAvgDeg(mid, avgDeg, nodes, dim, rescaledTemp, aCoeff, constantId, r, debug);
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

	return (sqrt(pi) * (tgamma((dim - k + 1) / 2.0) / tgamma(1 + (dim - k) / 2.0)));
}

double getHypDist(const int &dim, const double& zeta, const vector<double>& hidden, \
	const vector<double>& hidden_prime) {
	/* Compute hyperbolic distance between two points. Approximate trig functions inside acosh(). */

	// Obtain radial distances node pair
	double radial = hidden[0];
	double radial_prime = hidden_prime[0];

	// Calculate angle between the two points across the d dimensions
	double sineProduct = 1;
	double cosAngle = 0;
	for (int d = 0; d < dim; d++) {
		cosAngle += sineProduct * cos(hidden[d + 1]) * cos(hidden_prime[d + 1]);
		sineProduct *= sin(hidden[d + 1]) * sin(hidden_prime[d + 1]);
	}
	cosAngle += sineProduct;

	// Calculate hyperbolic distance
	//double coshDist = (cosh(zeta * radial) * cosh(zeta * radial_prime)) - (sinh(zeta * radial) * sinh(zeta * radial_prime) * cosAngle);
	double hypDist = (1.0 / zeta) * acosh((cosh(zeta * radial) * cosh(zeta * radial_prime)) - (sinh(zeta * radial) * sinh(zeta * radial_prime) * cosAngle));
	//double hypDist = acosh(exp(zeta * (radial + radial_prime)) * (1.0 - cosAngle) / 4.0) / zeta;  // Approximate hyperbolic trig functions inside acosh()

	return hypDist;
}

double getConnProb(const double& hypDist, const double& zeta, const double& temp, const double& mu) {
	/* Compute connection probability in the RHG model for two points with given hyperbolic distance. */

	return 1.0 / (1.0 + exp((zeta / (2.0 * temp)) * (hypDist - mu)));
}

bool parseInput(InputParser& input, long& idum, bool& exportCoordinates, int& nodes, int& dim, int& mode, \
	double& avgDeg, double& rescaledTemp, double& aCoeff, \
	double& gamma, double& nu, double& rescaledRadius, double& radiusHyp, double& mu, string& fname, \
	string& networkName, string& exportFile, string& metaFile) {
	/* Parse user input and perform input checks. */
	
	// File names
	string extension = ".dat";
	bool correctFileName = true;
	if (input.cmdOptionExists("-f")) {
		fname = input.getCmdOption("-f");
		if (!fname.empty()) {
			size_t dotIndex = fname.find_last_of(".");
			size_t extIndex = fname.find(extension);

			// Check if the last four characters are '.dat'
			if (dotIndex == fname.length() - extension.length() && extIndex != string::npos) {
				networkName = fname.substr(0, dotIndex);
				exportFile = networkName + ".coord" + extension;
				metaFile = networkName + ".meta" + extension;
			}
			else if (dotIndex == string::npos) {
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
		cout << "Incorrect filename provided. Use the -f option to provide a filename with " << extension \
			<< " extension." << endl;
		return true;
	}

	/* Parameters */
	// Network size n
	bool nSpecified = true;
	if (input.cmdOptionExists("-n")) {
		string nodesString = input.getCmdOption("-n");
		if (!nodesString.empty()) {
			if (nodesString.find(".") == string::npos) {
				nodes = stoi(nodesString);
				if (nodes <= 1) {
					cout << "The network size n must be an integer > 1. " << \
						"Use the -n option to provide the network size n." << endl;
					return true;
				}
			}
			else {
				cout << "The network size n must be an integer > 1. " << \
					"Use the -n option to provide the network size n." << endl;
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
		cout << "The network size n must be specified. Use the -n option to provide the network size n." << endl;
		return true;
	}

	// Dimensionality d
	bool dSpecified = true;
	if (input.cmdOptionExists("-d")) {
		string dimString = input.getCmdOption("-d");
		if (!dimString.empty()) {
			if (dimString.find(".") == string::npos) {
				dim = stoi(dimString);
				if (dim < 1) {
					cout << "The dimensionality d must be an integer > 0. " << \
						"Use the -d option to provide a valid dimensionality d" << endl;
					return true;
				}
			}
			else {
				cout << "The dimensionality d must be an integer > 0. " << \
					"Use the -d option to provide a valid dimensionality d." << endl;
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
		cout << "The dimensionality d must be specified. Use the -d option to provide the dimensionality d." << endl;
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
		cout << "No valid mode specified. Use one of the options: {-u, -h, -m}." << endl;
		return true;
	}
	else if (nModes == 1) {
		cout << "Mode: ";
		switch (mode) {
		case 1: cout << "user. "; break;
		case 2: cout << "hybrid. "; break;
		case 3: cout << "model-based. "; break;
		}
		cout << endl << endl;
	}
	else {
		cout << "Multiple modes selected. Use only one of the options: {-u, -h, -m}." << endl;
		return true;
	}

	// Mode dependent parameters
	if (mode == 1) {
		// User-based mode
		cout << "User-based mode not yet available." << endl;
		return true;
	}
	else if (mode == 2) {
		// Hybrid mode
		if (!((input.cmdOptionExists("-g") || input.cmdOptionExists("-a")) && input.cmdOptionExists("-k") && input.cmdOptionExists("-t"))) {
			cout << "Incorrect parameters provided. In hybrid mode, one must provide negative power-law exponent " \
				<< "gamma(-g) or rescaled radial component a(-a), average degree <k>(-k), and scaled temperature tau(-t).";
			return true;
		}

		// Check validity parameters hybrid mode
		if (input.cmdOptionExists("-g")) {
			gamma = stof(input.getCmdOption("-g"));
			if (gamma < 2) {
				cout << "The negative power-law exponent gamma must be >= 2. " \
					<< "Use the -g option to provide a valid negative power-law exponent gamma." << endl;
				return true;
			}
		}
		else {
			aCoeff = stof(input.getCmdOption("-a"));
			if (aCoeff < 1) {
				cout << "The rescaled radial component a must be >= 1. " \
					<< "Use the -a option to provide a valid rescaled radial component a." << endl;
				return true;
			}
		}
		
		avgDeg = stof(input.getCmdOption("-k"));
		rescaledTemp = stof(input.getCmdOption("-t"));
		if (rescaledTemp < 0) {
			cout << "The rescaled temperature tau must be >= 0. " \
				<< "Use the -t option to provide a valid rescaled temperature tau." << endl;
			return true;
		}
		if (avgDeg <= 0) {
			cout << "The average degree <k> must be > 0. " \
				<< "Use the -k option to provide a valid average degree <k>." << endl;
			return true;
		}

		// Check if gamma and tau are compatible
		if (input.cmdOptionExists("-g")) {
			if (rescaledTemp <= 1) {
				// tau <= 1
				aCoeff = gamma - 1.0;
			}
			else {
				if (gamma < rescaledTemp + 1) {
					// Invalid gamma and tau combination
					cout << "Chosen parameter gamma not possible at tau = " << rescaledTemp << ". " \
						<< "Use the -g and -t options to provide a valid combination of the parameter gamma " \
						<< "and rescaled temperature tau." << endl;
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
			if (rescaledTemp <= 1) {
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
			cout << "Incorrect parameters provided. In model-based mode, one must provide negative power-law exponent negative power-law exponent " \
				<< "gamma(-g) or rescaled radial component a(-a), scaling parameter nu(-nu) or rescaled radius(-radius), and scaled temperature tau(-t)." << endl;
			return true;
		}

		// Check validity parameters model-based mode
		rescaledTemp = stof(input.getCmdOption("-t"));
		if (rescaledTemp < 0) {
			cout << "The rescaled temperature tau must be >= 0. " \
				<< "Use the -t option to provide a valid rescaled temperature tau." << endl;
			return true;
		}
		if (input.cmdOptionExists("-g")) {
			gamma = stof(input.getCmdOption("-g"));
			if (gamma < 2) {
				cout << "The negative power-law exponent gamma must be >= 2. " \
					<< "Use the -g option to provide a valid negative power-law exponent gamma." << endl;
				return true;
			}
			if (rescaledTemp <= 1) {
				// tau <= 1
				aCoeff = gamma - 1.0;
			}
			else {
				if (gamma < rescaledTemp + 1) {
					// Invalid gamma and tau combination
					cout << "Chosen parameter gamma not possible at tau = " << rescaledTemp << ". " \
						<< "Use the -g and -t options to provide a valid combination of the parameter gamma " \
						<< "and rescaled temperature tau." << endl;
					return true;
				}
				else {
					// tau > 1
					aCoeff = (gamma - 1.0) / rescaledTemp;
				}
			}
		}
		else {
			aCoeff = stof(input.getCmdOption("-a"));
			if (aCoeff < 1) {
				cout << "The rescaled radial component a must be >= 1. " \
					<< "Use the -a option to provide a valid rescaled radial component a." << endl;
				return true;
			}

			// Compute gamma
			if (rescaledTemp <= 1) {
				// tau <= 1
				gamma = aCoeff + 1.0;
			}
			else {
				// tau > 1
				gamma = aCoeff * rescaledTemp + 1.0;
			}
			
		}

		if (input.cmdOptionExists("-nu")) {
			nu = stof(input.getCmdOption("-nu"));
			if (nu <= 0) {
				cout << "The scaling parameter nu must be > 0. " \
					<< "Use the -nu option to provide a valid scaling parameter nu." << endl;
				return true;
			}
			
			if (rescaledTemp <= 1) {
				// tau <= 1
				rescaledRadius = log(nodes / nu);
			}
			else {
				// tau > 1
				rescaledRadius = rescaledTemp * log(nodes / nu);
			}
			
		}
		else {
			rescaledRadius = stof(input.getCmdOption("-radius"));
			if (rescaledRadius <= 0) {
				cout << "The rescaled radius R_H a must be > 0. " \
					<< "Use the -radius option to provide a valid rescaled radius R_H." << endl;
				return true;
			}

			if (rescaledTemp <= 1) {
				// tau <= 1
				nu = nodes / exp(rescaledRadius);
			}
			else {
				// tau > 1
				nu = nodes / exp(rescaledRadius / rescaledTemp);
			}
		}
	}

	// Export coordinates
	if (input.cmdOptionExists("-v")) {
		exportCoordinates = true;
	}
	else {
		exportCoordinates = false;
	}

	/* Initialize random generator */
	random_device r;
	bool seedSpecified = true;
	if (input.cmdOptionExists("-seed")) {
		string seedString = input.getCmdOption("-seed");
		if (!seedString.empty()) {
			if (seedString.find(".") == string::npos) {
				idum = stol(seedString);
				if (idum == 0) {
					cout << "The random generator seed must be non-zero. " << \
						"Use the -seed option to provide a seed." << endl;
					return true;
				}
				else if (idum > 0) {
					// We need a negative seed
					idum = idum * ((long)-1);
				}
			}
			else {
				cout << "The random generator seed must be an integer. " << \
					"Use the -seed option to provide a seed." << endl;
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
		//cout << "No seed for the pseudorandom generator provided. Proceeding with an arbitrarily chosen seed." << endl;
		idum = r() * ((long)-1);
	}

	return false;
}

int main(int argc, char* argv[]) {

	/* Variable declaration */
	long idum;
	int nodes, dim, mode;
	double avgDeg, temp, rescaledTemp, alpha, aCoeff, gamma, mu, nu, zeta, \
		radiusHyp, rescaledRadius, constantId;
	string fname, networkName, exportFile, metaFile;
	bool exportCoordinates;

	/* Parse input */
	InputParser input(argc, argv);
	bool inputError = parseInput(input, idum, exportCoordinates, nodes, dim, mode, avgDeg, rescaledTemp, aCoeff, \
		gamma, nu, rescaledRadius, radiusHyp, mu, fname, networkName, exportFile, metaFile);
	if (inputError) {
		// Abort program
		return 0;
	}

	/* Timers */
	clock_t t_coords_start, t_connections_start, t_overall_start;
	t_overall_start = clock();
	

	/* RHG parameters */
	zeta = 1.0 , alpha = (zeta / 2.0) * aCoeff, temp = rescaledTemp / dim;
	constantId = sqrt(pi) * (tgamma(dim / 2.0) / tgamma((dim + 1) / 2.0));  // Constant I_{d,1}


	/* Numerical settings */
	double leftInit = 0.0, rightInit; 
	if (nodes > 1e9) {
		// Search for rescaled radius in (0, 50.0)
		rightInit = 50.0;
	}
	else {
		// Search for rescaled radius in (0, 40.0)
		rightInit = 40.0;
	}
	double TOL = 1e-3;
	bool debug = false;

	/* Numerical solving for radius hyperbolic ball */
	if (mode == 2) {
		// Hybrid mode: we need to find the radius numerically

		// GSL objects
		const gsl_rng_type* T;
		gsl_rng* r;
		gsl_rng_env_setup();
		T = gsl_rng_default;
		r = gsl_rng_alloc(T);

		cout << "Finding radius hyperbolic ball... ";
		cout.flush();
		rescaledRadius = bisection(avgDeg, leftInit, rightInit, nodes, dim, rescaledTemp, aCoeff, constantId, \
			r, TOL, debug);
		if (rescaledRadius < 0) {
			cout << "Chosen network cannot be generated. " << \
				"Please choose a different configuration of parameters." << endl;
			return 0;
		}
		else {
			cout << "Done." << endl << endl;
		}
		radiusHyp = (2.0 / (dim * zeta)) * rescaledRadius;
		mu = radiusHyp;

		if (rescaledTemp > 1) {
			// tau > 1
			nu = nodes / exp(rescaledRadius / rescaledTemp);
		}
		else {
			// 0 <= tau <= 1
			nu = nodes / exp(rescaledRadius);
		}

		// Clear GSL objects
		gsl_rng_free(r);
	}
	else {
		// Model-based mode: parameter nu or radius provided by user
		if (rescaledTemp > 1) {
			// tau > 1
			rescaledRadius = rescaledTemp * log(nodes/nu);
		}
		else {
			// 0 <= tau <= 1
			rescaledRadius = log(nodes / nu);
		}
		radiusHyp = (2.0 / (dim * zeta)) * rescaledRadius;
		mu = radiusHyp;
	}


	/* Print network information */
	cout << "------ RHG settings ------" << endl;
	cout << "File: " << fname << endl;
	cout << "n: " << nodes << endl << "d: " << dim << endl << "a: " << aCoeff << endl \
		<< "tau: " << rescaledTemp << endl << "nu: " << nu << endl << "radius: " << radiusHyp << endl \
		<< "rescaled radius: " << rescaledRadius << endl;


	/* Output file streams */
	ofstream ofile, expfile, mfile, logfile;

	/* Generate n points in a hyperbolic ball of dimensionality d + 1 */
	vector<vector<double>> hidden(nodes, vector<double>(dim + 1, 0));  // Each node has hidden variables r, \theta_1, ..., \theta_d
	double candidate;
	bool accepted;
	cout << endl << "Generating coordinates... ";
	cout.flush();
	t_coords_start = clock();
	for (int i = 0; i < nodes; i++) {

		// Radial coordinate
		hidden[i][0] = (1.0 / (alpha * dim)) * log(ran1(&idum) * (exp(alpha * dim * radiusHyp) - 1.0) + 1.0);

		// Angular coordinates dimensions 1 to d-1
		for (int d = 0; d < dim - 1; d++) {

			// Draw angular coordinate from [sin(theta)]^{d - k} with acceptance-rejection method
			accepted = false;
			while (!accepted) {
				// Draw candidate 
				candidate = acos(1.0 - 2.0 * ran1(&idum));

				if (ran1(&idum) <= pow(sin(candidate), (double)(dim - d - 2))) {
					// Accept candidate when unif random nr <= pdf

					hidden[i][d + 1] = candidate;
					accepted = true;
				}
			}
		}


		// Draw angular coordinate final dimension
		hidden[i][dim] = ran1(&idum) * 2.0 * pi;

	}
	float time_coords = (float)(clock() - t_coords_start) / CLOCKS_PER_SEC;


	/* Write coordinates */
	if (exportCoordinates) {
		expfile.open(exportFile, ios_base::out | ios_base::trunc);
		expfile << "id radial";
		for (int d = 0; d < dim; d++) {
			expfile << " theta" << d + 1;
		}
		expfile << endl;

		// Write coordinates for all nodes
		for (int i = 0; i < nodes; i++) {
			expfile << i;

			// Write coordinates of all dimensions
			for (int d = 0; d < dim + 1; d++) {
				expfile << ' ' << hidden[i][d];
			}

			expfile << endl;
		}
		expfile.close();
	}
	cout << "Done." << endl;
	
	/* Write meta info */
	mfile.open(metaFile, ios_base::out | ios_base::trunc);
	mfile << "name: " << networkName << endl;
	mfile << "nodes: " << nodes << endl << "dim: " << dim << endl << "alpha: " << alpha << endl \
		<< "a: " << aCoeff << endl << "gamma: " << gamma << endl << "temp: " << temp << endl \
		<< "tau: " << rescaledTemp << endl;
	mfile << "nu: " << nu << endl << "radius: " << radiusHyp  << endl  << "scaled radius: " << rescaledRadius << endl;
	mfile.close();

	/* Simulate connections */
	double hypDist, connectionProb;
	cout << "Generating links... ";
	cout.flush();
	t_connections_start = clock();
	ofile.open(fname, ios_base::out | ios_base::trunc);
	if (rescaledTemp == 0.0) {
		// At tau = 0 we have the step model
		for (int i = 0; i < nodes - 1; i++) {
			for (int j = i + 1; j < nodes; j++) {

				// Calculate hyperbolic distance between i and j
				hypDist = getHypDist(dim, zeta, hidden[i], hidden[j]);

				// Connect nodes in the step model if distance < radius
				if (hypDist < radiusHyp) {
					ofile << i << ' ' << j << endl;
				}
			}
		}
	}
	else {
		// Tau > 0
		for (int i = 0; i < nodes - 1; i++) {
			for (int j = i + 1; j < nodes; j++) {

				// Calculate hyperbolic distance between i and j
				hypDist = getHypDist(dim, zeta, hidden[i], hidden[j]);

				// Calculate connection probability nodes i and j
				connectionProb = getConnProb(hypDist, zeta, temp, mu); 

				// Simulate connection probability and assign link accordingly
				if (ran1(&idum) < connectionProb) {
					ofile << i << ' ' << j << endl;
				}
			}
		}
	}
	ofile.close();
	float time_connections = (float)(clock() - t_connections_start) / CLOCKS_PER_SEC;
	float time_execution = (float)(clock() - t_overall_start) / CLOCKS_PER_SEC;
	cout << "Done." << endl << endl;

	/* Print execution time */
	cout << "------ Execution times ------" << endl;
	cout << "Finding radius: " << time_execution - time_coords - time_connections << endl;
	cout << "Generating coordinates: " << time_coords << endl;
	cout << "Generating links: " << time_connections << endl;
	cout << "Overall: " << time_execution << endl;
	
	
	// Code for logging the generation time
	//logfile.open("time_log.dat", ios::out | ios::app);
	//logfile << networkName << ' ' << nodes << ' ' << dim << ' ' << time_execution << endl;
	//logfile.close();

	return 0;
}