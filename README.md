# RHG Generator #

Network generator for random hyperbolic graphs (RHGs) with a latent representation in a hyperbolic ball of any dimension *d*+1 >= 2. The generator returns the network in edgelist format and additionally a list of the node coordinates if desired. The generator can be operated in 2 different modes: hybrid and model-based. In hybrid mode, the user provides the desired average degree, temperature and power-law exponent of the resulting network. In model-based mode, the user provides all network parameters. See *Required Parameters* for the full lists of the required parameters per mode.

### Dependencies
- C++ 2011 compiler or newer (`std=c++11`)
- Library `gsl`, GNU Scientific Library (version >= `1.16`).

### Setup
- Install the required libraries using a package manager of choice.
- Clone the repository.
``` 
git clone https://gbudel@bitbucket.org/gbudel/rhg-generator.git 
```
- Navigate to the directory `rhg-generator` and invoke the `make` command.
```
cd rhg-generator
make
```
- If compiling fails, please check if the required compiler and libraries are present and can be found by the compiler.

### Usage
Call the executable `generate_rhg` with the mode of choice (model-based or hybrid) and provide the required parameters. Parameter flags and values should be separated by a whitespace, e.g., `-n 1000` sets the number of nodes to 1,000. A full list of the available parameters and which parameters are required in which mode can found below. The inputs are tested for the restriction within parentheses.

#### Parameters
* `-f` 		filename (either a `*.dat` filename or a filename without extension).
* `-n` 		network size (`int > 1`).
* `-d` 		dimensionality parameter *d* for the dimensionality *d + 1* of the hyperbolic ball (`int >= 1`).
* `-m` 		model-based mode (select one from `{-m, -h}`.
* `-h`		hybrid-mode (select one from `{-m, -h}`.
* `-t` 		rescaled temperature *tau* (`float >= 0`).
* `-g` 		negative power-law exponent *gamma* for *P(k) ~ k^(-gamma)* (`float >= 2`).
* `-a` 		radial component *a* (`float >= 1`).
* `-nu`		scaling parameter *nu* (`float > 0`).
* `-radius`	rescaled radius *R* of the hyperbolic ball (`float > 0`).
* `-v` 		(optional) a switch whether or not to export node coordinates after generation (`<filename>.coord.dat`).
* `-seed` 	(optional) the seed for the pseudorandom generator (`long != 0`).

#### Required Parameters
* Always required: `-f AND -n AND -d`.
* `-m` model-based mode: `-t AND [-nu OR -radius] AND [-g OR -a]`.
* `-h` hybrid mode: `-k AND -t AND [-g OR -a]`.

### Output
The following files will be created in the project directory.

* `<filename>.dat` 		the network in edge list format with node indices separated by a whitespace.
* `<filename>.meta.dat` 	a meta file with the chosen/computed network parameters.
* `<filename>.coord.dat`	(optional) a list of the node coordinates in the hyperbolic ball.

### Examples
Hybrid mode: generate a network with *n = 10,000* nodes for dimensionality *d = 3*, average degree *<k> = 10*, power-law exponent *gamma = 2.1* in the cold regime (*tau = 1/2*).
```
./generate_rhg -f example.dat -n 10000 -d 3 -h -g 2.1 -k 10 -t 0.5
```

Model-based mode: generate a network with *n = 1,000* nodes for dimensionality *d = 1*, radial component *a = 1*, scaling parameter *nu = 0.5* in the hot regime (*tau = 3/2*) and export the coordinates with the `-v` switch afterwards.
```
./generate_rhg -f example.dat -n 1000 -d 1 -m -a 1.0 -nu 0.5 -t 1.5 -v
```

### Numerical considerations
To prevent numerical issues within the computations, the following cut-off values in the parameters are applied based on the user inputs:
* When *tau < 0.05*, links are generated deterministically as if *tau = 0* (and when *tau < 0*, an error is thrown).
* When *|tau - 1| < 0.01*, the scaling of the critical regime *tau = 1* is invoked.
* When *tau >= 1.01*, the scaling of the hot regime *tau > 1* is invoked, and the parameters must adhere to *gamma <= tau + 1*, otherwise an error is thrown (because the parameters imply *a < 1*).
* When *a < 1.01*, the scaling of the case *a = 1* is invoked (and when *a < 1*, an error is thrown).
* When *a > 2400/R*, all radial coordinates will be very large and they cannot be reliably generated, therefore we resort to degenerate hyperbolic graphs where each radial coordinate *r = R*.
* In the critical regime *tau = 1*, *a > 1*, the Lambert function *W(x)* is evaluated using a look-up table of function values and linear intrapolation, closely resembling the true function values. 
* If *n <= 10^7*, trigonometric function values *sin(theta)* and *cos(theta)* are pre-calculated and stored in memory for all nodes in order to speed up calculations. If *n > 10^7*, these function values are not pre-calculated in order to save memory.
* In the current implementation, the maximum network size is *n = 10^9*. Note that for large dimensionality *d*, the memory can already get full for network sizes *n < 10^9*. The software contains memory exception handlers for this purpose and will throw an error in case the machine runs out of memory.