# RHG Generator #

Network generator for random hyperbolic graphs (RHGs) with a latent representation in a hyperbolic ball of *d+1* dimensions. The generator returns the network in edgelist format and additionally a list of the node coordinates if desired. The generator can be operator in 3 different modes, see *Required Parameters* for lists of the required parameters per mode.

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
Call the executable `generate_rhg` with the mode of choice (user-based, hybrid or model-based) and provide the required parameters. Parameter flags and values should be seperated by a whitespace, e.g., `-n 1000` sets the number of nodes to 1000. A full list of the available parameters and which parameters are required in which mode can found below.

#### Parameters
* `-f` 		filename (either a `*.dat` filename or a filename without extension).
* `-n` 		network size (integer > 1).
* `-d` 		dimensionality d of the hyperbolic ball with dimensionality d + 1 (integer >= 1).
* `-u`		user-based mode (select one from `{-u, -h, -m}`.
* `-h`		hybrid-mode (select one from `{-u, -h, -m}`.
* `-m` 		model-based mode (select one from `{-u, -h, -m}`.
* `-k`		average degree <k> (float > 0). 
* `-c`		clustering coefficient c (0 <= float <= 1). 
* `-t` 		rescaled temperature tau (float > 0).
* `-g` 		negative power-law exponent gamma P(k)~k^(-gamma) (float >= 2).
* `-a` 		radial component a (float >= 1).
* `-nu`		scaling parameter nu (float > 0).
* `-radius`	rescaled radius of the hyperbolic ball (float > 0).
* `-v` 		optional switch whether or not to export node coordinates after generation (`<filename>.coord.dat`).
* `-seed` 	option to provide pseudorandom generator seed (long != 0).

#### Required Parameters
* Required: `-f AND -n AND -d`.
* `-u` user-based mode: `-k AND -c AND [-g OR -a]`.
* `-h` hybrid mode: `-k AND -t AND [-g OR -a]`.
* `-m` model-based mode: `-t AND [-nu OR -radius] AND [-g OR -a]`.
* Optional: `-v`, `-seed`.

### Output
The following files will be created in the project directory.

* `<filename>.dat` 		the network in edge list format separated by a whitespace.
* `<filename>.meta.dat` 	a meta file with the chosen/computed network parameters.
* `<filename>.coord.dat`	(optional) a list of the node coordinates in the hyperbolic ball.

### Examples
Hybrid mode: generate a network with n = 1000 nodes for dimensionality d = 3, average degree <k> = 10, power-law exponent gamma = 2.1 in the cold regime (tau = 1/2).
```
./generate_rhg -f example.dat -n 1000 -d 3 -h -g 2.1 -k 10 -t 0.5
```

Model-based mode: generate a network with n = 1000 nodes for dimensionality d = 1, radial component a = 1.0, scaling parameter nu = 0.5 in the hot regime (tau = 3/2) and export the coordinates afterwards with the `-v` switch.
```
./generate_rhg -f example.dat -n 1000 -d 1 -m -a 1.0 -nu 0.5 -t 1.5 -v
```