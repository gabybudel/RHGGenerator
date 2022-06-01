# RHG Generator #

Network generator for random hyperbolic graphs (RHGs) with a latent representation in a hyperbolic ball of any dimension $d+1 \geq 2$. The generator returns the network in edgelist format and additionally a list of the node coordinates if desired. The generator can be operated in 2 different modes, see *Required Parameters* for lists of the required parameters per mode.

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
Call the executable `generate_rhg` with the mode of choice (hybrid or model-based) and provide the required parameters. Parameter flags and values should be separated by a whitespace, e.g., `-n 1000` sets the number of nodes to 1,000. A full list of the available parameters and which parameters are required in which mode can found below. The inputs are tested for the restriction within parentheses.

#### Parameters
* `-f` 		filename (either a `*.dat` filename or a filename without extension).
* `-n` 		network size (`int > 1`).
* `-d` 		dimensionality parameter *d* for the dimensionality *d + 1* of the hyperbolic ball (`int >= 1`).
* `-h`		hybrid-mode (select one from `{-h, -m}`.
* `-m` 		model-based mode (select one from `{-h, -m}`.
* `-t` 		rescaled temperature $\tau$ (`float >= 0`).
* `-g` 		negative power-law exponent $\gamma$ for $P(k) \sim k^{-\gamma}$ (`float >= 2`).
* `-a` 		radial component $\mathcal{a}$ (`float >= 1`).
* `-nu`		scaling parameter $\nu$ (`float > 0`).
* `-radius`	rescaled radius $\mathcal{R}$ of the hyperbolic ball (`float > 0`).
* `-v` 		(optional) a switch whether or not to export node coordinates after generation (`<filename>.coord.dat`).
* `-seed` 	(optional) the seed for the pseudorandom generator (`long != 0`).

#### Required Parameters
* Always required: `-f AND -n AND -d`.
* `-h` hybrid mode: `-k AND -t AND [-g OR -a]`.
* `-m` model-based mode: `-t AND [-nu OR -radius] AND [-g OR -a]`.

### Output
The following files will be created in the project directory.

* `<filename>.dat` 		the network in edge list format with node indices separated by a whitespace.
* `<filename>.meta.dat` 	a meta file with the chosen/computed network parameters.
* `<filename>.coord.dat`	(optional) a list of the node coordinates in the hyperbolic ball.

### Examples
Hybrid mode: generate a network with $n = 10^{4}$ nodes for dimensionality $d = 3$, average degree $\langle k \rangle = 10$, power-law exponent $\gamma = 2.1$ in the cold regime ($\tau = 1/2$).
```
./generate_rhg -f example.dat -n 10000 -d 3 -h -g 2.1 -k 10 -t 0.5
```

Model-based mode: generate a network with $n = 10^{3}$ nodes for dimensionality $d = 1$, radial component $\mathcal{a} = 1.0$, scaling parameter $\nu = 0.5$ in the hot regime ($\tau = 3/2$) and export the coordinates with the `-v` switch afterwards.
```
./generate_rhg -f example.dat -n 1000 -d 1 -m -a 1.0 -nu 0.5 -t 1.5 -v
```

### Numerical considerations
To prevent numerical issues within the computations, the following cut-off values in the parameters are applied based on the user inputs:
* When $\tau < 0.05$, links are generated deterministically as if $\tau = 0$ (and when $\tau < 0$, an error is thrown).
* When $|\tau - 1| < 0.01$, the scaling of the critical regime $\tau = 1$ is invoked.
* When $\tau \geq 1.01$, the scaling of the hot regime $\tau > 1$ is invoked, and the parameters must adhere to $\gamma \leq \tau + 1$, otherwise an error is thrown (because the parameters imply $\mathcal{a} < 1$).
* When $\mathcal{a} < 1.01$, the scaling of the case $\mathcal{a}=1$ is invoked (and when $\mathcal{a} < 1$, an error is thrown).
* When $\frac{\alpha}{2} \mathcal{R} > 600$, all radial coordinates will be very large and they cannot be reliably generated, therefore we resort to degenerate hyperbolic graphs where each radial coordinate $\mathcal{r}_{i} = \mathcal{R}$.
* In the critical regime $\tau = 1$, $\mathcal{a} > 1$, the Lambert function $W_{-1}(x)$ is evaluated using a look-up table of function values and linear interpolation, closely resembling the true function values. 