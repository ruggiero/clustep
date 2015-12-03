## Required libraries
(and the names of the packages in Ubuntu)
 
* NumPy (python-numpy)
* SciPy (python-scipy)
* Cython (cython), higher than 0.17.4
* Matplotlib (python-matplotlib)
* Argparse (python-argparse), in case you use Python < 3.2
* python-dev


## Installation

All the scripts can be executed without installation or compilation,
except for clustep.py, which uses a custom Cython library that has to
be compiled. For that, just cd to /cluster and type 'make'. A new file,
named optimized_funcions.so, will be created, and then clustep.py will
be ready for execution.


## Usage

### clustep.py

    clustep.py [-h] [--gas-core] [--dm-core] [--dm-only] [-o init.dat]

    Generate an initial conditions file for a galaxy cluster halo simulation.

    Optional arguments:
      -h, --help   show this help message and exit
      --gas-core   Sets the density profile for the gas to have a core.
      --dm-core    The same, but for the dark matter.
      --no-dm      No dark matter particles in the initial conditions. The dark
                   matter potential is still used when calculating the gas
                   temperatures
      --no-gas     No gas, only dark matter.
      -o init.dat  The name of the output file.

### profiles.py

    profiles.py [-h] [--gas-core] [--dm-core] file.dat

    Plots stuff.

    positional arguments:
      file.dat    The name of the input file.

    optional arguments:
      -h, --help  show this help message and exit
      --gas-core  Sets the density profile for the gas to have a core.
      --dm-core   The same, but for the dark matter.

Some analysis scripts are also included, you can try these out. I haven't
documented them because they are changed all the time and aren't all that
well written as of now.


## Author

    Rafael Ruggiero
    Ph.D student at Universidade de São Paulo (USP), Brazil
    Contact: rafael.ruggiero [at] usp.br

Credits for Dr. Rubens Machado (http://www.astro.iag.usp.br/~rgmachado/),
for the vital support and suggestions.

## Disclaimer

Feel free to use this code in your work, but please link this page
in your paper.
