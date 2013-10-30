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
except for cluster.py, which uses a custom Cython library that has to
be compiled. For that, just cd to /cluster and type 'make'. A new file,
named rejection.so, will be created, and then cluster.py will be ready
for execution.


## Usage

### cluster.py

cluster.py [-h] [--gas-core] [--dm-core] [--dm-only] [-o init.dat]

Generate an initial conditions file for a galaxy cluster halo simulation.

Optional arguments:
  -h, --help   show this help message and exit
  --gas-core   Sets the density profile for the gas to have a core.
  --dm-core    The same, but for the dark matter.
  --dm-only    Generates an initial conditions file containing only dark
               matter.
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


## Author:
    Rafael Ruggiero
    Undergraduate student at Universidade de São Paulo (USP), Brazil
    Contact: bluewhale [at] cecm.usp.br

Credits for Dr. Rubens Machado (http://www.astro.iag.usp.br/~rgmachado/),
for the vital support and suggestions.

