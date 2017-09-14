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

    usage: clustep.py [-h] [--gas-core] [--dm-core] [--no-dm] [--no-gas]
                      [-o init.dat]
    
    Generates an initial conditions file for a galaxy cluster halo simulation.
    
    optional arguments:
      -h, --help   show this help message and exit
      --gas-core   Sets gamma=0 in the Dehnen density profile assigned to the gas
                   component, causing it to feature a central core. By default
                   gamma=1, which is equivalent to a Hernquist density profile.
                   See 1993MNRAS.265..250D and 1990ApJ...356..359H.
      --dm-core    Exactly the same as above, but for the dark matter component.
      --no-dm      No dark matter particles in the initial conditions. The dark
                   matter potential is still used when calculating the gas
                   temperatures.
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


## Works which used this code

* [Ruggiero & Lima Neto (2017)](http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1703.08550)


## Author

    Rafael Ruggiero
    Ph.D student at Universidade de São Paulo (USP), Brazil
    Contact: rafael.ruggiero [at] usp.br

Credits for Prof. Dr. Rubens Machado (http://paginapessoal.utfpr.edu.br/rubensmachado),
for the vital support and suggestions.

## Disclaimer

Feel free to use this code in your work, but please link this page
in your paper.
