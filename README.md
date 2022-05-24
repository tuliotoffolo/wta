# Branch-and-adjust for a nonlinear Assignment Problem

Written by **Túlio Toffolo** and **Konstantin Pavlikov**.

Please address all contributions, suggestions, and inquiries to the current project administrator.

### More information

This repository contains a solver for the Weapon-Target Assignment Problem (WTA). 

The developed approach is detailed in the following [article](https://link.springer.com/article/10.1007/s10479-022-04525-6):

```bib

@article{AndersenPavlikovToffolo2020,
    author = {Andersen, Alexandre Colaers and Pavlikov, Konstantin and Toffolo, T{\'u}lio A. M.},
    date = {2022/01/13},
    doi = {10.1007/s10479-022-04525-6},
    isbn = {1572-9338},
    journal = {Annals of Operations Research},
    title = {Weapon-target assignment problem: exact and approximate solution algorithms},
    url = {https://doi.org/10.1007/s10479-022-04525-6},
    year = {2022},
    bdsk-url-1 = {https://doi.org/10.1007/s10479-022-04525-6}}
}
```

To cite the code or the paper, please use the bibtex above.

### Usage:

```text
Usage: python3 wta_cplex.py <input> [options]
    <input> : Path of the problem input file.

Options:
    -approach <approach>   : selected approach/formulation to execute; possible values are:
                             {branch-and-adjust, probchain, underapprox, upperapprox}
                             (default: branch-and-adjust).
    -branching <branching> : branching strategy from {cplex, probabilities} (default: probabilities).
    -delta <delta>         : delta value (default: 0.001).
    -emphasis <emphasis>   : cplex search emphasis (default: 3).
    -output <output>       : output folder to save the log and additional info (default: output).
    -solution <file>       : file to save final solution (default: null).
    -threads <threads>     : maximum number of threads (default: 1).
    -timelimit <timelimit> : runtime limit (default: inf).
```

Examples:
- ``python3 wta_cplex.py data/wta_50x100x1.txt``
- ``python3 wta_cplex.py data/wta_50x100x1.txt -delta 0.00001 -solution solutions/wta_50x100x1.sol -threads 1``

### Instance and solution files

Instance and solution files are available in the ``data`` folder of this repository.

### Requirements

The entire code was written in Python 3 using *Cplex 12.10*'s library. Therefore, *Cplex* must be installed in your system.

## Questions?

If you have any questions, please feel free to contact us.

