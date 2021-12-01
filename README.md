# Branch-and-adjust for an Assignment Problem

Written by **Túlio Toffolo** and **Konstantin Pavlikov**.

Please address all contributions, suggestions, and inquiries to the current project administrator.

### More information

The developed approach is detailed in an article currently under revision. More information will be available here soon.

<!--
To cite the code or the paper, please use the following bibtex:

```bib
@article{AndersenPavlikovToffolo2020,
  author = "Alexandre Colaers Andersen and Konstantin Pavlikov and T\'{u}lio Angelo Machado Toffolo",
  title = "Weapon-Target Assignment Problem: Exact and Approximate Solution Algorithms",
  journal = "?",
  volume = "?",
  number = "?",
  pages = "?--?",
  year = "2020",
  DOI ="?"
}
```
-->

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

