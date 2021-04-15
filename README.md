**MacroPy --- SBI 2021**
==============================

MacroPy is a command line program which reconstructs a whole biological macro-complex using PDB's of its pairwise interactions as input, either protein-protein, protein-DNA or protein-RNA.




# Requirements
MacroPy needs [Python 3.7](https://www.python.org/downloads/) or higher and the package [Biopython](https://biopython.org/). [Modeller](https://salilab.org/modeller/download_installation.html) is necessary to run the -min option. Biopython can be installed with the following command:
```
pip install biopython
```

# Installation
You can download and install MacroPy by cloning it and then executing the setup.py script:
```
git clone https://github.com/nicodr97/MacroPy.git
python setup.py install
```

# Program usage
To run the program you have to execute the following command with the proper arguments:

```
usage: Macro.py [-h] -i INPUT_DIRECTORY -o OUTPUT_DIRECTORY [-c COMPLEX_NAME]
                [-f] [-s STOICHIOMETRY] [-v] [-min [MINIMIZATION]] [-pdb]
                [-mc MAX_CHAINS] [-it IDENTITY_THRESHOLD] [-Rt RMSD_THRESHOLD]
                [-ns NEIGHBOR_SEARCH_DISTANCE] [-cd CLASHES_DISTANCE]
                [-nc NUMBER_CLASHES]

MacroPy 1.0 - Reconstruct a whole biological macro-complex using PDBs of its
pairwise interactions as input, either protein-protein, protein-DNA or
protein-RNA.

required arguments:
  -i INPUT_DIRECTORY, --input-directory INPUT_DIRECTORY
                        Directory containing the input structure files
                        (default: None)
  -o OUTPUT_DIRECTORY, --output-directory OUTPUT_DIRECTORY
                        Create the output directories (default: None)

optional arguments:
  -c COMPLEX_NAME, --complex-name COMPLEX_NAME
                        Reconstructed complex name (default: Complex)
  -f, --force           Force overwriting if the output directory already
                        exists (default: False)
  -s STOICHIOMETRY, --stoichiometry STOICHIOMETRY
                        File containing the stoichiometry (default: None)
  -v, --verbose         Program log will be printed to standard error while
                        running (default: False)
  -min [MINIMIZATION], --minimization [MINIMIZATION]
                        Perform an energy minimization by Conjugate Gradients
                        algorithm with the specified number of steps or until
                        convergence (default: False)
  -pdb, --save-pdb      Besides the .cif file, save a .pdb file with up to 25
                        chains (default: False)
  -mc MAX_CHAINS, --max-chains MAX_CHAINS
                        Number of chains of the complex at which to stop
                        adding new chains (default: 180)
  -it IDENTITY_THRESHOLD, --identity-threshold IDENTITY_THRESHOLD
                        Minimum percentage of sequence similarity (between 0
                        and 1) to consider two PDB chains the same (default:
                        0.95)
  -Rt RMSD_THRESHOLD, --RMSD-threshold RMSD_THRESHOLD
                        Maximum RMSD value to consider two (similar) PDB
                        chains the same (default: 2.5)
  -ns NEIGHBOR_SEARCH_DISTANCE, --Neighbor-Search-distance NEIGHBOR_SEARCH_DISTANCE
                        Minimum distance between two PDB chains to consider
                        that they are actually interacting (default: 5)
  -cd CLASHES_DISTANCE, --clashes-distance CLASHES_DISTANCE
                        Maximum distance between atoms of two chains to
                        consider that they have clashes between them (default:
                        1.8)
  -nc NUMBER_CLASHES, --number-clashes NUMBER_CLASHES
                        Maximum number of close atoms to consider that two
                        chains are clashing (default: 20)

```


For more details about MacroPy check the [report](https://github.com/nicodr97/MacroPy/blob/master/doc/Report.pdf).

# Authors
Nicolás Díaz Roussel  
Alejandra Omaira Morcillo Nieto  
Francho Nerín Fonz  
