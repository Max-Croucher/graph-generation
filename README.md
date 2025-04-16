# graph-generation
A repository of scripts and plugins for the generation and filtering of planar graphs.

This repository is currently incomplete. As-is and uncommented plugins are available, but more refined copies and other useful tools will be made available soon.

Messages printed to stdout aren't always useful and command-line switches aren't always intuitive.

Build the binaries with the makefile.

# Hamiltonian Filtering
To filter graphs for Hamiltonicity, use the Plantri plugin hamiltonian-filter.c

This version of hamiltonian-filter.c, which is included in this repository. Eventually, I will upload a version that runs on an unmodified version.

# Long Path Filtering
To filter graphs for longest path length, use the Plantri plugin long_path_filter.c

This version uses an unmodified Plantri, which you will need to get yourself, at https://users.cecs.anu.edu.au/~bdm/plantri/. I think there's an off-by-one error somewhere when manually specifying the length of paths to filter by. All will be fixed after my break.