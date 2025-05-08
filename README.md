# graph-generation
A repository of scripts and plugins for the generation and filtering of planar graphs.

This repository is currently incomplete. Included tools:
* plugin.c

# plugin.c

This plugin is designed for use with Brinkmann and McKay's graph generation tool Plantri, which is available at  https://users.cecs.anu.edu.au/~bdm/plantri/. This plugin implements a variety of graph filtering methods, to filter graphs for structural properties including Hamiltonicity and longest path length. It is implemented for the purposes of filtering polyhedral graphs, and thus other classes of planar graphs may not be correctly filtered. To compile Plantri with this plugin, the supplied makefile can be used.