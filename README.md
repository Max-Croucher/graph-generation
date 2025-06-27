# graph-generation
A repository of scripts and plugins for the generation and filtering of planar graphs.

This repository is currently incomplete. Included tools:
* plugin.c
* pathfinder_partial.c
* stellations/bin/all_ham_cycles.c
* stellations/bin/all_longest_paths.c
* stellations/read_pl.py
* stellations/draw_pl.py
* stellations/interactive_stellations.py

The C files can be compiled with the included makefile.

# plugin.c

This plugin is designed for use with Brinkmann and McKay's graph generation tool Plantri, which is available at  https://users.cecs.anu.edu.au/~bdm/plantri/. This plugin implements a variety of graph filtering methods, to filter graphs for structural properties including Hamiltonicity and longest path length. It is implemented for the purposes of filtering polyhedral graphs, and thus other classes of planar graphs may not be correctly filtered. To compile Plantri with this plugin, the supplied makefile can be used.
Full documentation for this plugin is included in the source code.

# pathfinder_partial.c

This program is written to iterate over the graphs in a planar code file (typically .pl), as output by the planar graph generation tool Plantri. For each graph in the given file, the program will compute the longest paths in this graph that can be generated with a given vertex as an endpoint. This tool also supports a range of constraints, including the blacklisting of vertices, or requirement that paths contain particular vertices or paths. By default, the program will display information on an encountered path if it is longer than any path found so far. This allows useful information to be obtained even if the program is prematurely halted.

Usage: `./pathfinder_partial <pl-filename> <start-vertex>`

The following optional arguments are also available:
* `--help` prints a help message and exits.
* `--halt-on-trace`: Halts the graph generation procedure for a graph if a trace is found.
* `--print`: Prints a generated path if it is equally as long as the longest path found so far, instead of only printing generated paths if they're longer.
* `--used-edges`: Prints an adjacency matrix of all edges used by any longest path. Vertex pairs that include a blacklisted vertex are displayed as an asterisk \"*\". Vertex pairs that do not correspond to an edge are displayed with a period \".\". Vertex pairs that correspond to an edge are displayed as a \"1\" or \"0\" based on whether or not there exists a longest path that uses this edge.
* `--blacklist-vertices [sequence of vertices]`: Specifies a sequence of vertices that cannot be present in a longest path.
* `--require-vertices [sequence of vertices]`: Specifies a sequence of vertices that must be present in a longest path.
* `--require-edges [sequence of hyphen-delimetered sequences of vertices]`: Specifies a sequence of paths that must be present in a longest path. These required paths must be vertex disjoint and cannot involve blacklisted vertices.

Example:

`./pathfinder_partial example_graph.pl 0 --blacklist-vertices 2 14 --require-vertices 3 --require-paths 0-4-15-29-7 11-27-16`

This command would find the longest paths in each graph given in `example_graph.pl` that has vertex 0 as an endpoint, provided each path does not use vertices 2 or 14, does use vertex 3, and uses the sub-paths 0-4-15-29-7 and 11-27-16.

# stellations/bin/all_ham_cycles.c

This program is used to obtain every Hamiltonian cycle of every graph in a planar `.pl` file. This program is required for `stellations/interactive_stellations.py`, and can be compiled using the supplied makefile. Full documentation is included in the source code.

# stellations/bin/all_longest_paths.c

This program is used to obtain every longest path of every graph in a planar `.pl` file. This program is required for `stellations/interactive_stellations.py`, and can be compiled using the supplied makefile. Full documentation is included in the source code.

# stellations/read_pl.py

This script is written to aid in the reading of planar `.pl` files output by Plantri. The script can be run in the in the following manner: `python3 read_pl.py <filename> [mode]`, where `mode` is one of the following instructions:
* `count`: In this mode, the number of graphs in the planar file is printed stdout.
* `print`: In this mode, each graph in the planar file is printed to stdout in the form of an adjacency list.
* `store`: In this mode, the adjacency list for each graph in the planar file is written to `filename.txt`. This is equivalent to piping the output of the `print` mode to a file.
* `split <integer n>`: In this mode, a directory `filename.split` is created, and the graphs in the given planar file are evenly partitioned into `n` numbered files in this directory.
* If no mode or filename is specified, a list of the available modes is printed.

# stellations/draw_pl.py

This script allows the graphs in a planar `.pl` file output by Plantri to be drawn as an image, which helps visualising the structure of a graph. The script can be run with the command `python3 draw_pl <filename>`, and will iterate through the graphs in the planar file, drawing each graph individually to the screen. If run with the optional `--save <directory>` argument, then each graph will be individually written to the specified directory as a `.png` file.

# stellations/interactive_stellations.py

This script lets planar graphs be interactively manipulated and stellated to better view some of the structural properties of these graphs. Through the somewhat poor decision making of the author, this script uses matplotlib interactivity to function, and this makes the script unstable. To use this script, the aforementioned programs `all_ham_cycles` and `all_longest_paths` are required, and thus must be compiled. This can be done using the makefile.

In its simplest form, the script can be run with the command `python3 interactive_stellations.py <filename>`, for which the first graph in the file will be displayed. If this graph is large (at least 15 vertices), then there may be a noticeable delay 

The edges in the graph are coloured to represent the configuration of Hamiltonian cycles in the graph. A black edge is used in no Hamiltonian cycle, a red edge is used in at least one Hamiltonian cycle, and a blue edge is used in every Hamiltonian cycle. If the graph is non-Hamiltonian, then all edges will be black.


Next, the following features are available:
* The interactive tool will begin in stellate mode. In this mode, clicking on a face of the graph (or outside the graph for the outer face) will cause that face to be stellated, and the edge colourings recalculated. In addition to this, the vertices can be dragged to different positions.
* Clicking the Delete button switches the tool to delete mode. In this mode, clicking a vertex or edge will delete it, and cause the edge colourings to be recalculated.
* Clicking the Change Outer Face button will let a face be selected as the new outer face. This will cause the graph to be redrawn, but the edge colours will not be recalculated.
* Clicking the Redraw Graph button will cause the graph to be redrawn.
* Clicking the Toggle Gradient Colouring button will cause the edge colours change to reflect the proportion of edges that use them.
* Clicking the Toggle Frequencies button will add a label to each non-black edge to display the proportion of Hamiltonian cycles that use this edge. Clicking the button again causes the labels to count the number of such cycles instead.
* The Forward and Back buttons allow stellations and deletions to be undone or redone.

In addition to this, the tool can be used with different modes or features:
* Running the script with `--index <integer n>` or `-i <integer n>` will cause the `n+1`th graph in the given planar file to be displayed instead of the default 0.
* Running the script with `--longest` or `-l` will cause the tool to run in Longest Path Mode, in which the edge colourings reflect the configurations of longest paths in the graph instead of Hamiltonian cycles. Vertices that are the endpoint of at least one longest path will be coloured green.
* Running the script with `--nothing` or `-n` will cause the tool to run without any path/cycle generation. This is useful for generating planar files for graphs too large to check for structural properties interactively.
* Running the script with `--output <filename>` or `-o <filename>` will cause the script to save the current graph onscreen to the provided file whenever a graph is stellated or has an element deleted. The index of the currently-displayed graph is displayed so the program can be easily restarted with the desired graph using the `--index` argument.