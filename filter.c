/* PLUGIN file to use with plantri.c
	Max Croucher, 2025
	This plugin is distributed at https://github.com/Max-Croucher/graph-generation
	Implements a variety of filtering methods to filter planar graphs for longest
	path length or Hamiltonicity. By default, this plugin will output graphs that
	have longest path length shorter than 2n/3. This length can be overwritten
	with -L To filter for non-Hamilonian graphs instead, the -H flag must be used.

	This plugin is tested only with graphs which are (mostly) polyhedral! i.e.
	polyhedral graphs, simplicial polyhedra, and Apollonian networks that aren't K3.
	
	Switches:
    * 			-L <int> determines the longest path length that graphs will be
					filtered for. No graph with a path this long will be output.
					Defaults to 2n/3. For example, to filter for nontraceable graphs
					(i.2. skip all graphs with longest paths with length n-1),
					set L to be n-1.
	*			-r <int> determines the maximum recursive depth of the snake heuristic.
	*			-W will enable the Floyd-Warshall Algorithm, which will compute
					the diameter of a graph to infer the existence of a path of
					length 3d-2. Assumes graphs are 3-connected.
    * 			-y will make Plantri output a table of longest path lengths. Is slow.
	* 			-S will set the splithint, if modulo splitting below the default
					is desired.
	*			-z will invert the filter to output graphs caught by the filter,
					instead of graphs that pass. Incompatible with pre-filtering
					techniques, such as -k.
	*			-D <int> will cause the plugin to only output graphs with the given
					diameter.
	*			-H will cause the plugin to output non-Hamiltonian graphs instead
					of graphs with sufficiently short longest paths.
    * 			-k <int> determines the maximum recursive depth of pre-filtering.
					This probably only works for polyhedral graphs, as it implements
					the PRE_FILTER_POLY filter.
	* 			-w will prevent the discarding of graphs at the recursive depth
					set by -k. Used for complex pre-filtering.

	Switch Compatibility:
	*	If -H and -W are not set, the algorithm will operate in longest-path mode,
		for which -k and -w are unusable
	*	If -r and -W are both set in longest-path mode, the switch which is set
		first will have its corresponding algorithm used first.
	*	If -D is set, the snake heuristic, -L and -r, -k, -w are unusable. -W is
		ignored.
	*	If -H is set, the algorithm will operate in Hamiltonian mode, for which
		-D, -L, -W are unusable.
	*	-z and -k are incompatible.


	Copyright for this plugin is held by the author
	Max Croucher, University of Canterbury, mpccroucher@gmail.com

	Copyright for Plantri is held jointly by its authors
	Gunnar Brinkmann, University of Gent, gunnar.brinkmann@ugent.be
	Brendan McKay, Australian National University, brendan.mcKay@anu.edu.au

	Licensed under the Apache License, Version 2.0 (the "License");
	you may not use this software except in compliance with the License.
	A copy of the License is included in the package and you can also
	view it at

		https://www.apache.org/licenses/LICENSE-2.0

	Unless required by applicable law or agreed to in writing, software
	distributed under the License is distributed on an "AS IS" BASIS,
	WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
	See the License for the specific language governing permissions and
	limitations under the License.
*/

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#define FILTER plugin_filter
#define PRE_FILTER_POLY prefilter()
#define SUMMARY summary
#define PLUGIN_INIT plugin_init()
#define SWITCHES "[-uagsETh -Ac#txm#P#bpe#f#qQ -odGVX -L#r#WySzD#Hk#w -v]"

static int current_graph_ham = 2;
static int ham_stack[256];
static int current_path_length=0;

static uint64_t total_found=0;
static uint64_t snake_found=0;
static uint64_t floyd_found=0;
static uint64_t graphs_seen=0;
static uint64_t graphs_checked=0;
static uint64_t path_length_table[256];
static uint64_t search_for_path=0;
static int64_t max_rec_depth = -1;
static float max_path_len;

static int zvalue = 0; //switch L
static int floyd_flag = 0; // switch W
static int do_diam = 0; // switch D
static int do_ham = 0; // switch H
static int do_truncate = 1; // switch w
static int invert_out = 0; // switch z
static int snake_depth = 0; // switch r
static int make_table = 0; // switch y
static int splitflag = -1; // switch Y

#define PLUGIN_SWITCHES \
	else if (arg[j] == 'L') zvalue = getswitchvalue(arg,&j); \
	else if (arg[j] == 'W') floyd_flag = 1; \
	else if (arg[j] == 'D') {do_diam = getswitchvalue(arg,&j); floyd_flag = 1;} \
	else if (arg[j] == 'H') do_ham = 1; \
	else if (arg[j] == 'z') invert_out = 1; \
	else if (arg[j] == 'r') snake_depth = getswitchvalue(arg,&j); \
	else if (arg[j] == 'k') max_rec_depth = getswitchvalue(arg,&j); \
	else if (arg[j] == 'w') do_truncate = 0; \
	else if (arg[j] == 'y') make_table = 1;  \
	else if (arg[j] == 'S') splitflag = getswitchvalue(arg,&j);


static void plugin_init()
/*  Check flag compatibility and print introductory messages
*/
{
INCOMPAT(invert_out&&(max_rec_depth!=-1),"-z","-k");
INCOMPAT(do_diam&&do_ham, "-D", "-H");
PERROR(!do_ham&&(max_rec_depth!=-1), "-k requires -H");
PERROR(!do_ham&&(!do_truncate), "-w requires -H");
INCOMPAT(do_diam&&(zvalue!=0), "-D", "-L");
INCOMPAT(do_diam&&snake_depth, "-D", "-r");
INCOMPAT(do_diam&&(max_rec_depth!=-1), "-D", "-k");
INCOMPAT(do_diam&&(!do_truncate), "-D", "-w");
INCOMPAT(do_ham&&(zvalue!=0), "-H", "-L");
INCOMPAT(do_ham&&(floyd_flag), "-H", "-W");
if (do_ham) {
	printf("Filtering for non-Hamiltonian graphs.\n");
} else if (do_diam) {
	printf("Filtering for graphs with diameter %d.\n", do_diam);
} else {
	if (floyd_flag) printf("Using the Floyd-Warshall Algorithm.\n");
	if (zvalue) {
		max_path_len = zvalue;
		printf("Filtering for graphs with longest path length < %.3f (overwritten by flag)\n",
			   max_path_len);
	} else {
		max_path_len = ((float)(2 * maxnv) / 3);
		printf("Filtering for graphs with longest path length < %.3f.\n",
			   max_path_len);
	}
}
if (snake_depth) printf("Using the Snake Heuristic with maximum depth %d.\n",
	     				snake_depth);
if (max_rec_depth != -1) {
	if (do_truncate) {
		printf("Truncating generation at depth %ld.\n", max_rec_depth);
	} else {
		printf("Using complex filtering at depth %ld.\n", max_rec_depth);
	}
};
if (splitflag != -1) {splithint = splitflag; printf("Set splithint to %d\n.", splithint);}
}


static void build_path(uint64_t *path_vertices, int *current_vertex, int *parents,
					   int backwards)
/*  Used by the Snake Heuristic to abitrarily construct a maximal path. Paths are
	stored by tracking vertex membership bitwise with path_vertices, and the order
	of vertices in a path is stored in the parent array. current_vertex is the
	endpoint of the array, and backwards indicates the direction to build the path.
*/
{
	int keep_adding = 1;
	while (keep_adding) {
		int found_next = 0;
		int next_vertex = (backwards) ? -1 : *current_vertex;
		EDGE *e, *elast;
		e = elast = firstedge[*current_vertex];
		do {
			if (!(*path_vertices & (1ULL << e->end))) {
				found_next = 1;
				next_vertex = e->end;
				if (backwards) {
					parents[*current_vertex] = next_vertex; parents[next_vertex] = -1;
				} else {
					parents[next_vertex] = *current_vertex;
				}
				*path_vertices |= 1ULL << next_vertex;
				*current_vertex = next_vertex;
			}
			e = e->next;
		} while ((!found_next) && (e != elast));
		if (!found_next) keep_adding = 0;
	}
}


static int pivot(int* parents, int last_vertex, int pivot)
/*  Implements the Pivot() operation. A path (encoded in the parent array) with
	endpoint last_vertex is pivoted to a vertex adjacent to last_vertex. This
	operation reorders the vertices in the path by reversing the order of vertices
	after pivot. The resulting vertex sequence is a valid path since last_vertex
	and pivot are adjacent.
*/
{
	int intermediary_vertices[maxnv];
	int num_swaps=0;
	while (last_vertex != pivot) {
		intermediary_vertices[num_swaps] = last_vertex;
		last_vertex = parents[last_vertex];
		num_swaps++;
	}
	if (num_swaps == 1) return last_vertex;
	int curr_parent = pivot;
	for (int i=0; i<num_swaps; i++) {
		parents[intermediary_vertices[i]] = curr_parent;
		curr_parent = intermediary_vertices[i];
	}
	return curr_parent;
}


static int hamitlonian_snake(uint64_t *path_vertices, int first_vertex,
							 int last_vertex, int* parents, int max_depth,
							 int previous_pivot)
/*  Main recursive function of the Snake Heuristic. This function begins with an
	end maximal path, and pivots on every possible internal vertex, calling
	hamitlonian_snake() on each path encountered. A maximum depth is specified,
	and to prevent wasted computation from pivoting onto the previously pivoted
	vertex, a previous pivot is specified.
*/
{
	int length = __builtin_popcount(*path_vertices); //count bits in integer
	if (length > current_path_length) current_path_length = length;
	EDGE *e, *elast;
	if ((*path_vertices + 1) >> maxnv) { //If path is on n vertices: check if ham
		e = elast = firstedge[last_vertex];
		do {
			if (first_vertex == e->end) return 1; // check if hamiltonian
			e = e->next;
		} while (e != elast);
	}
	if (!max_depth) return 0;
	
	e = elast = firstedge[last_vertex];
	do {
		uint64_t curr_path = *path_vertices;
		int v_piv = e->end;
		if ((v_piv != parents[last_vertex]) && (v_piv != previous_pivot)) {
			int pivot_parents[maxnv];
			memcpy(pivot_parents, parents, sizeof(int)*maxnv);
			int new_last = pivot(pivot_parents, last_vertex, v_piv);
			build_path(&curr_path, &new_last, pivot_parents, 0); //attempt to grow
			int result = hamitlonian_snake(&curr_path, first_vertex, new_last,
										   pivot_parents, max_depth-1, v_piv);
			if (result) { //Return with success if Hamiltonian
				*path_vertices = curr_path;
				return 1;
			}
		}
		e = e->next;
	} while (e != elast);
	return 0;
}


static int hamiltonian_loop(uint64_t path_vertices, int current_vertex)
/*  Brute-Force Hamiltonian classifier. This function recursively builds paths
until a Hamiltonian cycle is found. If the current graph is non-Hamiltonian, all
paths with an abitrary vertex as an endpoint must be checked. Is slow on large graphs.
*/
{
	int length = __builtin_popcount(path_vertices);
	if (length > current_path_length) current_path_length = length;
	
	EDGE *e, *elast;
	e = elast = firstedge[current_vertex];
	do {
		int vertex = e->end;
		if (((path_vertices + 1) >> maxnv) && vertex == 0) return 1;
		if (!(path_vertices & (1ULL << vertex))) { // if edge 'vertex' is not in path
			int recursive_result = hamiltonian_loop(path_vertices | 1ULL << vertex, vertex);
			if (recursive_result) return 1;
		}
		e = e->next;
	} while (e != elast);
	return 0;
}


static int is_hamiltonian()
/*  This function determines if the current graph is Hamiltonian by running the
	Snake Heuristic, and resorting to is-hamiltonian if the heuristic fails. This
	function is the fuction that should be called to determine Hamiltonicity.
*/
{
	graphs_checked++;
	if (snake_depth == 0) { //Skip snake if snake is not used.
		uint64_t recursive_path = 1;
		int recursive_start_vertex = 0;
		if (hamiltonian_loop(recursive_path, recursive_start_vertex)) {
			total_found++;
			return 1;
		}
		return 0;
	}
	
	uint64_t path_vertices = 1; //Set up empty path
	int current_vertex = 0;
	int parent[maxnv];
	parent[current_vertex] = -1;
	build_path(&path_vertices, &current_vertex, parent, 0); //0 means forwards
	int last_vertex = current_vertex;
	int first_vertex = 0;
	build_path(&path_vertices, &first_vertex, parent, 1);
	int is_snake_ham = hamitlonian_snake(&path_vertices, first_vertex, last_vertex,
										 parent, snake_depth, first_vertex);
	if (is_snake_ham) { // If the Snake Heuristic worked
		snake_found ++;
		total_found++;
		return 1;
	} else { // Otherwise call the brute-force
		uint64_t recursive_path = 1;
		int recursive_start_vertex = 0;
		int result = hamiltonian_loop(recursive_path, recursive_start_vertex);
		if (result) total_found++;
		return result;
	}
}


static int longest_path(uint64_t *path, int last_vertex, int curr_len)
/* Similar to hamiltonian_loop(). Returns the length of the longest path in the current
	graph.
*/
{
	if (maxnv == curr_len) {
		return curr_len;
	}
	uint64_t current_longest = 0ULL;
	int current_length = curr_len;
	EDGE *e, *elast;
	e = elast = firstedge[last_vertex];
	
	do {
		int vertex = e->end;
		if (!(*path & (1ULL << vertex))) {
			uint64_t next_path = *path + (1ULL << vertex);
			int next_path_len = longest_path(&next_path, vertex, curr_len + 1);
			if (next_path_len == maxnv) {
				*path = next_path;
				return next_path_len;
			}
			if (next_path_len > current_length) {
				current_longest = next_path;
				current_length = next_path_len;
			}
		}
		e = e->next;
	} while (e != elast);
	*path = current_longest;
	return current_length;
}


static int get_longest_path(void)
/*  Initialises the recursive longest_path() to determine the length of the longest
	path in the current graph
*/
{
	search_for_path++;
	int path_length = 0;
	for (int i=0; i<(maxnv); i++) {
		uint64_t p = 1ULL << i;
		path_length = longest_path(&p, i, 1);
		if (path_length == maxnv) return path_length;
	}
	return path_length;
}


static int floyd_warshall()
/*  Implements the Floyd Warshall Algorithm, which constructs the distance array
of a graph and computes the diameter of this graph
*/
{
	int8_t distance_matrix[maxnv][maxnv];
	for (int v=0; v<maxnv; v++) {
		for (int u=0; u<maxnv; u++) {
			distance_matrix[v][u] = (u == v) ? 0 : -1;
		}
		EDGE *e, *elast;
		e = elast = firstedge[v];
		do {
			distance_matrix[v][e->end] = 1;
			e = e->next;
		} while (e != elast);
	}
	
	for (int k=0; k<maxnv; k++) {
		for (int i=0; i<maxnv; i++) {
			for (int j=0; j<maxnv; j++) {
				if ((distance_matrix[i][k] >= 0 && distance_matrix[k][j] >= 0)
					&& ((distance_matrix[i][j] < 0) || (distance_matrix[i][j] >
					distance_matrix[i][k] + distance_matrix[k][j])))
				{
					distance_matrix[i][j] = distance_matrix[i][k] + distance_matrix[k][j];
				}
			}
		}
	}
	
	int diameter = 0;
	for (int v=0; v<maxnv; v++) {
		for (int u=0; u<maxnv; u++) {
			if (diameter < distance_matrix[u][v]) diameter = distance_matrix[u][v];
		}
	}
	return diameter;
}


static int snake_recurse_path(uint64_t *path_vertices, int first_vertex,
							  int last_vertex, int* parents, int max_depth,
							  int previous_pivot)
/*  Implements the Snake Heuristic for finding a long path in a graph. Similar
to hamitlonian_snake()
*/
{
	EDGE *e, *elast;
	
	int path_length = 0;
	for (int i = 0; i<nv; i++) {
		if (*path_vertices & 1ULL<<i) path_length++;
	}
	if (path_length-1 >= max_path_len) return 0;
	if (!max_depth) return 0; // halt if at max depth
	
	e = elast = firstedge[last_vertex];
	do {
		uint64_t curr_path = *path_vertices;
		int v_piv = e->end;
		if ((v_piv != parents[last_vertex]) && (v_piv != previous_pivot)) {
			int pivot_parents[nv];
			memcpy(pivot_parents, parents, sizeof(int)*nv);
			int new_last = pivot(pivot_parents, last_vertex, v_piv);
			build_path(&curr_path, &new_last, pivot_parents, 0);
			int result = snake_recurse_path(&curr_path, first_vertex, new_last,
											pivot_parents, max_depth-1, v_piv);
			if (result) {
				*path_vertices = curr_path;
				return 1;
			}
		}
		e = e->next;
	} while (e != elast);
	return 0;
}


static int long_path_loop(uint64_t path_vertices, int current_vertex,
						 int current_length)
/* Recursively Determines if the current graph has a sufficiently long path.*/
{
	if (current_length-1 >= max_path_len) {
		//check if current path is too long. current_length counts the
		//NUMBER_OF_VERTICES, and thus current_length-1 countes the number of edges
		return 1;
	}
	EDGE *e, *elast;
	e = elast = firstedge[current_vertex];
	
	do {
		int vertex = e->end;
		if (!(path_vertices & (1ULL << vertex))) { // if edge 'vertex' not in path
			int recursive_result = long_path_loop(path_vertices | (1ULL << vertex),
												 vertex, current_length+1);
			if (recursive_result) return 1;
		}
		e = e->next;
	} while (e != elast);
	return 0;
}


static int has_long_path()
/*  Main procedure for determining if a graph has a sufficiently long path. Calls
	snake_recurse_path() if the snake heuristic is to be used, and then the exhaustive
	long_path_loop if it fails.*/
{
	graphs_checked++;
	if (snake_depth == 0) {
		for (int i = 0; i<nv; i++) {
			uint64_t recursive_path = 1ULL << i;
			if (long_path_loop(recursive_path, i, 1)) {
				total_found++;
				return 1;
			}
		}
		return 0;
	}
	
	uint64_t path_vertices = 1;
	int current_vertex = 0;
	int parent[nv];
	parent[current_vertex] = -1;
	build_path(&path_vertices, &current_vertex, parent, 0); //0 means forwards
	int last_vertex = current_vertex;
	int first_vertex = 0;
	build_path(&path_vertices, &first_vertex, parent, 1);
	int finally_ham = snake_recurse_path(&path_vertices, first_vertex, last_vertex,
										 parent, snake_depth, first_vertex);
	if (finally_ham) {
		snake_found++;
		total_found++;
		return 1;
	} else {
		for (int i = 0; i<((float)nv-1); i++) {
			uint64_t recursive_path = 1ULL << i;
			if (long_path_loop(recursive_path, i, 1)) {
				total_found++;
				return 1;
			}
		}
		return 0;
	}
}


static int plugin_filter(int nbtot, int nbop, int doflip)
/*  Main entrypoint for the filter, and uses switches (-H, -D) to determine which
	filtering type is used. If -H (non-Hamiltonian) is used, and pre-filtering is
	enabled, this function starts by inheriting information about the graph that
	genenerated the current graph, and can skip classifying the current graph if
	the predecessor is non-Hamiltonian, or if the complex-filtering depth is exceeded.
	For other filtering types, relevent filters are called, and finally this function
	updates some stats
*/
{
	graphs_seen++;
	int result;
	if (do_ham) { // Hamiltonian mode
		if (max_rec_depth != -1) {
			int degeneracy = (6*nv-12-ne)/2;
			if (degeneracy == 0) {
				current_graph_ham = is_hamiltonian();
				if (do_truncate && max_rec_depth==1 && current_graph_ham==1) {
					ham_stack[degeneracy]=current_graph_ham;
					return !(1 ^ invert_out);
				}
			} else if (do_truncate && degeneracy >= max_rec_depth
					   && ham_stack[degeneracy-1]==1)
			{
				ham_stack[degeneracy]=1;
				return !(1 ^ invert_out);
			} else if (do_truncate && degeneracy > max_rec_depth) {
				current_graph_ham = ham_stack[degeneracy];
			} else if (degeneracy>0 && ham_stack[degeneracy-1]==0) {
					current_graph_ham = 0;
			} else {
				current_graph_ham = is_hamiltonian();
			}
			ham_stack[degeneracy]=current_graph_ham;
		} else {
			current_graph_ham = is_hamiltonian();
		}
		result = !(current_graph_ham ^ invert_out);
	} else if (do_diam) { // Diameter mode
		graphs_checked++;
		int diameter = floyd_warshall();
		if (diameter == do_diam) {
			total_found++;
			floyd_found++;
		}
		result = diameter == do_diam;
	} else { // Longest path mode
		int long_path_loop = 0;
		if (floyd_flag) {
			int diameter = floyd_warshall();

			if (3*diameter - 2 >= max_path_len) {
				//Path formed by exploiting 3-connectivity
				floyd_found++;
				long_path_loop = 1;
			}
		}
		if (!long_path_loop) long_path_loop = has_long_path();
		result = (!long_path_loop) ^ invert_out;
	}

	current_path_length = 0;
	if (make_table && result) path_length_table[get_longest_path()]++;
	if (graphs_seen % 10000000 == 0) {
		printf("Seen %ld graphs. Checked %ld graphs. Found %ld graphs.\n", graphs_seen, graphs_checked, 
			   invert_out? total_found : graphs_checked - total_found);
		fflush(stdout);
	}
    return result;
}


static int prefilter()
/*  This function is invoved by PRE_FILTER_POLY, and is only used if pre-filtering
	in Hamiltonian mode (-k and -H). This function acts as flow control for
	scanpoly_c3(). If simple pre-filtering is being used and the recursive depth
	is met, this function halts generation from the current graph.
*/
{
	if (do_ham && max_rec_depth != -1 && do_truncate) {
		if ((6*nv-12-ne)/2 >= max_rec_depth) {
			return current_graph_ham==0;
		}
	}
	return 1;
}


static void summary()
/*  This function runs after graph generation and prints some stats.
*/
{
	printf("\n");
	printf("Seen %ld graphs.\n", graphs_seen);
	printf("Filters applied to %ld graphs.\n", graphs_checked);
	if (snake_depth > 0) printf("Graphs found using Snake: %ld.\n", snake_found);
	if (floyd_flag) printf("Graphs found using Floyd: %ld.\n", floyd_found);
	if (make_table && (oswitch ? totalout_op : totalout)) {
		printf("\nPath Length Table:\nlen : Number of graphs with this length\n");
		for (int i=0; i < 256; i++) {
			if (path_length_table[i] > 0) printf("%3d : %ld\n", i-1,
				path_length_table[i]);
		}
		printf("\n");
	}
	if (splitflag != -1) printf("Splitting level was %d.\n", splitlevel);
}
