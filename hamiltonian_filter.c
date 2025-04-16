/* PLUGIN file to use with plantri.c
	Max Croucher, 2025
	This plugin is distributed at https://github.com/Max-Croucher/graph-generation
	Implements a variety of filtering methods to filter planar graphs for Hamiltonicity.

	This plugin is tested only with graphs which are (mostly) polyhedral! i.e. polyhedral graphs, simplicial polyhedra an Apollonian networks

    * switches: -z will invert the filter to only output Hamiltonian graphs.
    * 			-r <int> determines the maximum recursive depth of the snake heuristic.
    * 			-k <int> determines the maximum recursive depth of pre-filtering. This has only been tested with polyhedral graphs.
	* 			-w will disable strict filtering at the recursive depth set by -k for complex pre-filtering.
    * 			-y will make Plantri output a table of longest path lengths.
	* 			-W will set the splithint, if modulo splitting below the default is desired.

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


//define global vars if needed

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#define FILTER plugin_filter
#define SUMMARY summary

//default switch behaviour

#define DEFAULT_TRUNCATE TRUE
#define DEFAULT_KEEP_HAM FALSE
#define DEFAULT_SNAKE_DEPTH 0
#define DEFAULT_MAKE_TABLE FALSE
#define DEFAULT_MAX_REC_DEPTH -1

//int paths_searched = 0;

uint16_t current_path_length=0;

uint64_t total_found=0;
uint64_t snake_found=0;
uint64_t graphs_checked=0;
uint64_t end_graphs_seen=0;
uint64_t end_graphs_checked=0;

#ifdef PATHCOUNT
uint64_t paths_expanded=0;
#endif

uint64_t path_length_table[256];
uint64_t search_for_path=0;

static int do_truncate = DEFAULT_TRUNCATE;
static int keep_ham = DEFAULT_KEEP_HAM;
static int snake_depth = DEFAULT_SNAKE_DEPTH;
static int make_table = DEFAULT_MAKE_TABLE;
static int splitflag = -1;
int64_t max_rec_depth=DEFAULT_MAX_REC_DEPTH;
#define PLUGIN_SWITCHES else if (arg[j] == 'z') keep_ham = TRUE; else if (arg[j] == 'r') snake_depth = getswitchvalue(arg,&j); else if (arg[j] == 'k') max_rec_depth = getswitchvalue(arg,&j); else if (arg[j] == 'w') do_truncate = FALSE; else if (arg[j] == 'y') make_table = TRUE;  else if (arg[j] == 'W') splitflag = getswitchvalue(arg,&j);

#define HAM_STR keep_ham ? "Hamiltonian" : "Non-Hamiltonian"

#define PLUGIN_INIT printf("Only saving %s graphs with recursive depth %d\n", HAM_STR, snake_depth); if (max_rec_depth != -1) {do_truncate ? printf("Truncating generation after %ld levels\n", max_rec_depth) : printf("Inferring hamiltonicity of output from filtering up to level %ld\n", max_rec_depth);}; if (splitflag != -1) splithint = splitflag; 

static void build_path(uint64_t *path_vertices, int *current_vertex, int *parents, int backwards) {
	int keep_adding = 1;
	while (keep_adding) {
		int found_next = 0;
		int next_vertex = (backwards) ? -1 : *current_vertex;
		EDGE *e, *elast;
		e = elast = firstedge[*current_vertex];
		do {
#ifdef PATHCOUNT
paths_expanded++;
#endif
			if (!(*path_vertices & (1ULL << e->end))) {
				found_next = 1;
				next_vertex = e->end;
				if (backwards) {parents[*current_vertex] = next_vertex; parents[next_vertex] = -1;} else {parents[next_vertex] = *current_vertex;}
				*path_vertices |= 1ULL << next_vertex;
				*current_vertex = next_vertex;
			}
			e = e->next;
		} while ((!found_next) && (e != elast));
		if (!found_next) keep_adding = 0;
	}
}

static int pivot(int* parents, int last_vertex, int pivot) {
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

static int snake_recurse(uint64_t *path_vertices, int first_vertex, int last_vertex, int* parents, int max_depth, int previous_pivot) {
	uint16_t length = __builtin_popcount(*path_vertices);
	if (length > current_path_length) current_path_length = length;
	
	EDGE *e, *elast;
	if ((*path_vertices + 1) >> maxnv) { //If path is on n vertices: check if ham
		e = elast = firstedge[last_vertex];
		do {
			if (first_vertex == e->end) return 1; // check if ham
			e = e->next;
		} while (e != elast);
	}
	if (!max_depth) return 0; // halt if at max depth
	
	e = elast = firstedge[last_vertex];
	do {
		// do the snake thing
		uint64_t curr_path = *path_vertices;
		int vertex_to_pivot = e->end;
		if ((vertex_to_pivot != parents[last_vertex]) && (vertex_to_pivot != previous_pivot)) { // no point pivoting with itself or the pivot in a previous iteration
			int pivot_parents[maxnv];
			memcpy(pivot_parents, parents, sizeof(int)*maxnv);
			int new_last = pivot(pivot_parents, last_vertex, vertex_to_pivot);
			build_path(&curr_path, &new_last, pivot_parents, FALSE); //attempt to grow
			int result = snake_recurse(&curr_path, first_vertex, new_last, pivot_parents, max_depth-1, vertex_to_pivot);
			if (result) {
				*path_vertices = curr_path;
				return 1;
			}
		}
		e = e->next;
	} while (e != elast);
	return 0;
}

static int is_hamiltonian(uint64_t path_vertices, int current_vertex) {// returns 1 if hamiltonian
	uint16_t length = __builtin_popcount(path_vertices);
	if (length > current_path_length) current_path_length = length;
	
	EDGE *e, *elast;
	e = elast = firstedge[current_vertex];
	
	do {
#ifdef PATHCOUNT
paths_expanded++;
#endif
		int vertex = e->end;
		if (((path_vertices + 1) >> maxnv) && vertex == 0) return 1; // current path is on n vertices and last vertex is adjacent with first vertex (i.e. first n bits are set and vertex == 0)
		if (!(path_vertices & (1ULL << vertex))) { // if edge 'vertex' is not in path (i.e. bit 'vertex' is not set)
			int recursive_result = is_hamiltonian(path_vertices | 1ULL << vertex, vertex);
			if (recursive_result) return 1;
		}
		e = e->next;
	} while (e != elast);
	
	return 0;
}

static int ham_snake() { //Returns 1 if hamiltonian, 0 otherwise
	graphs_checked++;
	
	if (snake_depth == 0) {
		uint64_t recursive_path = 1;
		int recursive_start_vertex = 0;
		return is_hamiltonian(recursive_path, recursive_start_vertex);
	}
	
	uint64_t path_vertices = 1;
	int current_vertex = 0;
	int parent[maxnv];
	parent[current_vertex] = -1;
	build_path(&path_vertices, &current_vertex, parent, FALSE); //FALSE marks pathing forwards
	int last_vertex = current_vertex;
	int first_vertex = 0;
	build_path(&path_vertices, &first_vertex, parent, TRUE); // TRUE marks pathing backwards
	int finally_ham = snake_recurse(&path_vertices, first_vertex, last_vertex, parent, snake_depth, first_vertex);
	if (finally_ham) {
		snake_found ++;
		total_found++;
		return 1;
	} else {
		uint64_t recursive_path = 1;
		int recursive_start_vertex = 0;
		int result = is_hamiltonian(recursive_path, recursive_start_vertex);
		if (result) total_found++;
		return result;
	}
}

uint8_t longest_path(uint64_t *path, uint8_t last_vertex, uint8_t curr_len) {
	if (maxnv == curr_len) {
		return curr_len;
	}
	uint64_t current_longest = 0ULL;
	uint8_t current_length = curr_len;
	EDGE *e, *elast;
	e = elast = firstedge[last_vertex];
	
	do {
		int vertex = e->end;
		if (!(*path & (1ULL << vertex))) {
			uint64_t next_path = *path + (1ULL << vertex);
			uint8_t next_path_len = longest_path(&next_path, vertex, curr_len + 1);
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

uint8_t get_longest_path(void) {
	search_for_path++;
	uint8_t path_length = 0;
	for (int i=0; i<(maxnv); i++) {
		uint64_t p = 1ULL << i;
		path_length = longest_path(&p, i, 1);
		if (path_length == maxnv) return path_length;
	}
	return path_length;
}

static int plugin_filter(int nbtot, int nbop, int doflip, int check_ham, int just_checked) {
#ifdef PATHCOUNT
paths_expanded=0;
uint64_t x = snake_found;
#endif
	
	current_path_length = 0;
	end_graphs_seen++;
	int result;
	if (just_checked) {
		result = !(check_ham ^ keep_ham);
	} else if (!check_ham) {
		result = !keep_ham;
	} else {
		end_graphs_checked++;
		result = !(ham_snake() ^ keep_ham);
	}
	
	if (make_table && result) path_length_table[get_longest_path()]++;
	if (end_graphs_seen % 10000000 == 0) {
		printf("Checked %ld graphs. %ld %s graph(s) found\n", end_graphs_seen, keep_ham? total_found : end_graphs_seen - end_graphs_checked, HAM_STR);
		fflush(stdout);
	}
#ifdef PATHCOUNT
printf("HAM: %s, SNAKE FOUND: %s, PATHS EXPANDED: %ld\n", result ? "NO" : "YES", x == snake_found ? "NO" : "YES", paths_expanded);
#endif
    return result;
}

static void summary() {
	printf("Graphs found using Snake: %ld\n", snake_found);
	printf("Checked %ld graphs\n", graphs_checked);
	printf("Seen %ld graphs on final filter\n", end_graphs_seen);
	printf("Checked %ld graphs on final filter\n", end_graphs_checked);
	if (make_table) {
		uint64_t checksum = 0;
		printf("Path Length Table (length is number of vertices):\nlen : Number of graphs with this length\n");
		for (uint16_t i=0; i < 256; i++) {
			if (path_length_table[i] > 0) {printf("%3d : %ld\n", i, path_length_table[i]); checksum += path_length_table[i];}
		}
		uint64_t out_count = oswitch ? totalout_op : totalout;
		if (checksum != out_count) printf("Checksum failed. %ld != %ld\n", checksum, out_count);
		printf("Reverted to exhastive path checking for %ld graphs\n", search_for_path);
		if (splitflag != -1) printf("Splitting level is %d.\n", splitlevel);
	}
}
