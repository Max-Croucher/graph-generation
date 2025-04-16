/* PLUGIN file to use with plantri.c
	Max Croucher, 2025
	This plugin is distributed at https://github.com/Max-Croucher/graph-generation
	Implements a variety of filtering methods to filter planar graphs by path lengths.

	This plugin is tested only with graphs which are (mostly) polyhedral! i.e. polyhedral graphs, simplicial polyhedra an Apollonian networks

    * switches: -z<int> will overwrite maximum path length to <int>. by default this limit is 2n/3.
    * 			-r <int> determines the maximum recursive depth of the snake heuristic.
    * 			-w will enable the Floyd-Warshall Algorithm.
	* 			-y <int> will filter graphs to output only graphs with diameter set with this flag.

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


float max_path_len;
uint64_t num_skipped=0;
uint64_t num_found=0;
uint64_t snake_found=0;
uint64_t triangs_seen=0;
uint64_t triangs_expanded=0;
uint64_t floyd_found=0;

uint64_t diameters[256];

static int snake_depth = 0;
static int zvalue = 0;
static int floyd_flag = 0;
static int do_diam = 0;
#define PLUGIN_SWITCHES else if (arg[j] == 'z') zvalue = getswitchvalue(arg,&j);\
else if (arg[j] == 'r') snake_depth = getswitchvalue(arg,&j);\
else if (arg[j] == 'w') floyd_flag = 1;\
else if (arg[j] == 'y') {do_diam = getswitchvalue(arg,&j); floyd_flag = 1;}


#define PLUGIN_INIT if (zvalue && do_diam) {fprintf(stderr, "Error: -z and -y may not be set together.\n"); exit(EXIT_FAILURE);}\
if (do_diam) printf("Only outputting graphs with diameter %d.\n", do_diam); else {if (zvalue) {max_path_len = zvalue; printf("Skipping all graphs with path length >= %.3f (overwritten by flag)\n", max_path_len);} else {max_path_len = ((float)(2 * maxnv) / 3); printf("Skipping all graphs with path length >= %.3f.\n", max_path_len);}}

static void build_path(uint64_t *path_vertices, int *current_vertex, int *parents, int backwards) {
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
	int intermediary_vertices[nv];
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
		int vertex_to_pivot = e->end;
		if ((vertex_to_pivot != parents[last_vertex]) && (vertex_to_pivot != previous_pivot)) { // no point pivoting with itself or the pivot in a previous iteration
			int pivot_parents[nv];
			memcpy(pivot_parents, parents, sizeof(int)*nv);
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

static int has_long_path(uint64_t path_vertices, int current_vertex, int current_length) {
	if (current_length-1 >= max_path_len) { //check if current path is too long. current_length counts the NUMBER_OF_VERTICES, and thus current_length-1 countes the number of edges
		return 1;
	}
	EDGE *e, *elast;
	e = elast = firstedge[current_vertex];
	
	do {
		int vertex = e->end;
		if (!(path_vertices & (1ULL << vertex))) { // if edge 'vertex' is not in path (i.e. bit 'vertex' is not set)
			int recursive_result = has_long_path(path_vertices | (1ULL << vertex), vertex, current_length+1);
			if (recursive_result) return 1;
		}
		e = e->next;
	} while (e != elast);
	return 0;
}

static int long_path_snake() {

	if (snake_depth == 0) {
		for (int i = 0; i<nv; i++) {
			uint64_t recursive_path = 1ULL << i;
			if (has_long_path(recursive_path, i, 1)) return 1;
		}
		return 0;
	}
	
	uint64_t path_vertices = 1;
	int current_vertex = 0;
	int parent[nv];
	parent[current_vertex] = -1;
	build_path(&path_vertices, &current_vertex, parent, FALSE); //FALSE marks pathing forwards
	int last_vertex = current_vertex;
	int first_vertex = 0;
	build_path(&path_vertices, &first_vertex, parent, TRUE); // TRUE marks pathing backwards
	int finally_ham = snake_recurse(&path_vertices, first_vertex, last_vertex, parent, snake_depth, first_vertex);
	if (finally_ham) {
		snake_found++;
		return 1;
	} else {
		for (int i = 0; i<((float)nv-1); i++) {
			uint64_t recursive_path = 1ULL << i;
			if (has_long_path(recursive_path, i, 1)) return 1;
		}
		return 0;
	}
}

static void print_matrix(int8_t matrix[][maxnv]) {
	for (int v=0; v<maxnv; v++) {
		for (int u=0; u<maxnv; u++) {
			switch (matrix[v][u]) {
				case -1:
					printf(" x ");
					break;
				case 0:
					printf(" 0 ");
					break;
				default:
					printf(" %-2d", matrix[v][u]);
			}
		}
		printf("\n");
	}
}

static int floyd() {
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
				if ((distance_matrix[i][k] >= 0 && distance_matrix[k][j] >= 0) && ((distance_matrix[i][j] < 0) || (distance_matrix[i][j] > distance_matrix[i][k] + distance_matrix[k][j]))) {
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

static int is_permissable() {
	int diameter;
	if (floyd_flag) {
		diameter = floyd();
		diameters[diameter]++;
		
		if (do_diam) return diameter == do_diam;

		if (3*diameter - 2 >= max_path_len) {//Path formed by exploiting 3-connectivity, measured by number of edges.
			floyd_found++;
			return 0;
		}
	}

	return !long_path_snake();
	
}

static int plugin_filter(int nbtot, int nbop, int doflip) {
	//paths_searched = 0;
    if (!is_permissable()) {
		num_skipped++;
		if ((num_skipped + num_found) % 10000000 == 0) {
			printf("Checked %ld graphs. %ld counterexamples found\n", num_skipped + num_found, num_found);
			fflush(stdout);
		}
		//printf("FOUND %d\n", paths_searched);
		return 0;
		}
    num_found++;
    if ((num_skipped + num_found) % 10000000 == 0) {
		printf("Checked %ld graphs. %ld counterexamples found\n", num_skipped + num_found, num_found);
		fflush(stdout);
	}
	//printf("NOT FOUND %d\n", paths_searched);
	//printf("%ld\n", num_skipped + num_found);
    return 1;
}

static void summary() {
	if (floyd_flag) printf("Graphs found using Floyd: %ld\n", floyd_found);
	printf("Graphs found using Snake: %ld\n", snake_found);
	printf("Graphs skipped: %ld\n", num_skipped);
	printf("Graphs saved: %ld\n", num_found);
	if (floyd_flag) {
		printf("Diameters of graphs:\nd   : Number of graphs with this diameter:\n");
		for (uint16_t i=0; i < maxnv; i++) {
			if (diameters[i] > 0) printf("%3d : %ld\n", i, diameters[i]);
		}
	}
}
