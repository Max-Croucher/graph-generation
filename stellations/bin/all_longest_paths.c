/* all_longest_paths.c
	Max Croucher, 2025
	This program is written to iterate over the graphs in a planar code file (typically .pl),
	as output by the planar graph generation tool Plantri. For each graph in the
	given file, this program will print to stdout a textual representation of the
	graph as an adjacency matrix, and then print every longest path in this graph.
	Note that a particular Longest path will be written twice - once for each direction,
	so further processing of the output is required to count the number of longest
	paths a graph may have.
	For example, if the given file contains only the complete graph K4
	(as is given in example_graphs/k4.pl), the following will be output:

	graph = {
	  0: [0, 1, 1, 1],
	  1: [1, 0, 1, 1],
	  2: [1, 1, 0, 1],
	  3: [1, 1, 1, 0]
	}
	0 - 1 - 2 - 3 
	0 - 1 - 3 - 2 
	0 - 2 - 1 - 3 
	0 - 2 - 3 - 1 
	0 - 3 - 1 - 2 
	0 - 3 - 2 - 1 
	1 - 0 - 2 - 3 
	1 - 0 - 3 - 2 
	1 - 2 - 0 - 3 
	1 - 2 - 3 - 0 
	1 - 3 - 0 - 2 
	1 - 3 - 2 - 0 
	2 - 0 - 1 - 3 
	2 - 0 - 3 - 1 
	2 - 1 - 0 - 3 
	2 - 1 - 3 - 0 
	2 - 3 - 0 - 1 
	2 - 3 - 1 - 0 
	3 - 0 - 1 - 2 
	3 - 0 - 2 - 1 
	3 - 1 - 0 - 2 
	3 - 1 - 2 - 0 
	3 - 2 - 0 - 1 
	3 - 2 - 1 - 0 

	Optional Parameters:
	In addition to the filename, the program will accept an optional integer parameter,
	which indicates how many graphs should be output.

	Copyright for this plugin is held by the author
	Max Croucher, University of Canterbury, mpccroucher@gmail.com
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h> 
#include <string.h>

#define HEADER_SIZE 15
#define MAX_ORDER 63
int64_t max_graphs;
uint64_t graphs_seen = 0ULL;
uint64_t path_lengths[MAX_ORDER] = {0ULL};
uint8_t all_path_vertices[MAX_ORDER] = {0ULL};
int current_remaining = 0;
int target_vertex = 0;
int graph_index = 0;
uint64_t filter_count = 0ULL;
uint64_t filter_checksum = 0ULL;
uint64_t num_paths_found = 0ULL;
uint64_t vertex_mask = 0;


int check_header(unsigned char *header) {
	/* Verifies the file header */
	unsigned char true_header[HEADER_SIZE] = ">>planar_code<<";
	for (int i = 0; i < HEADER_SIZE; i++) {
		if (true_header[i] != header[i]) {
			return 0;
		}
	}
	return 1;
}


void print_graph(uint8_t *graph, uint8_t order) {
	/* Prints the adjacency matrix of the graph */
	printf("graph = {\n");
	for (uint8_t i=0; i<order; i++) {
		printf("  %d: [", i);
		for (uint8_t j=0; j<order; j++) {
			if (j > 0) printf(", ");
			printf("%d", graph[i*order+j]);
		}
		if (i < order-1) {printf("],\n");} else {printf("]\n");}
	}
	printf("}\n");
}


uint8_t longest_path(uint8_t *graph, uint8_t order, uint64_t *path, uint8_t last_vertex, uint8_t curr_len) {
	/* Recursively finds the length of the longest path in the current graph */
	if (order == curr_len) {
		return curr_len;
	}
	uint64_t current_longest = 0ULL;
	uint8_t current_length = curr_len;
	for (uint8_t i=0; i<order; i++) {
		if ((graph[last_vertex*order+i]) && !(*path & (1ULL << i))) {
			uint64_t next_path = *path + (1ULL << i);
			uint8_t next_path_len = longest_path((uint8_t *)graph, order, &next_path, i, curr_len + 1);
			if (next_path_len == order) {
				*path = next_path;
				return next_path_len;
			}
			if (next_path_len > current_length) {
				current_longest = next_path;
				current_length = next_path_len;
			}
		}
	}
	*path = current_longest;
	return current_length;
}


int get_path_length(uint8_t order, uint8_t* graph) {
	/* Determines the length of the longest path in the current graph */
	uint64_t current_longest = 0ULL;
	uint8_t current_length = 0;
	for (uint8_t v=0; v<order; v++) {
		uint64_t path = 1ULL << v;
		uint8_t path_len = longest_path(graph, order, &path, v, 1);
		if (path_len > current_length) {
			current_length = path_len;
			current_longest = path;
		}
		if (current_length == order) break;
	}
	return current_length;
}


uint8_t find_longest_paths_recurse(uint8_t *graph, uint8_t order, uint64_t *path, uint8_t last_vertex, uint8_t curr_len, uint8_t target_length) {
	/* Recursively builds paths, printing any that are of length 'target_length' */
	if (target_length == curr_len) {
		num_paths_found++;
		for (int i=0; i<target_length; i++) {
			printf("%d ", all_path_vertices[i]);
			if (i+1 < target_length) printf("- ");
		}
		printf("\n");
		vertex_mask |= *path;
	}
	uint64_t current_longest = 0ULL;
	uint8_t current_length = curr_len;
	for (uint8_t i=0; i<order; i++) {
		if ((graph[last_vertex*order+i]) && !(*path & (1ULL << i))) {
			uint64_t next_path = *path + (1ULL << i);
			all_path_vertices[curr_len] = i;
			uint8_t next_path_len = find_longest_paths_recurse((uint8_t *)graph, order, &next_path, i, curr_len + 1, target_length);
			if (next_path_len > current_length) {
				current_longest = next_path;
				current_length = next_path_len;
			}
		}
	}
	*path = current_longest;
	return current_length;
}


void find_longest_paths(uint8_t order, uint8_t* graph, uint8_t pathlen) {
	/* Triggers the recursive procedure for finding longest paths, starting from every vertex. */
	vertex_mask = 0ULL;
	num_paths_found = 0ULL;
	filter_checksum++;
	for (uint8_t v=0; v<order; v++) {
		if (all_path_vertices || !((vertex_mask + 1) >> order)) {
			uint64_t path = 1ULL << v;
			all_path_vertices[0] = v;
			find_longest_paths_recurse(graph, order, &path, v, 1, pathlen);
		}
	}
}


int scan_file(char filename[]) {
	/* Read and process a file */
	FILE *f;
	f = fopen(filename, "rb");
	if (f == NULL) {
		return 1;
	}
	
	unsigned char header[HEADER_SIZE];
	fread(header, HEADER_SIZE, 1, f);
	if (!check_header(header)) return 2;
	
	int c;
	while (!(max_graphs == 0) & ((c = fgetc(f)) != EOF)) {
		uint8_t order = (uint8_t)c;
		
		uint8_t graph[order][order];
		for (uint8_t i=0; i<order; i++) {
			for (uint8_t j=0; j<order; j++) {
				graph[i][j] = 0;
			}
			int go = 1;
			do {
				uint8_t c = fgetc(f);
				if (c) {
					graph[i][c-1] = 1;
				} else {
					go = 0;
				}
			} while (go);
		}
		
		uint8_t length = get_path_length(order, *graph);

		graphs_seen++;
		path_lengths[length]++;
		print_graph(*graph, order);
		find_longest_paths(order, *graph, length);
		
		if (max_graphs > 0) max_graphs--;
	}
	fclose(f);
	return 0;
}


int main(int argc, char *argv[]) {
	/* Main function that reads parameters and starts file processing */
	if (argc != 2 && argc != 3) {
		printf("Usage: %s <pl filename> [max_graphs]\n", argv[0]);
		exit(1);
	}
	
	max_graphs = (argc == 3) ? atoll(argv[2]) : -1;

	int code = scan_file(argv[1]);
	if (code == 1) {
		printf("Cannot open the file.\n");
		exit(1);
	} else if (code == 2) {
		printf("Invalid pl file.\n");
		exit(1);
	}
	return 0;
}
