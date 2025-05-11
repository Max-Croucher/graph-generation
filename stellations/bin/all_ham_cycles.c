/* all_ham_cycles.c
	Max Croucher, 2025
	This program is written to iterate over the graphs in a planar code file (typically .pl),
	as output by the planar graph generation tool Plantri. For each graph in the
	given file, this program will print to stdout a textual representation of the
	graph as an adjacency matrix, and then print every Hamiltonian cycle in this
	graph, with vertex 0 as the first vertex. Note that a particular Hamiltonian
	cycle will be written in both directions, so further processing of the output
	is required to count the number of Hamiltonian cycles a graph may have.
	For example, if the given file contains only the complete graph K4
	(as is given in example_graphs/k4.pl), the following will be output:

	graph = {
      0: [0, 1, 1, 1],
      1: [1, 0, 1, 1],
      2: [1, 1, 0, 1],
      3: [1, 1, 1, 0]
	}
	3 - 2 - 1 - 0
	2 - 3 - 1 - 0
	3 - 1 - 2 - 0
	1 - 3 - 2 - 0
	2 - 1 - 3 - 0
	1 - 2 - 3 - 0

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
uint64_t graphs_seen = 0ULL;
uint64_t path_lengths[MAX_ORDER] = {0ULL};
uint8_t order;
int64_t max_graphs;


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


void print_graph(uint8_t *graph) {
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


void print_path(int *parents, int end_vertex) {
	/* Print the current path*/
	printf("%d", end_vertex);
	do {
		end_vertex = parents[end_vertex];
		printf(" - %d", end_vertex);
	} while (end_vertex);
	printf("\n");
}


int get_ham_cycles(uint8_t *graph, int current_vertex, uint64_t path_vertices, int *parents, int current_length) {
	/* Recursively find Hamiltonian cycles in the current graph*/
	if (current_length == order) {
		return graph[order*current_vertex];
	}
	for (int i = 0; i < order; i++) {
		if (!(path_vertices & (1ULL << i)) && graph[order*current_vertex+i]) {
			uint64_t new_path = path_vertices + (1ULL << i);
			parents[i] = current_vertex;
			if (get_ham_cycles(graph, i, new_path, parents, current_length+1)) print_path(parents, i);
		}
	}
	return 0;
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
	while (!(max_graphs == 0) && (c = fgetc(f)) != EOF) {

		order = (uint8_t)c;
		
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
		
		print_graph((uint8_t *)graph);
		
		uint64_t vertices = 1ULL;
		int parents[order];
		get_ham_cycles(*graph, 0, vertices, parents, 1);
		if (max_graphs > 0) max_graphs--;
		
	}
	fclose(f);
	return 0;
}


int main(int argc, char *argv[]) {
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
