/* all_longest_paths.c
	Max Croucher, 2025
	This program is written to iterate over the graphs in a planar code file (typically .pl),
	as output by the planar graph generation tool Plantri. For each graph in the
	given file, the program will compute the longest paths in this graph that can
	be generated with a given vertex as an endpoint. This tool also supports a range
	of constraints, including the blacklisting of vertices, or requirement that paths
	contain particular vertices or paths. More information on these constraints
	is available in the readme.
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
#define MAX_ORDER 63 // path membership is tracked using the bits of a uint64_t.
// this could be changed to an array, which would increase maximum order to the largest
// graph size supported by planar code files, which is 255 ish.

uint64_t graphs_seen = 0ULL;
uint64_t paths_seen = 0ULL;
uint64_t blacklisted_vertices = 0ULL;
uint64_t required_vertices = 0ULL;
uint64_t required_path_endpoints = 0ULL;
uint64_t required_path_internal_vertices = 0ULL;
uint64_t non_traversable_vertices = 0ULL;
uint8_t required_paths[MAX_ORDER][MAX_ORDER] = {};
uint8_t required_path_lengths[MAX_ORDER] = {};
uint8_t path_order[MAX_ORDER] = {0};
uint8_t halt_on_trace = 0;
uint8_t global_longest = 0;
uint8_t print_verbose = 0;
uint8_t do_used_edges = 0;
uint8_t edge_membership[MAX_ORDER][MAX_ORDER] = {0};
uint8_t graph[MAX_ORDER][MAX_ORDER] = {0};
#define ADMISSIBLE(path_arr) (!(required_vertices&~(path_arr)))


inline static int check_header(unsigned char *header) {
	/* Verify the header of a planar code file */
	unsigned char true_header[HEADER_SIZE] = ">>planar_code<<";
	for (int i = 0; i < HEADER_SIZE; i++) {
		if (true_header[i] != header[i]) {
			return 0;
		}
	}
	return 1;
}


static void print_graph(uint8_t order) {
	/* print an adjacency matrix for a graph */
	printf("graph = {\n");
	for (uint8_t i=0; i<order; i++) {
		printf("  %2d: [", i);
		for (uint8_t j=0; j<order; j++) {
			if (j > 0) printf(", ");
			printf("%d", graph[i][j]);
		}
		if (i < order-1) {printf("],\n");} else {printf("]\n");}
	}
	printf("}\n");
}


inline static void print_used_edges(uint8_t order) {
	/* Print the set of edges used by any of the longest paths generated */
	printf("Used Edges: \n ");
	for (int v=0; v<order; v+=5) {
		printf("   %2d", v);
	}
	printf("\n ");
	for (int v=0; v<order; v+=5) {
		printf("    |");
	}
	printf("\n");
	for (int u=0; u<order; u++) {
		if (!(u%5)) {
			printf("%2d - ", u);
		} else {
			printf("     ");
		}
		for (int v=0; v<order; v++) {
			if (v<u) {
				printf(" ");
			} else if (blacklisted_vertices & ((1ULL << u) | (1ULL << v))) {
				printf("*");
			} else if (!graph[u][v]) {
				printf(".");
			} else {
				printf("%d", edge_membership[u][v]);
			}
		}
		printf("\n");
	}
}


inline static void print_path(uint64_t path_vertices, uint8_t order) {
	/* print useful information about a particular path */
	uint8_t length = __builtin_popcount(path_vertices);
	printf("\n\nCurrent length/Longest valid path length: %d/%d\n", length-1, global_longest-1);
	printf("Current Path: %d", path_order[0]);
	for (int i=1; i<length; i++) {
		printf("-%d", path_order[i]);
	}
	printf("\nAdmissible? %s\n", ADMISSIBLE(path_vertices) ? "Yes" : "No");
	printf("Path Membership: ");
	for (int i=0; i<order; i++) {
		if (blacklisted_vertices & (1ULL << i)) {
			printf("*");
		} else if (path_vertices & (1ULL << i)) {
			printf("1");
		} else {
			printf("0");
		}
	}
	printf("\n                 |");
	printf("%*c\n", order-1, '|');
	printf("             v = 0");
	printf("%*d\n", order-1, order-1);
}

inline static void update_edge_membership(uint8_t curr_len) {
	/* add an encountered edge to the path membership table */
	for (int i=0; i<curr_len-1; i++) {
		if (path_order[i] < path_order[i+1]) {
			edge_membership[path_order[i]][path_order[i+1]] = 1;
		} else {
			edge_membership[path_order[i+1]][path_order[i]] = 1;
		}
	}
}


static uint8_t longest_path(uint8_t order, uint64_t *path, uint8_t last_vertex, uint8_t curr_len) {
	/* recursively find the longest path with the given constraints */
	paths_seen++;
	if (ADMISSIBLE(*path) && curr_len > global_longest) { // does a path meet all the requirements
		global_longest = curr_len;
		for (int u=0; u<order; u++) {
			for (int v=u; v<order; v++) {
				edge_membership[u][v] = 0;
			}
		}
		if (do_used_edges) update_edge_membership(curr_len);
		print_path(*path, order);
	} else if (ADMISSIBLE(*path) && curr_len == global_longest) {
		if (do_used_edges) update_edge_membership(curr_len);
		if (print_verbose) print_path(*path, order);
	}
	if (do_used_edges && !(paths_seen%100000000)) print_used_edges(order); // periodically report
	if (halt_on_trace && order == curr_len) {
		return curr_len;
	}
	uint64_t current_longest = 0ULL;
	uint8_t current_length = curr_len;
	for (uint8_t i=0; i<order; i++) {
		if ((graph[last_vertex][i]) && !(*path & (1ULL << i)) && !(non_traversable_vertices & (1ULL << i))) {
			uint64_t next_path = *path + (1ULL << i);
			path_order[curr_len] = i;
			uint8_t next_path_len;
			if (required_path_endpoints & (1ULL << i)) {
				for (int j=1; j<required_path_lengths[i]; j++) {
					next_path |= (1ULL << required_paths[i][j]);
					path_order[curr_len+j] = required_paths[i][j];
				}
				next_path_len = longest_path(order, &next_path, required_paths[i][required_path_lengths[i]-1], curr_len + required_path_lengths[i]);
			} else {
				next_path_len = longest_path(order, &next_path, i, curr_len + 1);
			}
			if (halt_on_trace && next_path_len == order) {
				*path = next_path;
				return next_path_len;
			}
			if (ADMISSIBLE(*path) && next_path_len > current_length) {
				current_longest = next_path;
				current_length = next_path_len;
			}
		}
	}
	*path = current_longest;
	return current_length;
}


inline static int scan_file(char filename[], int start_vertex) {
	/* read a file and process every graph */
	FILE *f;
	f = fopen(filename, "rb");
	if (f == NULL) {
		return 1;
	}
	
	unsigned char header[HEADER_SIZE];
	fread(header, HEADER_SIZE, 1, f);
	if (!check_header(header)) return 2;
	
	int c;
	while ((c = fgetc(f)) != EOF) {
		uint8_t order = (uint8_t)c;

		if (blacklisted_vertices >> order) {printf("Error: Blacklisted vertices must be less than the order of graph %ld, which is %d.\n", graphs_seen, order); exit(EXIT_FAILURE);}
		if ((required_path_endpoints >> order) || (required_path_endpoints >> order)) {printf("Error: Vertices in required paths must be less than the order of graph %ld, which is %d.\n", graphs_seen, order); exit(EXIT_FAILURE);}
		if (order > MAX_ORDER) {printf("Error: Graph is too large. Maximum order is %d.\n", MAX_ORDER); exit(EXIT_FAILURE);}
		if (order <= start_vertex) {printf("Error: Start vertex %d must be less than the order of graph %ld, which is %d.\n", start_vertex, graphs_seen, order); exit(EXIT_FAILURE);}

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
		for (int u=0; u<order; u++) {
			for (int v=u; v<order; v++) {
				edge_membership[u][v] = 0;
			}
		}

		if (required_path_endpoints) {
			for (int i=0; i<order; i++) {
				if (required_path_endpoints & (1ULL << i)) {
					for (int j=0; j<required_path_lengths[i]-1; j++) {
						if (!graph[required_paths[i][j]][required_paths[i][j+1]]) {printf("Error: required edge (%d,%d) is not present in graph %ld.\n", required_paths[i][j], required_paths[i][j+1], graphs_seen); exit(EXIT_FAILURE);}
					}
				}
			}
		}

		uint8_t length;

		paths_seen = 0ULL;
		uint64_t path = 1ULL << start_vertex;
		path_order[0] = start_vertex;
		global_longest = 1;

		if (required_path_endpoints & (1ULL << start_vertex)) {
			for (int j=1; j<required_path_lengths[start_vertex]; j++) {
				path |= (1ULL << required_paths[start_vertex][j]);
				path_order[j] = required_paths[start_vertex][j];
			}
			global_longest = required_path_lengths[start_vertex];
			length = longest_path(order, &path, required_paths[start_vertex][required_path_lengths[start_vertex]-1], required_path_lengths[start_vertex]);
		} else {
			length = longest_path(order, &path, start_vertex, 1);
		}
		if (do_used_edges) print_used_edges(order);
		printf("(Graph %ld): The longest path with vertex %d as endpoint has length %d.\n", graphs_seen, start_vertex, global_longest-1);
		graphs_seen++;
	}
	fclose(f);
	return 0;
}


inline static void help_msg(char* name) {
	/* print help */
	printf("Usage: %s <pl-filename> <start-vertex>\n\nOptional Arguments:\n\n", name);
	printf("--help\n    Prints this message and exits.\n\n");
	printf("--halt-on-trace\n    Halts the graph generation procedure for a graph if a trace is found.\n\n");
	printf("--print\n    Prints a generated path if it is equally as long as the longest path found so far,\n");
	printf("    instead of only printing generated paths if they're longer.\n\n");
	printf("--used-edges\n    Prints an adjacency matrix of all edges used by any longest path.\n");
	printf("    Vertex pairs that include a blacklisted vertex are displayed as an asterisk \"*\".\n");
	printf("    Vertex pairs that do not correspond to an edge are displayed with a period \".\".\n");
	printf("    Vertex pairs that correspond to an edge are displayed as a \"1\" or \"0\" based on\n");
	printf("    whether or not there exists a longest path that uses this edge.\n\n");
	printf("--blacklist-vertices [sequence of vertices]\n");
	printf("    Specifies a sequence of vertices that cannot be present in a longest path.\n\n");
	printf("--require-vertices [sequence of vertices]\n");
	printf("    Specifies a sequence of vertices that must be present in a longest path.\n\n");
	printf("--require-edges [sequence of hyphen-delimetered sequences of vertices]\n");
	printf("    Specifies a sequence of paths that must be present in a longest path.\n");
	printf("    These required paths must be vertex disjoint and cannot involve blacklisted vertices.\n");
	exit(1);
}


inline static void parse_required_paths(char* path_string) {
	/* process an input string of paths */
	uint8_t length = 0;
	uint8_t current_path[MAX_ORDER];
	char* token = strtok(path_string, "-");
	while (token)
	{
		uint8_t v = atoi(token);
		if (v >= MAX_ORDER) {printf("Error: Vertex in required path is too large. Maximum is %d.\n", MAX_ORDER+1); exit(EXIT_FAILURE);}
		if ((1ULL << v) & (required_path_endpoints | required_path_internal_vertices)) {printf("Error: Required paths must be vertex disjoint.\n"); exit(EXIT_FAILURE);}
		if ((1ULL << v) & blacklisted_vertices) {printf("Error: Required paths cannot use blacklisted vertices.\n"); exit(EXIT_FAILURE);}
		current_path[length] = v;
		length++;
		if (length >= MAX_ORDER) {printf("Error: required path cannot be longer than the maximum graph order %d.\n", MAX_ORDER); exit(EXIT_FAILURE);}
		token = strtok(NULL, "-");
	}
	if (length == 1) {
		printf("WARNING: Ignoring required path of length 0.\n");
	} else {
		required_path_endpoints |= 1ULL << current_path[0];
		required_path_endpoints |= 1ULL << current_path[length-1];
		required_path_lengths[current_path[0]] = length;
		required_path_lengths[current_path[length-1]] = length;
		for (int j=0; j<length; j++) {
			required_paths[current_path[0]][j] = current_path[j];
			required_paths[current_path[length-1]][j] = current_path[length-1-j];
			if (j && (j<length-1)) required_path_internal_vertices |= 1ULL << current_path[j];
		}
	}
}


int main(int argc, char *argv[]) {
	/* read optional arguments and set up program */
	uint8_t input_state = 0;
	if (argc < 3) help_msg(argv[0]);
	for (int i = 3; i < argc; i++) {
		if (!strcmp(argv[i], "--help")) {
			input_state = 0;
			help_msg(argv[0]);
		}
		else if (!strcmp(argv[i], "--halt-on-trace")) {
			input_state = 0;
			halt_on_trace = 1;
		} else if (!strcmp(argv[i], "--print")) {
			input_state = 0;
			print_verbose = 1;
		} else if (!strcmp(argv[i], "--used-edges")) {
			input_state = 0;
			do_used_edges = 1;
		} else if (!strcmp(argv[i], "--blacklist-vertices")) {
			input_state = 1;
		} else if (!strcmp(argv[i], "--require-vertices")) {
			input_state = 2;
		} else if (!strcmp(argv[i], "--require-paths")) {
			input_state = 3;
		} else if (input_state == 1) {
			int v = atoi(argv[i]);
			if (v >= MAX_ORDER) {printf("Error: Blacklisted vertex is too large. Maximum is %d.\n", MAX_ORDER); exit(EXIT_FAILURE);}
			blacklisted_vertices |= (1ULL << v);
		} else if (input_state == 2) {
			int v = atoi(argv[i]);
			if (v >= MAX_ORDER) {printf("Error: Required vertex is too large. Maximum is %d.\n", MAX_ORDER); exit(EXIT_FAILURE);}
			required_vertices |= (1ULL << v);
		} else if (input_state == 3) {
			parse_required_paths(argv[i]);
		}
	}

	int start_vertex = atoi(argv[2]);
	if (start_vertex >= MAX_ORDER) {printf("Error: Starting vertex is too large. Maximum is %d.\n", MAX_ORDER+1); exit(EXIT_FAILURE);}
	if (blacklisted_vertices & (1ULL << start_vertex)) {printf("Error: Starting vertex cannot be blacklisted.\n"); exit(EXIT_FAILURE);}
	if (required_path_internal_vertices & (1ULL << start_vertex)) {printf("Error: Starting vertex cannot be an internal vertex of a required path.\n"); exit(EXIT_FAILURE);}

	required_vertices |= required_path_endpoints | required_path_internal_vertices;
	non_traversable_vertices = blacklisted_vertices | required_path_internal_vertices;

	int code = scan_file(argv[1], start_vertex);
	if (code == 1) {
		printf("Cannot open the file.\n");
		exit(1);
	} else if (code == 2) {
		printf("Invalid pl file.\n");
		exit(1);
	}
	
	return 0;
}
