""" Classes and methods for storing, handling, reading and writing planar graphs.
Max Croucher, February 2024
"""

from copy import deepcopy
from itertools import combinations
from random import random, shuffle
import numpy as np
import matplotlib.pyplot as plt
from draw_pl import draw
from read_pl import read_graphs


class AdjacencyList():
    """ Defines an AdjacencyList Class, to hold the structure of a Graph Object
    """
    def __init__(self):
        self._adj = dict()


    def __str__(self):
        return f"{type(self).__name__}: {self._adj}"


    def __deepcopy__(self, _):
        new = type(self)()
        new._adj = deepcopy(self._adj)
        return new


    def get_adjacency(self):
        """Return a copy of the adjacency list"""
        return deepcopy(self._adj)


    @property
    def vertices(self):
        """Return a list containing the vertices"""
        return list(self._adj.keys())


    @property
    def edges(self):
        """Return a list containing the edges"""
        e = []
        for v in self.vertices:
            for u in self.adjacent_vertices(v):
                if [u, v] not in e:
                    e.append([u, v])
        return e


    @property
    def order(self):
        """Return the number of vertices"""
        return len(self.vertices)


    @property
    def size(self):
        """Return the number of edges"""
        return len(self.edges)


    def import_dict(self, adj):
        """Set the adjacency list"""
        self._adj = deepcopy(adj)


    def degree(self, v):
        """Return the degree of a given vertex"""
        return len(self.adjacent_vertices(v))


    def min_degree(self):
        """compute the minimum degree of a vertex"""
        return min(self.degree(v) for v in self.vertices)


    def valid_vertex_name(self):
        """Find the smallest positive integer that is not a vertex name"""
        for i in range(1, self.order+2):
            if i not in self.vertices:
                return i
        raise ValueError("Unable to name a new vertex.")


    def add_vertex(self, v=None):
        """Add a vertex to the AdjacencyList"""
        if v is None:
            v = self.valid_vertex_name()
        if v not in self.vertices:
            self._adj[v] = []
        return v


    def adjacent_vertices(self, v):
        """Return a list of all vertices adjacent to a given list"""
        return list(self._adj[v])


    def add_edge(self, u, v=None):
        """Add an edge to the AdjacencyList, adding endpoints as vertices if not present"""
        if v is None:
            u, v = u
        if u == v:
            raise ValueError("Cannot add an edge from a vertex to itself.")
        self.add_vertex(u)
        self.add_vertex(v)
        if u not in self._adj[v]:
            self._adj[v].append(u)
        if v not in self._adj[u]:
            self._adj[u].append(v)


    def is_adjacent(self, u, v=None):
        """Check if two vertices are adjacent"""
        if v is None:
            u, v = u
        return u in self._adj[v]


    def delete_vertex(self, v):
        """Deletes a vertex from an AdjacencyList, including any incident edges"""
        for u in self.adjacent_vertices(v):
            self.delete_edge(u, v)
        self._adj.pop(v)


    def delete_edge(self, u, v=None):
        """Deletes an edge from an AdjacencyList"""
        if v is None:
            u, v = u
        if not self.is_adjacent(u, v):
            raise ValueError(f"There is no edge between {u} and {v}.")
        self._adj[u].remove(v)
        self._adj[v].remove(u)


    def rename_vertex(self, old_name, new_name):
        """Rename a vertex, adjusting both inbound and outbound edges"""
        if old_name not in self.vertices:
            raise KeyError(f"Vertex {old_name} is not in the graph")
        if new_name in self.vertices:
            raise KeyError(f"Vertex {new_name} is an existing vertex")
        self._adj[new_name] = self._adj.pop(old_name, None)
        for adjacent_vertex in self._adj[new_name]:
            i = self._adj[adjacent_vertex].index(old_name)
            self._adj[adjacent_vertex] = self._adj[adjacent_vertex][:i]+[new_name]+self._adj[adjacent_vertex][i+1:]


    def permute(self, ordering=None):
        """Permute the names of vertices. If no ordering is specified, the vertex
           names are shuffled at random"""
        if ordering is None:
            ordering = deepcopy(self.vertices)
            shuffle(ordering)
        else:
            if sorted(self.vertices) != sorted(ordering):
                raise ValueError("ordering must be a permutation of the vertices")
        new_adjacency = {}
        for i, v in enumerate(self.vertices):
            new_adjacency[ordering[i]] = deepcopy(self._adj[v])
            for j in range(len(new_adjacency[ordering[i]])):
                new_adjacency[ordering[i]][j] = ordering[self.vertices.index(self._adj[v][j])]

        new_keys = list(new_adjacency.keys())
        shuffle(new_keys)
        new_graph = Graph()
        new_graph.import_dict({x: new_adjacency[x] for x in new_keys})
        return new_graph


class Graph(AdjacencyList):
    """ Defines a Graph Class, to represent a graph and hold useful methods
        Inherits from AdjacencyList
    """
    def __init__(self, vertices=None, edges=None):
        """Initialise a Graph object"""
        super().__init__()
        if vertices is not None:
            for v in vertices:
                self.add_vertex(v)
        if edges is not None:
            for e in edges:
                self.add_edge(e)


    def reachable_vertices(self, v, target=None):
        """List all vertices that lie in the same component as a given vertex.
            Optional argument target halts the search if a given vertex is reached"""
        found = []
        expanded = [v]
        while len(expanded) != 0:
            current = expanded.pop()
            found.append(current)
            if target is not None and target == current:
                return found
            new_vertices = self.adjacent_vertices(current)
            for u in new_vertices:
                if u not in found and u not in expanded:
                    expanded.append(u)
        return found


    @property
    def components(self):
        """List all components of a graph, including singletons"""
        comps = []
        visited = set()
        for v in self.vertices:
            if v not in visited:
                current = self.reachable_vertices(v)
                visited.update(current)
                comps.append(current)
        return comps


    def is_reachable(self, u, v):
        """checks if there exists a path between two vertices"""
        return u in self.reachable_vertices(v, target=u)


    def is_face(self, vertices, direction='cw'):
        """Determines if a set of vertices is a face"""
        if len(vertices) < 3:
            raise ValueError("Graph.is_face requires at least 3 positional arguments")

        if direction.lower() == 'ccw':
            offset = 1
        elif direction.lower() == 'cw':
            offset = -1
        else:
            raise ValueError("Optional parameter offset must be cw or ccw")

        for i in range(len(vertices)):
            u = vertices[i]
            v = vertices[(i+1)%len(vertices)]
            w = vertices[(i+2)%len(vertices)]
            try:
                if self.adjacent_vertices(v)[(self.adjacent_vertices(v).index(u) + offset) % self.degree(v)] != w:
                    return False
            except ValueError:
                return False
        return True


    def draw(self, *mode_args):
        """draw a graph using networkx and matplotlib"""
        if len(mode_args) > 0 and mode_args[0].lower() == "wheel":
            self.draw_wheel(*mode_args[1:])
            return
        draw(self._adj)


    def floyd(self):
        """Implements the Floyd-Warshall Algorithm, which returns an all-pairs
           shortest paths distance matrix"""
        distance = {v: {u: float('inf') for u in self.vertices} for v in self.vertices}
        for v in self.vertices:
            distance[v][v] = 0
            for u in self.adjacent_vertices(v):
                distance[v][u] = 1
        
        for k in self.vertices:
            for i in self.vertices:
                for j in self.vertices:
                    if distance[i][k] + distance[k][j] < distance[i][j]:
                        distance[i][j] = distance[i][k] + distance[k][j]
        return distance


    def diameter(self, d=None):
        """Use the floyd warshall algorithm to determine the diameter of the graph."""
        if d is None:
            d = self.floyd()
        return max(max(d[x].values()) for x in d)


    def print_distances(self, d=None):
        """Print a distance matrix for debugging"""
        if d is None:
            d = self.floyd()
        char_size = max(len(str(x)) for x in self.vertices)+1
        v_seq = sorted(self.vertices)
        print(' '.join([' '*char_size] + list(map(lambda x: f"{x:^{char_size}}", v_seq))))
        for v in v_seq:
            print(' '.join(
                [f"{v:^{char_size}}"] + list(map(lambda x: f"{x:^{char_size}}", [d[v][u] for u in v_seq]))
            ))


    def is_k_vertex_connected(self, k, return_cut=False):
        """A slow method for determinining if a graph is k vertex connected"""
        candidate_vertices = None
        for candidate_vertices in combinations(self.vertices, k):
            current_graph = Graph()
            current_graph.import_dict(self._adj)
            for v in candidate_vertices:
                current_graph.delete_vertex(v)
                if len(current_graph.components) > 1:
                    return (False, candidate_vertices) if return_cut else False
        return (True, candidate_vertices) if return_cut else True


    def get_vertex_connectivity(self, return_cut=False):
        """A very slow method to find the vertex connectivity of a graph"""
        k = 0
        while True:
            is_k, cut = self.is_k_vertex_connected(k, return_cut=True)
            if not is_k:
                return (k, cut) if return_cut else k
            k += 1


    def is_k_edge_connected(self, k, return_cut=False):
        """A slow method for determining if a graph is k edge connected"""
        candidate_edges = None
        for candidate_edges in combinations(self.edges, k):
            current_graph = Graph()
            current_graph.import_dict(self._adj)
            for e in candidate_edges:
                current_graph.delete_edge(e)
                if len(current_graph.components) > 1:
                    return (False, candidate_edges) if return_cut else False
        return (True, candidate_edges) if return_cut else True


    def get_edge_connectivity(self, return_cut=False):
        """A very slow method to find the edge connectivity of a graph"""
        k = 0
        while True:
            is_k, cut = self.is_k_edge_connected(k, return_cut=True)
            if not is_k:
                return (k, cut) if return_cut else k
            k += 1


    def augment(self, v, e1, e2):
        """Implements a graph operation which increases the order of a simplicial
           polyhedron by stretching a vertex into two, where for some planar embedding,
           the vertices between e1 and e2 (inclusive) are joined to one of the
           stretched vertices, and the vertices between e2 and e1 (inclusive)
           are joined to the other."""
        if len(set([v, e1, e2])) != 3:
            raise ValueError("Vertices must be distinct")
        wheel_vertices = self.adjacent_vertices(v)
        u = self.add_vertex()
        self.add_edge(u, e1)
        self.add_edge(u, v)
        self._adj[u].append(e2)
        self._adj[e2].insert(self._adj[e2].index(0), u)
        curr_vert_index = (wheel_vertices.index(e1)-1) % len(wheel_vertices)
        edge_indices_to_remove = []
        edge_indices_to_add = []
        while wheel_vertices[curr_vert_index] != e2:
            if curr_vert_index == 0:
                curr_vert_index = self.degree(v)
            edge_indices_to_remove.append(curr_vert_index)
            self._adj[wheel_vertices[curr_vert_index]][self._adj[wheel_vertices[curr_vert_index]].index(v)] = u
            edge_indices_to_add.append(wheel_vertices[curr_vert_index])

            curr_vert_index -= 1

        self._adj[u] += edge_indices_to_add[::-1]
        for index in sorted(edge_indices_to_remove, reverse=True):
            self._adj[v].pop(index)
        return u


    def stellate(self, vertices):
        """Implements a graph operation which adds a new vertex to the face of a
        planar graph, and joins this vertex to every vertex on the boundary of this face."""
        if not self.is_face(vertices, direction='cw'):
            if self.is_face(vertices, direction='ccw'):
                vertices = list(reversed(vertices))
            else:
                raise ValueError(f"{vertices} is not a face")
        new_vertex = self.add_vertex()
        self._adj[new_vertex] = list(vertices)
        for i, v in enumerate(vertices):
            self._adj[v].insert(self._adj[v].index(vertices[(i-1) % len(vertices)]), new_vertex)
        return new_vertex


    def draw_wheel(self, v, u=None):
        """draw a wheel graph using matplotlib such that the vertices are drawn
        around a central vertex v"""
        hull_vertices = self.adjacent_vertices(v)
        if u is not None:
            hull_vertices.remove(u)
            start_u_index = self.adjacent_vertices(u).index(v)
            insert_value = self.adjacent_vertices(u)[start_u_index + 1]
            start_split_index = hull_vertices.index(insert_value)
            vertices_to_add = []
            for i in range(2, self.degree(u)-1):
                vertices_to_add.append(self.adjacent_vertices(u)[i%self.degree(u) + start_u_index])
            i=0
            for i, vertex in enumerate(vertices_to_add):
                hull_vertices.insert(start_split_index+i+1, vertex)
            end_split_index = (start_split_index+i+2) % len(hull_vertices)
            middle = (end_split_index - start_split_index) / 2
            if len(vertices_to_add) == 0:
                middle -= 0.5
            if middle < 0:
                middle += len(hull_vertices) / 2
            middle += start_split_index

        size = 5
        positions = {}
        plt.figure(figsize=(10, 10))
        for i, vertex in enumerate(hull_vertices):
            x = np.sin(i * 2 * np.pi / len(hull_vertices)) * size
            y = np.cos(i * 2 * np.pi / len(hull_vertices)) * size
            positions[vertex] = (x,y)
        if u is None:
            positions[v] = (0, 0)
        else:
            positions[v] = (
                -np.sin(middle * 2 * np.pi / len(hull_vertices)) * size*0.3,
                -np.cos(middle * 2 * np.pi / len(hull_vertices)) * size*0.3
            )
            positions[u] = (
                np.sin(middle * 2 * np.pi / len(hull_vertices)) * size*0.3,
                np.cos(middle * 2 * np.pi / len(hull_vertices)) * size*0.3
            )

        for k, (x, y) in positions.items():
            plt.scatter(x ,y, s=500, c='b', edgecolors='k', zorder=10)
            plt.annotate(k, (x, y), size=15, color='w', ha='center', va='center', zorder=15)

        for v1, v2 in self.edges:
            if v1 in positions.keys() and v2 in positions.keys():
                plt.plot(
                    (positions[v1][0], positions[v2][0]),
                    (positions[v1][1], positions[v2][1]),
                    'k',
                    zorder=5
                )

        plt.xlim((-size*1.1, size*1.1))
        plt.ylim((-size*1.1, size*1.1))
        plt.show()


def compare_paths(p, q):
    """Return True if two paths are identical or the same path in different directions"""
    return len(p) == len(q) and (p == q or all([x == y for x, y in zip(p, reversed(q))]))


def generate_paths(graph, current_path=None):
    """Recursively generates a list of all paths in a graph starting from a vertex or path"""
    if current_path is None:
        current_path = [graph.vertices[0]]
    if type(current_path) != list:
        current_path = [current_path] # if a starting vertex is given, convert to singleton path
    neighbours = graph.adjacent_vertices(current_path[-1])
    yield current_path
    for n in neighbours:
        if n not in current_path:
            yield from generate_paths(graph, current_path + [n])


def longest_paths(graph, clean=False, max_path_len=None):
    """Returns a list containing each of the longest paths in a graph"""
    curr_len = 0
    curr_paths = []
    done_vertices = []
    for v in graph.vertices:
        for p in generate_paths(graph, v):
            if len(p) > curr_len:
                if max_path_len is not None and len(p) >= max_path_len:
                    return []
                curr_len = len(p)
                curr_paths = []
            if len(p) == curr_len and p[-1] not in done_vertices:
                curr_paths.append(p)
        done_vertices.append(v)

    if not clean:
        return curr_paths
    paths = []
    for p in curr_paths:
        if p not in paths and p[::-1] not in paths:
            paths.append(p)
    return paths

def generate_random_uniform_graph(n, p):
    """Generate a random uniform graph"""
    vertices = range(n)
    edges = []
    for u, v in combinations(vertices, 2):
        if random() <= p:
            edges.append({u, v})
    return Graph(vertices, edges)


def test_graph_file(num_vertices, display=True, trim=False):
    """load a pl file and test all stored graphs"""
    filename = f'../Databases/convexpolytopes_{num_vertices}.pl'
    #filename = f'/media/max/New Volume/pl/{num_vertices}.pl'
    graphs = read_graphs(filename)
    if trim:
        path_limit = (2 * num_vertices / 3) + 1
    else:
        path_limit = None
    graphs_done = 0
    for adj in graphs:
        G = Graph()
        G.import_dict(adj)
        paths = []
        paths = longest_paths(G, clean=True, max_path_len=path_limit)
        graphs_done += 1
        if graphs_done % 1000000 == 0:
            #exit()
            print(f'Checked {graphs_done}')
        if display and not paths:
            print(f"{adj=}")
            print(f"Found {len(paths)} longest paths of length"
                  +f"{len(paths[0]) - 1} ({len(paths[0])} vertices)")
            print(f"{paths=}")
            print()


def make_wheel(n):
    """make a wheel graph with n spokes"""
    if n < 3:
        raise ValueError("n must be at least 3")
    G = Graph()
    G.add_vertex(0)
    for i in range(1, n+1):
        G.add_edge(i, i%n+1)
        G.add_edge(0, i)
    return G
