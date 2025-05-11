""" Read a generated .pl file, to count, print or store results.
Arguments: <filename> [mode]
where:
    filename is the .pl file to read
    mode is either "count", "print" or "store". If not provided, the argument
        defaults to "store". If mode is "count", the program will count the
        number of graphs in the file. If mode is "print", each graphs will be
        printed to console as a dictionary. If mode is "store", each graph will
        be saved to filename.txt as a textual representation of each graph as a
        dict.

Max Croucher, March 2024
"""


from sys import argv
from pathlib import Path


class PLWriter():
    """A file handling class to implement writing to a .pl file in the planar_code
       file format. Designed to be used with 'with'"""
    HEADER = b">>planar_code<<"
    def __init__(self, filename, mode='wb', ignore_header=False):
        self.filename = filename
        self.mode = mode
        self.ignore_header=ignore_header
        self.obj = None


    def __enter__(self):
        self.obj = open(self.filename, self.mode)
        if not self.ignore_header:
            self.obj.write(self.HEADER)
        return self.obj


    def __exit__(self, *exceptions):
        self.obj.close()


def load_graph_file(filename):
    """Load a file in binary mode"""
    return open(filename, 'rb')


def read_graphs_from_file(filename):
    """Read a file and return each graph as a dictionary"""
    file = load_graph_file(filename)
    if file.read(15) != b'>>planar_code<<':
        yield from []
    else:
        vertex_count = 1
        while True:
            b = file.read(1)
            if not b:
                break # eof
            num_vertices = int.from_bytes(b, 'big')
            vertex_count = 1
            adj_list = dict()
            while vertex_count <= num_vertices:
                b = file.read(1)
                adj_list[vertex_count] = []
                while b and b != b'\x00':
                    adj_list[vertex_count].append(int.from_bytes(b, 'big'))
                    b = file.read(1)
                vertex_count += 1
            yield adj_list
        file.close()


def read_graphs_from_dir(dirname):
    """Recursively read graphs in a directory by attempting to read graphs from
    every file in the directory"""
    for file in dirname.iterdir():
        yield from read_graphs(file)


def read_graphs(filename):
    """Read graphs from a file or from contents of a directory"""
    if filename.is_dir():
        yield from read_graphs_from_dir(filename)
    else:
        yield from read_graphs_from_file(filename)


def print_graphs(filename):
    """Print all generated graphs to stdout"""
    g = read_graphs(filename)
    for graph in g:
        for v, a in graph.items():
            print(f"{v}: {a}")
        print()


def count_graphs(filename):
    """Count all graphs in a file or directory"""
    g = read_graphs(filename)
    count = 0
    for _ in g:
        count += 1
    return count


def store_graphs(filename):
    """Read all graphs from a file or directory and store results in a text file"""
    g = read_graphs(filename)
    with open(f'{filename}.txt', 'w', encoding='utf-8') as f:
        for graph in g:
            for v, a in graph.items():
                f.write(f"{v}: {a}\n")
            f.write('\n')


def write_graph(file_obj, adj):
    """Write an adjacency list to a graph file"""
    bytes_to_write = [len(adj)]
    if 0 in adj.keys():
        raise ValueError("A vertex name of 0 is forbidden for writing")
    for i in range(1, len(adj)+1):
        if i not in adj.keys():
            raise ValueError("Vertex names must be consecutive integers beginning with 1")
        bytes_to_write += adj[i] + [0]
    file_obj.write(bytearray(bytes_to_write))


def split_pl(input_name, output_name, count):
    """Partition a planar grapf file into count disjoint parts. Useful for parallelisation"""
    output_name.mkdir(parents=True, exist_ok=True)
    pl_files = [PLWriter(output_name.joinpath(f"{i}.pl")) for i in range(count)]
    [pl.__enter__() for pl in pl_files]
    for i, g in enumerate(read_graphs(input_name)):
        write_graph(pl_files[i%count].obj, g)
    [pl.__exit__() for pl in pl_files]


def main(*arguments):
    """main for cli"""
    if len(arguments) == 3 and arguments[2] == "count":
        print(count_graphs(Path(arguments[1])))
    elif len(arguments) == 3 and arguments[2] == "print":
        print_graphs(Path(arguments[1]))
    elif len(arguments) == 3 and arguments[2] == "store":
        store_graphs(Path(arguments[1]))
    elif len(arguments) == 4 and arguments[2] == "split":
        split_pl(Path(arguments[1]), Path(arguments[1]+".split"), int(arguments[3]))
    else:
        print(f"Usage: python3 {arguments[0]} <pl filename> [count | print | store | split <n>")


if __name__ == "__main__":
    main(*argv)
