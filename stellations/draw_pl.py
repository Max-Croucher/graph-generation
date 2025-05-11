""" Render all graphs from a generated .pl file.
Arguments: <filename> [--save output_dir]
where:
    filename is the .pl file to draw
    if --save is set, instead of displaying the images, an output directory must
    be provided, in which each image will be saved.

Max Croucher, March 2024
"""


from sys import argv
import os
import networkx as nx
import matplotlib.pyplot as plt
from read_pl import read_graphs_from_file


WIDTH = 20
HEIGHT = 20


def draw(n, fname=None, edge_colourings=None, title=None, text=None):
    """Draw a planar drawing of a graph"""
    plt.figure(figsize=(WIDTH, HEIGHT))
    g = nx.PlanarEmbedding()
    g.set_data(n)

    if edge_colourings is not None:
        colourings = []
        for i, j in g.edges():
            if (i, j) in edge_colourings:
                colourings.append(edge_colourings[i, j])
            elif (j, i) in edge_colourings:
                colourings.append(edge_colourings[j, i])
            else:
                colourings.append('k')

    pos = nx.planar_layout(g)
    nx.draw_networkx(
        g, pos,
        with_labels=True,
        edge_color=(colourings if edge_colourings is not None else 'k'),
        node_size=100,
        linewidths=2
    )
    if title is not None:
        plt.title(title)
    if text is not None:
        plt.text(0, -0.55, text, horizontalalignment='center', wrap=True )
    if fname is None:
        plt.show()
    else:
        plt.savefig(os.path.join(fname))
    plt.close()


def main(filename, outdir=None):
    """main for cli"""
    if outdir is not None:
        try:
            os.mkdir(outdir)
        except FileExistsError:
            pass

    for i, graph_dict in enumerate(read_graphs_from_file(filename)):
        if outdir is not None:
            draw(graph_dict, os.path.join(outdir, f"{i}.png"))
        else:
            draw(graph_dict)


if __name__ == '__main__':
    if len(argv) == 2:
        main(argv[1])
    if len(argv) == 4 and argv[2].lower() in ['-s', '--save']:
        main(argv[1], argv[3])
    else:
        print(f"Usage: python {argv[0]} <pl filename> [--save ouptut_directory]")
