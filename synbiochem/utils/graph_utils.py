'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=invalid-name
from igraph import Graph
import matplotlib.pyplot as plt


def add_vertex(graph, name, kwds=None):
    '''Add vertex.'''
    try:
        return graph.vs.find(name)
    except ValueError:
        if not kwds:
            kwds = {}

        graph.add_vertex(name, **kwds)
        return graph.vs.find(name)


def add_edge(graph, vertex_from, vertex_to, kwds=None):
    '''Add edge.'''
    if not kwds:
        kwds = {}

    graph.add_edge(vertex_from.index, vertex_to.index, **kwds)


def get_roots(graph):
    '''Get roots.'''
    return [vs for vs, outdg in zip(graph.vs, graph.outdegree())
            if not outdg]


def plot_graph(graph, outfile=None, layout_name='kk'):
    '''Plot labelled graph.'''
    positions = _get_positions(graph, layout_name)

    xs, ys = zip(*positions.values())
    ax = plt.gca()
    ax.set_xlim(min(xs), max(xs))
    ax.set_ylim(min(ys), max(ys))

    # Plot nodes:
    labels = [vertex['name'] for vertex in graph.vs]

    colours = ['lightgreen' if vertex['is_reagent'] else 'lightblue'
               for vertex in graph.vs]

    patches = _plot_vertices(positions, labels, colours)

    # Plot edges:
    _plot_edges(positions, [e.tuple for e in graph.es], patches)

    plt.axis('off')

    if outfile:
        plt.savefig(outfile)
    else:
        plt.show()


def _get_positions(graph, layout_name):
    '''Gets positions.'''
    if layout_name == 'tree':
        roots = [v.index for v in get_roots(graph)]

        if not roots:
            roots = [0]

        layout = graph.layout(layout_name, mode='in', root=roots)
    else:
        layout = graph.layout(layout_name)

    return {k: layout[k] for k in range(len(layout.coords))}


def _plot_vertices(positions, labels, colours):
    '''Plot vertices.'''
    patches = []

    xs, ys = zip(*positions.values())

    for x, y, label, colour in zip(xs, ys, labels, colours):
        bbox_props = dict(boxstyle='circle,pad=0.5',
                          fc=colour,
                          ec='grey',
                          lw=1)

        patches.append(plt.text(x, y, label, ha='center', va='center',
                                size=10,
                                zorder=-1,
                                bbox=bbox_props))
    return patches


def _plot_edges(positions, edges, patches):
    '''Plot edges.'''
    for edge in edges:
        x = [positions[edge[1]][0], positions[edge[0]][0]]
        y = [positions[edge[1]][1], positions[edge[0]][1]]

        plt.annotate('',
                     xy=(x[0], y[0]), xycoords='data',
                     xytext=(x[1], y[1]), textcoords='data',
                     zorder=1,
                     arrowprops=dict(arrowstyle='->',
                                     color='0.5',
                                     patchA=patches[edge[0]],
                                     patchB=patches[edge[1]],
                                     ),
                     )


def _get_graph(vertices=16, children=2):
    '''Gets graph.'''
    labels = map(str, range(vertices))
    graph = Graph.Tree(vertices, children)
    return labels, graph
