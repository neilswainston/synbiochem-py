'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=invalid-name


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
