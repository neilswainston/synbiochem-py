'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''


class Vertex():
    '''Class to represent a vertex.'''

    def __init__(self, name, attributes):
        self.__attributes = attributes
        self.__attributes['name'] = name
        self.__in = []
        self.__out = []

    def predecessors(self):
        '''Get predecessors.'''
        return self.__in

    def attributes(self):
        '''Get attributes.'''
        return self.__attributes

    def add_out(self, vertex_to, attributes):
        '''Add edge out.'''
        self.__out.append((vertex_to, attributes))

    def add_in(self, vertex_from, attributes):
        '''Add edge in.'''
        self.__in.append((vertex_from, attributes))

    def indegree(self):
        '''Get indegree.'''
        return len(self.__in)

    def is_root(self):
        '''Is the Vertex a root?'''
        return not self.__out

    def __repr__(self):
        return self.__attributes['name']


class Graph():
    '''Class to represent a graph.'''

    def __init__(self):
        self.__vertices = {}

    def find_vertex(self, name):
        '''Find vertex.'''
        if name in self.__vertices:
            return self.__vertices[name]

        raise ValueError(name)

    def add_vertex(self, name, attributes):
        '''Add vertex.'''
        self.__vertices[name] = Vertex(name, attributes)

    def add_edge(self, vertex_from, vertex_to, attributes):
        '''Add edge.'''
        vertex_from.add_out(vertex_to, attributes)
        vertex_to.add_in(vertex_from, attributes)

    def get_roots(self):
        '''Get roots.'''
        return [vtx for vtx in self.__vertices.values() if vtx.is_root()]


def add_vertex(graph, name, attributes=None):
    '''Add vertex.'''
    try:
        return graph.find_vertex(name)
    except ValueError:
        if not attributes:
            attributes = {}

        graph.add_vertex(name, attributes)
        return graph.find_vertex(name)


def add_edge(graph, vertex_from, vertex_to, attributes=None):
    '''Add edge.'''
    if not attributes:
        attributes = {}

    graph.add_edge(vertex_from, vertex_to, attributes)


def get_roots(graph):
    '''Get roots.'''
    return graph.get_roots()
