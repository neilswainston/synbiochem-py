'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  pablocarbonell / neilswainston / alanwilliams
'''
from xml.etree import ElementTree
from xml.etree.ElementTree import Element, SubElement, tostring

_PATH_TO_SEQ = './sbol:DnaComponent/sbol:dnaSequence/sbol:DnaSequence' + \
    '/sbol:nucleotides'


class SBOLDocument(object):
    '''Class representing an SBOL Document.'''

    def __init__(self, filename=None):
        self.__ns = {'rdf': 'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
                     'sbol': 'http://sbols.org/v1#'}

        for prefix, uri in self.__ns.iteritems():
            ElementTree.register_namespace(prefix, uri)

        if filename is None:
            self.__root = Element(self.__get_ns('rdf') + 'RDF')
            self.__init_sequence()
        else:
            with open(filename) as fle:
                self.__root = ElementTree.parse(fle).getroot()

    def get_sequence(self):
        '''Gets the DNA sequence.'''
        return self.__get_seq_element().text.upper()

    def set_sequence(self, seq):
        '''Sets the DNA sequence.'''
        return self.__get_seq_element().set_text(seq)

    def __init_sequence(self, seq=''):
        '''Initialises a SBOLDocument.'''
        sbolns = self.__get_ns('sbol')
        child_elem = SubElement(self.__root, sbolns + 'DnaComponent')
        child_elem = SubElement(child_elem, sbolns + 'dnaSequence')
        child_elem = SubElement(child_elem, sbolns + 'DnaSequence')
        child_elem = SubElement(child_elem, sbolns + 'nucleotides')
        child_elem.text = seq

    def __get_seq_element(self):
        '''Gets the DNA sequence element.'''
        elements = self.__root.findall(_PATH_TO_SEQ, self.__ns)
        return elements[0]

    def __get_ns(self, prefix):
        '''Returns namespace string.'''
        return '{' + self.__ns[prefix] + '}'

    def __repr__(self):
        return tostring(self.__root)

    def __add__(self, other):
        return '%s plus %s' % (self, other)
