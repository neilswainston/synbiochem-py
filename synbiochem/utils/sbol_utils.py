'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  pablocarbonell / neilswainston / alanwilliams
'''
from xml.etree import ElementTree
from xml.etree.ElementTree import Element, SubElement, tostring
import uuid

_PATH_TO_DNA_COMP = './sbol:DnaComponent'

_PATH_TO_NAME = './sbol:DnaComponent/sbol:name'

_PATH_TO_SEQ = './sbol:DnaComponent/sbol:dnaSequence/sbol:DnaSequence' + \
    '/sbol:nucleotides'


class SBOLDocument(object):
    '''Class representing an SBOL Document.'''

    def __init__(self, uri_prefix=None, filename=None):
        self.__ns = {'rdf': 'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
                     'sbol': 'http://sbols.org/v1#'}

        for prefix, uri in self.__ns.iteritems():
            ElementTree.register_namespace(prefix, uri)

        if filename is None:
            self.__root = Element(self.__get_ns('rdf') + 'RDF')
            self.__uri_prefix = uri_prefix
            self.__init_sequence(str(uuid.uuid4()))
        else:
            with open(filename) as fle:
                self.__root = ElementTree.parse(fle).getroot()
                self.__uri_prefix = self.__get_uri_prefix()

    def get_sequence(self):
        '''Gets the DNA sequence.'''
        elem = self.__get_element(_PATH_TO_SEQ)
        return elem.text.upper()

    def get_name(self):
        '''Gets the name.'''
        elem = self.__get_element(_PATH_TO_NAME)
        return elem.text

    def set_sequence(self, seq):
        '''Sets the DNA sequence.'''
        elem = self.__get_element(_PATH_TO_SEQ)
        elem.text = seq

    def set_name(self, seq):
        '''Sets the DNA sequence.'''
        elem = self.__get_element(_PATH_TO_NAME)
        return elem.set_text(seq)

    def __init_sequence(self, doc_id, seq=''):
        '''Initialises a SBOLDocument.'''
        rdf_ns = self.__get_ns('rdf')
        sbol_ns = self.__get_ns('sbol')
        dna_elem = SubElement(self.__root,
                              sbol_ns + 'DnaComponent',
                              {rdf_ns + 'about': self.__uri_prefix + doc_id})
        child_elem = SubElement(dna_elem, sbol_ns + 'displayId')
        child_elem.text = doc_id
        child_elem = SubElement(dna_elem, sbol_ns + 'dnaSequence')
        child_elem = SubElement(child_elem,
                                sbol_ns + 'DnaSequence',
                                {rdf_ns + 'about': self.__uri_prefix +
                                 str(uuid.uuid4())})
        child_elem = SubElement(child_elem, sbol_ns + 'nucleotides')
        child_elem.text = seq

    def __get_element(self, path, element=None):
        '''Gets the element from the path.'''
        if element is None:
            element = self.__root

        elements = element.findall(path, self.__ns)

        if len(elements) > 0:
            return elements[0]
        else:
            return self.__create_element(path)

    def __get_ns(self, prefix):
        '''Returns namespace string.'''
        return '{' + self.__ns[prefix] + '}'

    def __get_uri_prefix(self):
        '''Returns uri prefix.'''
        elem = self.__get_element(_PATH_TO_DNA_COMP)
        uri = elem.get(self.__get_ns('rdf') + 'about')
        return uri[:uri.rfind('#')]

    def __create_element(self, path, element=None):
        tokens = path.split('/')
        element = self.__get_element(tokens[0], element)

    def __repr__(self):
        return tostring(self.__root)

    def __add__(self, other):
        concat = SBOLDocument(self.__uri_prefix)
        concat.set_sequence(self.get_sequence() + other.get_sequence())
        concat.set_name(self.get_name() + ' + ' + other.get_name())
        return concat
