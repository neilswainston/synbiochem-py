'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
from xml.etree import ElementTree


class SbolDocument(object):
    '''Simple, lightweight SBOL v1.1 reader / writer?'''

    def __init__(self, filename):
        with open(filename) as fle:
            self.__root = ElementTree.parse(fle).getroot()

        self.__ns = {'rdf': 'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
                     'sbol': 'http://sbols.org/v1#'}

    def get_dna_sequence(self):
        '''Gets the DNA sequence.'''
        path = './sbol:DnaComponent/sbol:dnaSequence/sbol:DnaSequence' + \
            '/sbol:nucleotides'
        elements = self.__root.findall(path, self.__ns)
        return elements[0].text.upper()

    def __add__(self, other):
        return '%s plus %s' % (self, other)
