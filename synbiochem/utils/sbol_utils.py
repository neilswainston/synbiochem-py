'''
Created on 9 Mar 2016

@author: neilswainston
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
