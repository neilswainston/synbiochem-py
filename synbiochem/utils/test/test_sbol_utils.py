'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import unittest

from django.core.files.temp import NamedTemporaryFile
from sbol.sbol import Document

import synbiochem.utils.sbol_utils as sbol_utils


class Test(unittest.TestCase):
    '''Test class for chemical_utils.'''

    def test_clone(self):
        '''Tests clone method.'''
        sbol_doc = Document()
        sbol_doc.read('sbol.xml')

        out = NamedTemporaryFile()
        clone_doc = sbol_utils.clone(sbol_doc)
        clone_doc.write(out.name)

        clone_doc = Document()
        clone_doc.read(out.name)

        self.assertEqual(sbol_doc.num_sbol_objects, clone_doc.num_sbol_objects)

    def test_concatenate(self):
        '''Tests concatenate method.'''
        sbol_doc1 = Document()
        sbol_doc1.read('sbol.xml')
        sbol_doc2 = Document()
        sbol_doc2.read('sbol2.xml')

        out = NamedTemporaryFile()
        concat_doc = sbol_utils.concatenate([sbol_doc1, sbol_doc2])
        concat_doc.write(out.name)

        concat_doc = Document()
        concat_doc.read(out.name)

        self.assertEqual(sbol_doc1.num_sbol_objects +
                         sbol_doc2.num_sbol_objects,
                         concat_doc.num_sbol_objects)

if __name__ == "__main__":
    unittest.main()
