'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=too-many-public-methods
from tempfile import NamedTemporaryFile
import unittest

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
                         sbol_doc2.num_sbol_objects - 2,
                         concat_doc.num_sbol_objects)

    def test_app_restrict_site_match(self):
        '''Tests apply_restriction_site method.'''
        docs = _get_apply_restrict_site_docs('(?<=GAGTC.{5}).*')
        self.assertEquals(len(docs), 3)

    def test_app_restrict_site_nomatch(self):
        '''Tests aplly_restriction_site method.'''
        docs = _get_apply_restrict_site_docs('(?<=JJJJJ.{5}).*')
        self.assertEquals(len(docs), 1)


def _get_apply_restrict_site_docs(restrict):
    '''Tests aplly_restriction_site method.'''
    parent_doc = Document()
    parent_doc.read('sbol.xml')
    docs = []

    for doc in sbol_utils.apply_restrict(parent_doc, restrict):
        out = NamedTemporaryFile()
        doc.write(out.name)
        read_doc = Document()
        read_doc.read(out.name)
        docs.append(read_doc)

    return docs


if __name__ == "__main__":
    unittest.main()
