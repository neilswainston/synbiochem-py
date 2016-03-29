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

        clone_doc = _round_trip(sbol_utils.clone(sbol_doc))

        self.assertEqual(sbol_doc.num_sbol_objects, clone_doc.num_sbol_objects)

    def test_concat(self):
        '''Tests concat method.'''
        sbol_doc2 = Document()
        sbol_doc2.read('sbol2.xml')
        sbol_doc3 = Document()
        sbol_doc3.read('sbol3.xml')

        concat_doc = _round_trip(sbol_utils.concat([sbol_doc2,
                                                    sbol_doc3]))

        self.assertEqual(sbol_doc2.num_sbol_objects +
                         sbol_doc3.num_sbol_objects - 2,
                         concat_doc.num_sbol_objects)

    def test_app_restrict_site_match(self):
        '''Tests apply_restriction_site method.'''
        _, docs = _get_apply_restrict_site_docs('(?<=gagtc.{5}).*')
        self.assertEquals(len(docs), 3)

    def test_app_restrict_site_nomatch(self):
        '''Tests aplly_restriction_site method.'''
        parent, docs = _get_apply_restrict_site_docs('(?<=JJJJJ.{5}).*')
        self.assertEquals(len(docs), 1)

        for par_annot, doc_annot in zip(parent.components[0].annotations,
                                        docs[0].components[0].annotations):
            self.assertEquals(par_annot.start, doc_annot.start)
            self.assertEquals(par_annot.end, doc_annot.end)


def _round_trip(doc):
    '''Saves a sbol Document to a temp file and re-reads it.'''
    out = NamedTemporaryFile()
    doc.write(out.name)
    doc = Document()
    doc.read(out.name)
    return doc


def _get_apply_restrict_site_docs(restrict):
    '''Tests aplly_restriction_site method.'''
    parent = Document()
    parent.read('sbol.xml')
    return parent, [_round_trip(doc)
                    for doc in sbol_utils.apply_restrict(parent, restrict)]


if __name__ == "__main__":
    unittest.main()
