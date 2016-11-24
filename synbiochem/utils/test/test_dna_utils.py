'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=no-member
# pylint: disable=too-many-public-methods
import os
import tempfile
import unittest

import synbiochem.utils.dna_utils as dna_utils


class Test(unittest.TestCase):
    '''Test class for dna_utils.'''

    def test(self):
        '''Tests round trip equality.'''
        directory = os.path.dirname(os.path.realpath(__file__))
        dna1 = dna_utils.read(os.path.join(directory, 'sbol.xml'))
        dna2 = _round_trip(dna1)
        self.assertEqual(dna1, dna2)

    def test_clone(self):
        '''Tests clone method.'''
        directory = os.path.dirname(os.path.realpath(__file__))
        dna1 = dna_utils.read(os.path.join(directory, 'sbol.xml'))
        dna2 = _round_trip(dna1.clone())
        self.assertEqual(dna1, dna2)

    def test_concat(self):
        '''Tests concat method.'''
        directory = os.path.dirname(os.path.realpath(__file__))
        dna2 = dna_utils.read(os.path.join(directory, 'sbol2.xml'))
        dna3 = dna_utils.read(os.path.join(directory, 'sbol3.xml'))
        concat_dna = _round_trip(dna_utils.concat([dna2, dna3]))

        self.assertFalse(concat_dna.features[0].forward)

        self.assertEqual(len(dna2.features) + len(dna3.features),
                         len(concat_dna.features))

    def test_app_restrict_site_linear(self):
        '''Tests apply_restriction_site method.'''
        _, dnas = _get_apply_restrict_site_dnas('(?<=gagtc.{5}).*', False)
        self.assertEquals([len(dna.seq) for dna in dnas], [25, 25, 831])

    def test_app_restrict_site_circular(self):
        '''Tests apply_restriction_site method.'''
        _, dnas = _get_apply_restrict_site_dnas('(?<=gagtc.{5}).*', True)
        self.assertEquals([len(dna.seq) for dna in dnas], [856, 25])

    def test_app_restrict_site_nomatch(self):
        '''Tests aplly_restriction_site method.'''
        parent, dnas = _get_apply_restrict_site_dnas('(?<=JJJJJ.{5}).*', False)
        self.assertEquals(len(dnas), 1)
        self.assertEquals(parent, dnas[0])


def _round_trip(dna):
    '''Writes / reads DNA object, via SBOL export / import.'''
    tmp = tempfile.NamedTemporaryFile()
    dna_utils.write(dna, tmp.name)
    return dna_utils.read(tmp.name)


def _get_apply_restrict_site_dnas(restr, circ):
    '''Tests apply_restriction_site method.'''
    directory = os.path.dirname(os.path.realpath(__file__))
    par = dna_utils.read(os.path.join(directory, 'restrict.xml'))
    return par, [_round_trip(dna)
                 for dna in dna_utils.apply_restricts(par, [restr], circ)]

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
