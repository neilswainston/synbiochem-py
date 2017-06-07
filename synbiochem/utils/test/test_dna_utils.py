'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=no-member
# pylint: disable=too-many-public-methods
import json
import os
import unittest

from Bio import Restriction
from synbiochem.utils import dna_utils, sbol_utils

import synbiochem.utils.test.test_sbol_utils as test_sbol_utils


class Test(unittest.TestCase):
    '''Test class for dna_utils.'''

    def test_copy(self):
        '''Tests copy method.'''
        directory = os.path.dirname(os.path.realpath(__file__))
        dna1 = sbol_utils.read(os.path.join(directory, 'sbol.xml'))
        dna2 = test_sbol_utils.round_trip(dna1.copy())
        self.assertEqual(dna1, dna2)

    def test_concat(self):
        '''Tests concat method.'''
        directory = os.path.dirname(os.path.realpath(__file__))
        dna2 = sbol_utils.read(os.path.join(directory, 'sbol2.xml'))
        dna3 = sbol_utils.read(os.path.join(directory, 'sbol3.xml'))
        concat_dna = test_sbol_utils.round_trip(dna_utils.concat([dna2, dna3]))

        self.assertFalse(concat_dna['features'][0]['forward'])

        self.assertEqual(len(dna2['features']) + len(dna3['features']),
                         len(concat_dna['features']))

    def test_json(self):
        '''Tests json roundtrip.'''
        directory = os.path.dirname(os.path.realpath(__file__))
        dna1 = sbol_utils.read(os.path.join(directory, 'sbol.xml'))
        params = json.loads(json.dumps(dna1))
        dna2 = dna_utils.DNA(**params)
        self.assertEqual(dna1, dna2)

    def test_app_restrict_site_linear(self):
        '''Tests apply_restriction_site method.'''
        _, dnas = _get_apply_restrict_site_dnas(Restriction.MlyI, False)
        self.assertEquals([len(dna['seq']) for dna in dnas], [25, 831, 25])

    def test_app_restrict_site_circular(self):
        '''Tests apply_restriction_site method.'''
        _, dnas = _get_apply_restrict_site_dnas('MlyI', True)
        self.assertEquals([len(dna['seq']) for dna in dnas], [50, 831])

    def test_app_restrict_site_nomatch(self):
        '''Tests aplly_restriction_site method.'''
        parent, dnas = _get_apply_restrict_site_dnas(Restriction.HgaI, False)
        self.assertEquals(len(dnas), 1)
        self.assertEquals(parent, dnas[0])


def _get_apply_restrict_site_dnas(restr, circ):
    '''Tests apply_restriction_site method.'''
    directory = os.path.dirname(os.path.realpath(__file__))
    par = sbol_utils.read(os.path.join(directory, 'restrict.xml'))
    return par, [test_sbol_utils.round_trip(dna)
                 for dna in dna_utils.apply_restricts(par, [restr], circ)]


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
