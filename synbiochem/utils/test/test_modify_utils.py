'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=too-many-public-methods
import unittest

from synbiochem.utils.modify_utils import CodonSelector


class Test(unittest.TestCase):
    '''Test class for modify_utils.'''

    def test_optimise_codons(self):
        '''Tests isclose method.'''
        cod_sel = CodonSelector()
        codons = cod_sel.optimise_codons({'aminoAcids': 'FLIMV',
                                          'organism': {'id': '37762'}})

        self.assertEqual(len(codons), 12)
        self.assertEqual(codons[0][0], 'DTK')


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
