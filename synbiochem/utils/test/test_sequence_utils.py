'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import unittest

import synbiochem.utils.sequence_utils as seq_utils


class CodonOptimiserTest(unittest.TestCase):
    '''Test class for CodonOptimiser.'''

    def test_get_random_codon(self):
        '''Tests get_random_codon method.'''
        cod_opt = seq_utils.CodonOptimiser('83333')
        self.assertEquals('CTA', cod_opt.get_random_codon('L', ['CTG', 'TTA',
                                                                'CTT', 'TTG',
                                                                'CTC']))

    def test_get_random_codon_fail(self):
        '''Tests get_random_codon method.'''
        cod_opt = seq_utils.CodonOptimiser('83333')
        self.assertRaises(
            ValueError, cod_opt.get_random_codon, 'M', ['ATG'])


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
