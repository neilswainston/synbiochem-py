'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=too-many-public-methods
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


class Test(unittest.TestCase):
    '''Test class for sequence_utils.'''

    def test_get_uniprot_values(self):
        '''Tests get_uniprot_values method.'''
        self.assertEquals(seq_utils.get_uniprot_values(['P19367', 'P46882'],
                                                       ['organism-id'], 1),
                          {'P19367':
                           {'Organism ID': '9606', 'Entry': 'P19367'},
                           'P46882':
                           {'Organism ID': '5061', 'Entry': 'P46882'}})

    def test_find_orfs(self):
        '''Tests find_orfs method.'''
        results = seq_utils.find_orfs(
            'agcgtgcgat', min_prot_len=1)
        self.assertIn('ACD', [tokens[4] for tokens in results])

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
