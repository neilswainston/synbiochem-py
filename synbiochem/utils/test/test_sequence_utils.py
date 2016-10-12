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

    def test_get_codon_prob(self):
        '''Tests get_codon_prob method.'''
        cod_opt = seq_utils.CodonOptimiser('83333')
        self.assertAlmostEquals(0.46, cod_opt.get_codon_prob('CTG'), 2)

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

    def test_get_codon_usage_organisms(self):
        '''Tests get_codon_usage_organisms() method.'''
        organisms = seq_utils.get_codon_usage_organisms()
        expected = {'\'Arthromyces ramosus\'': '5451',
                    'thiosulfate-reducing anaerobe SRL 4198': '267367'}
        self.assertDictContainsSubset(expected, organisms)

    def test_get_uniprot_values(self):
        '''Tests get_uniprot_values method.'''
        self.assertEquals(seq_utils.get_uniprot_values(['P19367', 'P46882'],
                                                       ['organism-id'], 1),
                          {'P19367': {'Organism ID': '9606'},
                           'P46882': {'Organism ID': '5061'}})

    def test_translate(self):
        '''Tests translate method.'''
        results = seq_utils.translate('agcgtgcgatcc', min_prot_len=1)
        self.assertIn('ACD', [tokens[5] for tokens in results])

    def test_do_blast(self):
        '''Tests do_blast method.'''
        id_seq = {'test': seq_utils.get_random_dna(1024)}
        results = seq_utils.do_blast(id_seq, id_seq, evalue=10, word_size=4)

        alignments = []

        for result in results:
            for alignment in result.alignments:
                for hsp in alignment.hsps:
                    alignments.append(hsp)
                    print hsp

        self.assertGreater(len(alignments), 1)

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
