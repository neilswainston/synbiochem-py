'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=too-many-public-methods
import random
import unittest

from synbiochem.utils import seq_utils


class CodonOptimiserTest(unittest.TestCase):
    '''Test class for CodonOptimiser.'''

    def test_get_codon_prob(self):
        '''Tests get_codon_prob method.'''
        cod_opt = seq_utils.CodonOptimiser('83333')
        self.assertAlmostEquals(0.46, cod_opt.get_codon_prob('CTG'), 2)

    def test_get_codon_optim_seq(self):
        '''Tests get_codon_optim_seq method.'''
        cod_opt = seq_utils.CodonOptimiser('83333')
        aa_codes = seq_utils.AA_CODES
        aa_codes.pop('End')

        aa_seq = ''.join([random.choice(aa_codes.values())
                          for _ in range(random.randint(100, 2500))])

        max_repeat_nuc = 5
        restr_enzyms = ['BsaI']
        dna_seq = cod_opt.get_codon_optim_seq(aa_seq,
                                              max_repeat_nuc=max_repeat_nuc,
                                              restr_enzyms=restr_enzyms)

        self.assertFalse(seq_utils.is_invalid(dna_seq, max_repeat_nuc,
                                              restr_enzyms))

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
        '''Tests get_codon_usage_organisms method.'''
        organisms = seq_utils.get_codon_usage_organisms()
        expected = {'\'Arthromyces ramosus\'': '5451',
                    'thiosulfate-reducing anaerobe SRL 4198': '267367'}
        self.assertDictContainsSubset(expected, organisms)

    def test_find_invalid(self):
        '''Tests find_invalid method.'''
        seq = 'ggtctaaaaatttttttaaaaaccagagtttttt'

        self.assertEquals(seq_utils.find_invalid(seq, 5, ['BsaI']),
                          [10, 11, 28])

    def test_is_invalid(self):
        '''Tests is_invalid method.'''
        seq_inv = 'ggtctaaaaatttttttaaaaaccagagtttttt'
        self.assertTrue(seq_utils.is_invalid(seq_inv, 5, ['BsaI']))

        seq_val = 'ggtctaaaa'
        self.assertFalse(seq_utils.is_invalid(seq_val, 5, ['BsaI']))

    def test_get_random_dna(self):
        '''Tests get_random_dna method.'''
        lngth = random.randint(1000, 10000)
        self.assertEqual(lngth,
                         len(seq_utils.get_random_dna(lngth, 4, ['BsaI'])))

    def test_get_seq_by_melt_temp(self):
        '''Tests get_seq_by_melt_temp method.'''
        seq, _ = seq_utils.get_seq_by_melt_temp('agcgtgcgaagcgtgcgatcctcc', 70)
        self.assertEqual(seq, 'agcgtgcgaagcgtgcgatc')

    def test_get_rand_seq_by_melt_temp(self):
        '''Tests get_rand_seq_by_melt_temp method.'''
        target_temp = random.randint(50, 100)
        _, temp = seq_utils.get_rand_seq_by_melt_temp(target_temp, 4, ['BsaI'])
        self.assertTrue(abs(target_temp - temp) / target_temp < 0.025)

    def test_get_uniprot_values(self):
        '''Tests get_uniprot_values method.'''
        result = seq_utils.get_uniprot_values(['P19367', 'P42212'],
                                              ['organism-id',
                                               'protein names'], 1)

        expected = {'P19367': {'Organism ID': '9606',
                               'Protein names': ['Hexokinase-1',
                                                 'EC 2.7.1.1',
                                                 'Brain form hexokinase',
                                                 'Hexokinase type I',
                                                 'HK I']},
                    'P42212': {'Organism ID': '6100',
                               'Protein names': ['Green fluorescent protein']}}

        self.assertEquals(result, expected)

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
