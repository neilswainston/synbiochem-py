'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import random
import unittest

from synbiochem.utils import mut_utils, seq_utils


class MutationTest(unittest.TestCase):
    '''Test class for Mutation.'''

    def test_init_err_wt_seq_aa(self):
        '''Tests __init__.'''
        self.assertRaises(AssertionError, mut_utils.Mutation, 'J', 12, 'A')

    def test_init_err_mut_seq_aa(self):
        '''Tests __init__.'''
        self.assertRaises(AssertionError, mut_utils.Mutation, 'A', 12, 'J')

    def test_init_err_wt_seq_dna(self):
        '''Tests __init__.'''
        self.assertRaises(AssertionError, mut_utils.Mutation, 'D', 12, 'A',
                          'dna')

    def test_init_err_mut_seq_dna(self):
        '''Tests __init__.'''
        self.assertRaises(AssertionError, mut_utils.Mutation, 'A', 12, 'D',
                          'dna')

    def test_init_err_pos_zero(self):
        '''Tests __init__.'''
        self.assertRaises(AssertionError, mut_utils.Mutation, 'A', 0, 'A',
                          'dna')

    def test_init_err_pos_neg(self):
        '''Tests __init__.'''
        self.assertRaises(AssertionError, mut_utils.Mutation, 'A', -5, 'A',
                          'dna')

    def test_init_err_pos_non_int(self):
        '''Tests __init__.'''
        self.assertRaises(AssertionError, mut_utils.Mutation, 'A', 5.6, 'A',
                          'dna')

    def test_repr(self):
        '''Tests __repr__.'''
        mut = mut_utils.Mutation('A', 12, '-')
        self.assertEqual(mut.__repr__(), 'A12-')

    def test_str(self):
        '''Tests __str__.'''
        mut = mut_utils.Mutation('A', 12, 'V')
        self.assertEqual(mut.__str__(), 'A12V')

    def test_eq(self):
        '''Tests __eq__.'''
        mut_a = mut_utils.Mutation('A', 12, 'V')
        mut_b = mut_utils.Mutation('A', 12, 'V')
        mut_c = mut_utils.Mutation('A', 13, 'V')
        mut_d = mut_utils.Mutation('A', 12, 'T')

        self.assertEqual(mut_a, mut_b)
        self.assertNotEqual(mut_a, mut_c)
        self.assertNotEqual(mut_a, mut_d)

    def test_cmp(self):
        '''Tests __cmp__.'''
        mut_a = mut_utils.Mutation('A', 12, 'V')
        mut_b = mut_utils.Mutation('A', 5, 'V')
        mut_c = mut_utils.Mutation('A', 12, 'T')
        lst = [mut_a, mut_b, mut_c]
        sort_lst = [mut_b, mut_c, mut_a]

        self.assertNotEqual(lst, sorted(lst))
        self.assertEqual(sort_lst, sorted(lst))


class Test(unittest.TestCase):
    '''Test class for sequence_utils.'''

    def test_parse_mut_str(self):
        '''Tests parse_mut_str method.'''
        mut_a = mut_utils.Mutation('A', 12, 'V')
        mut_b = mut_utils.Mutation('P', 5, '-')
        mut_lst = [mut_a, mut_b]
        mut_str = ' '.join([str(mut) for mut in mut_lst])

        self.assertEqual(mut_str, 'A12V P5-')
        self.assertEqual(mut_utils.parse_mut_str(mut_str), mut_lst)

    def test_get_mutations(self):
        '''Tests get_mutations method.'''
        wt_seq = 'QWERTY'
        mt_seq = 'QTER-Y'
        mut_a = mut_utils.Mutation('W', 2, 'T')
        mut_b = mut_utils.Mutation('T', 5, '-')
        mut_lst = [mut_a, mut_b]
        self.assertEqual(mut_utils.get_mutations(wt_seq, mt_seq), mut_lst)

    def test_apply_mutations(self):
        '''Tests apply_mutations method.'''
        aa_len = 1024
        wt_seq = [random.choice(list(seq_utils.AA_CODES.values()))
                  for _ in range(aa_len)]
        mt_seq = list(wt_seq)

        mutations = []

        for _ in range(16):
            pos = random.randint(0, aa_len)
            wt_res = wt_seq[pos]
            aa_codes = list(seq_utils.AA_CODES.values())
            aa_codes.remove('*')
            mt_res = random.choice(aa_codes)
            mt_seq[pos] = mt_res
            mutations.append(mut_utils.Mutation(wt_res, pos + 1, mt_res))

        self.assertEqual(mut_utils.apply_mutations(wt_seq, mutations),
                         ''.join(mt_seq))


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
