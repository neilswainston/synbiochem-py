'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import unittest

from synbiochem.utils import mut_utils


class MutationTest(unittest.TestCase):
    '''Test class for CodonOptimiser.'''

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
        mut = mut_utils.Mutation('A', 12, 'V')
        self.assertEquals(mut.__repr__(), 'A12V')

    def test_str(self):
        '''Tests __str__.'''
        mut = mut_utils.Mutation('A', 12, 'V')
        self.assertEquals(mut.__str__(), 'A12V')

    def test_eq(self):
        '''Tests __eq__.'''
        mut_a = mut_utils.Mutation('A', 12, 'V')
        mut_b = mut_utils.Mutation('A', 12, 'V')
        mut_c = mut_utils.Mutation('A', 13, 'V')
        mut_d = mut_utils.Mutation('A', 12, 'T')

        self.assertEquals(mut_a, mut_b)
        self.assertNotEquals(mut_a, mut_c)
        self.assertNotEquals(mut_a, mut_d)

    def test_cmp(self):
        '''Tests __cmp__.'''
        mut_a = mut_utils.Mutation('A', 12, 'V')
        mut_b = mut_utils.Mutation('A', 5, 'V')
        mut_c = mut_utils.Mutation('A', 12, 'T')
        lst = [mut_a, mut_b, mut_c]
        sort_lst = [mut_b, mut_c, mut_a]

        self.assertNotEquals(lst, sorted(lst))
        self.assertEquals(sort_lst, sorted(lst))


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
