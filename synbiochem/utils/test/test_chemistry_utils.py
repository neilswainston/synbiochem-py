'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import unittest

import synbiochem.utils.chemistry_utils as chm_util


class Test(unittest.TestCase):
    '''Test class for chemical_utils.'''

    def test_get_molecular_mass(self):
        '''Tests get_molecular_mass method.'''
        self.assertAlmostEqual(chm_util.get_molecular_mass('H2O'),
                               18.010564684)

    def test_get_elem_comp(self):
        '''Tests get_elem_comp method.'''
        self.assertEqual(chm_util.get_elem_comp('Fe4S'), {'Fe': 4, 'S': 1})

    def test_balance_unbalanced(self):
        '''Tests get_elem_comp method.'''
        unbalanced = [('CO2', 0, -1.0), ('C5H7O4', -1, -1.0),
                      ('C3H3O3', -1, 1.0)]
        is_balanced, was_balanced, _ = chm_util.balance(unbalanced)

        self.assertTrue(is_balanced)
        self.assertFalse(was_balanced)

    def test_balance_balanced(self):
        '''Tests get_elem_comp method.'''
        balanced = [('C5H7O4', -1, -1.0), ('H', 1, 1.0),
                    ('C3H3O3', -1, 2.0), ('CO2', 0, -1.0)]
        is_balanced, was_balanced, _ = chm_util.balance(balanced)

        self.assertTrue(is_balanced)
        self.assertTrue(was_balanced)


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
