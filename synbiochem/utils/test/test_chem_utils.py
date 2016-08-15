'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=too-many-public-methods
import unittest

import synbiochem.utils.chem_utils as chm_util


class Test(unittest.TestCase):
    '''Test class for chemical_utils.'''

    def test_get_molecular_mass(self):
        '''Tests get_molecular_mass method.'''
        self.assertAlmostEqual(chm_util.get_molecular_mass('H2O'),
                               18.010564684)

    def test_get_elem_comp(self):
        '''Tests get_elem_comp method.'''
        self.assertEqual(chm_util.get_elem_comp('Fe4S'), {'Fe': 4, 'S': 1})

    def test_parse_equation(self):
        '''Tests parse_equation method.'''
        eqn = '5.6 Fe4S + -3.2 water = 17.8 SiO2'
        self.assertEqual(chm_util.parse_equation(eqn),
                         [('Fe4S', -5.6), ('water', 3.2), ('SiO2', 17.8)])

    def test_parse_equation_error(self):
        '''Tests parse_equation method (with error).'''
        eqn = '5.6 Fe4S + -3.2 water = n+m SiO2'
        self.assertRaises(ValueError, chm_util.parse_equation, eqn)

    def test_balance_unbalanced(self):
        '''Tests get_elem_comp method for unbalanced reaction.'''
        unbalanced = [('CO2', 0, -1.0, 'CO2'),
                      ('C5H7O4', -1, -1.0, 'C5H7O4'),
                      ('C3H3O3', -1, 1.0, 'C3H3O3')]
        is_balanced, was_balanced, balanced_def = chm_util.balance(unbalanced)

        self.assertTrue(is_balanced)
        self.assertFalse(was_balanced)
        self.assertEqual(sorted([('C3H3O3', -1, 2.0, 'C3H3O3'),
                                 ('C5H7O4', -1, -1.0, 'C5H7O4'),
                                 ('CO2', 0, -1.0, 'CO2'),
                                 ('H', 1, 1.0, 'CHEBI:24636')]),
                         sorted(balanced_def))

    def test_balance_balanced(self):
        '''Tests get_elem_comp method for balanced reaction.'''
        balanced = [('C5H7O4', -1, -1.0, 'C5H7O4'),
                    ('H', 1, 1.0, 'H'),
                    ('C3H3O3', -1, 2.0, 'C3H3O3'),
                    ('CO2', 0, -1.0, 'CO2')]
        is_balanced, was_balanced, balanced_def = chm_util.balance(balanced)

        self.assertTrue(is_balanced)
        self.assertTrue(was_balanced)
        self.assertEqual(sorted(balanced), sorted(balanced_def))

    def test_balance_problem(self):
        '''Tests get_elem_comp method for problematic reaction.'''
        problem = [('C144H238N2O57P2', -2.0, -1.0, 'C144H238N2O57P2'),
                   ('C86H142O9P', -1.0, -1.0, 'C86H142O9P'),
                   ('H', 1.0, 1.0, 'H'),
                   ('C150H248N2O62P2', -2.0, 1.0, 'C150H248N2O62P2'),
                   ('C80H131O4P', -2.0, 1.0, 'C80H131O4P')]

        is_balanced, was_balanced, balanced_def = chm_util.balance(
            problem, optional_comp=[])

        self.assertTrue(is_balanced)
        self.assertTrue(was_balanced)
        self.assertEqual(sorted(problem), sorted(balanced_def))

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
