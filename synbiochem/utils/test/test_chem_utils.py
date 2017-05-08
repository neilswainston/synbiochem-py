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
                         [['Fe4S', -5.6], ['water', 3.2], ['SiO2', 17.8]])

    def test_parse_equation_error(self):
        '''Tests parse_equation method (with error).'''
        eqn = '5.6 Fe4S + -3.2 water = n+m SiO2'
        self.assertRaises(ValueError, chm_util.parse_equation, eqn)


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
