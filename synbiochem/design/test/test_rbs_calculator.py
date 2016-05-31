'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=too-many-public-methods
import unittest

from synbiochem.design.rbs_calculator import RbsCalculator


class TestRbsCalculator(unittest.TestCase):
    '''Test class for RbsCalculator.'''

    def test_calc_kinetic_score(self):
        '''Tests calc_kinetic_score method.'''
        calc = RbsCalculator()
        m_rna = 'TTCTAGAGGGGGGATCTCCCCCCAAAAAATAAGAGGTACACATGACTAAAACTTTCA' + \
            'AAGGCTCAGTATTCCCACTGAG'
        start_pos = 41
        self.assertAlmostEqual(calc.calc_kinetic_score(m_rna, start_pos),
                               0.628571428571)

    def test_get_calc_dgs(self):
        '''Tests calc_dgs method.'''
        calc = RbsCalculator()
        m_rna = 'TTCTAGAGGGGGGATCTCCCCCCAAAAAATAAGAGGTACACATGACTAAAACTTTCA' + \
            'AAGGCTCAGTATTCCCACTGAG'
        r_rna = 'acctcctta'
        dgs = calc.calc_dgs(r_rna, m_rna)
        self.assertEqual(dgs[0], [41, 74])
        self.assertAlmostEqual(dgs[1][0], -8.070025836938619)
        self.assertAlmostEqual(dgs[1][1], 3.312588580920539)


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
