'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=too-many-public-methods
import unittest

import synbiochem.utils.math_utils as math_utils


class Test(unittest.TestCase):
    '''Test class for math_utils.'''

    def test_isclose(self):
        '''Tests isclose method.'''
        self.assertTrue(math_utils.isclose(18.015054684, 18.015054684))
        self.assertFalse(math_utils.isclose(18.015054684, 18.015054))


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
