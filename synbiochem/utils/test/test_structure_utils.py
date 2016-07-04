'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=too-many-public-methods
import unittest

import synbiochem.utils.structure_utils as struct_utils


class Test(unittest.TestCase):
    '''Test class for structure_utils.'''

    def test_get_seq_struct(self):
        '''Tests isclose method.'''
        seq_structs = struct_utils.get_seq_structs(['1DWI'])
        self.assertTrue(all([len(v[0]) == len(v[1])
                             for v in seq_structs.values()]))


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
