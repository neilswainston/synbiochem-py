'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=wrong-import-order
import unittest

from synbiochem.utils import neo4j_utils

import pandas as pd


class Test(unittest.TestCase):
    '''Test class for neo4j_utils.'''

    def test_type_files(self):
        '''Tests type_files method.'''
        data = [[[1.0, 3.7, -34.8e-17],
                 ['some', 'random', 'string'],
                 [4, -3],
                 [True, False, True],
                 32.124,
                 -34,
                 'random string',
                 False],
                [[1.6, 3.77, -374.8e-17],
                 ['another', 'random', 'string'],
                 [47, -31],
                 [True, True, False, True],
                 132.124,
                 -314,
                 'another random string',
                 True]]

        columns = ['float_array',
                   'string_array',
                   'int_array',
                   'boolean_array',
                   'float',
                   'int',
                   'string',
                   'boolean']

        df = pd.DataFrame(data, columns=columns)
        new_dfs = neo4j_utils.type_dfs([df])

        expected = ['float_array:float[]',
                    'string_array:string[]',
                    'int_array:int[]',
                    'boolean_array:boolean[]',
                    'float:float',
                    'int:int',
                    'string:string',
                    'boolean:boolean']

        self.assertEqual(sorted((list(new_dfs[0].columns))), sorted(expected))


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
