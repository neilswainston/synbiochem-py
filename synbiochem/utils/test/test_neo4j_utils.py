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
                 False,
                 13,
                 56.3,
                 [23.5, 'Tree', False]],
                [[1.6, 3.77, -374.8e-17],
                 ['another', 'random', 'string'],
                 [47, -31],
                 [True, True, False, True],
                 132.124,
                 -314,
                 'another random string',
                 True,
                 34,
                 78.3,
                 [-123.5, 6, 'Trees', True]]]

        columns = ['float_array',
                   'string_array',
                   'int_array',
                   'boolean_array',
                   'float',
                   'int',
                   'string',
                   'boolean',
                   ':START_ID',
                   ':END_ID(Label)',
                   ':LABEL']

        df = pd.DataFrame(data, columns=columns)
        new_df = neo4j_utils.type_df(df, array_delimiter='|')

        expected = ['float_array:float[]',
                    'string_array:string[]',
                    'int_array:int[]',
                    'boolean_array:boolean[]',
                    'float:float',
                    'int:int',
                    'string:string',
                    'boolean:boolean',
                    ':START_ID',
                    ':END_ID(Label)',
                    ':LABEL']

        self.assertEqual(sorted((list(new_df.columns))), sorted(expected))


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
