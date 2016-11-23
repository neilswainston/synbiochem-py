'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=too-many-public-methods
import tempfile
import unittest

import synbiochem.utils.dna_utils as dna_utils


class Test(unittest.TestCase):
    '''Test class for dna_utils.'''

    def test(self):
        '''Tests round trip equality.'''
        dna1 = dna_utils.read('sbol.xml')
        tmp = tempfile.NamedTemporaryFile()
        dna_utils.write(dna1, tmp.name)
        dna2 = dna_utils.read(tmp.name)
        self.assertEqual(dna1, dna2)

if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
