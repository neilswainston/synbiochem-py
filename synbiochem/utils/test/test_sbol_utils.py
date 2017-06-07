'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=no-member
# pylint: disable=too-many-public-methods
import os
import tempfile
import unittest

from synbiochem.utils import sbol_utils


class Test(unittest.TestCase):
    '''Test class for dna_utils.'''

    def test(self):
        '''Tests round trip equality.'''
        directory = os.path.dirname(os.path.realpath(__file__))
        dna1 = sbol_utils.read(os.path.join(directory, 'sbol.xml'))
        dna2 = round_trip(dna1)
        self.assertEqual(dna1, dna2)


def round_trip(dna):
    '''Writes / reads DNA object, via SBOL export / import.'''
    tmp = tempfile.NamedTemporaryFile()
    sbol_utils.write(dna, tmp.name)
    return sbol_utils.read(tmp.name)


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
