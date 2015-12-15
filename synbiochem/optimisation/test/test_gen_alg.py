'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=no-member
# pylint: disable=too-many-public-methods
import unittest
import synbiochem.optimisation.gen_alg as gen_alg


class TestChromosome(unittest.TestCase):

    def test_breed(self):
        length = 16
        chrom1 = gen_alg.Chromosome(0, length)
        chrom2 = gen_alg.Chromosome(2 ** length - 1, length)

        for _ in range(1024):
            chrom1.breed(chrom2)

        self.assertEqual(
            chrom1.get_chromosome() + chrom2.get_chromosome(), 2 ** length - 1)
