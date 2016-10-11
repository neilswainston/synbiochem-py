'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=bad-builtin
# pylint: disable=too-few-public-methods
from collections import defaultdict, Counter
from operator import mul
import itertools
import operator
import sys

from synbiochem.utils import sequence_utils
import Bio.Data.CodonTable as CodonTable


class CodonSelector(object):
    '''Class to optimise codon selection.'''

    def __init__(self, table_id=1):
        self.__codon_to_aa = \
            CodonTable.unambiguous_dna_by_id[table_id].forward_table
        self.__aa_to_codon = defaultdict(list)

        for codon, amino_acid in self.__codon_to_aa.iteritems():
            self.__aa_to_codon[amino_acid].append(codon)

    def optimise_codon(self, amino_acids):
        '''Optimises codon selection.'''
        codons = self.__get_codons(set(amino_acids))
        combos = [combo for combo in itertools.product(*codons)]
        analyses = list(set([self.__analyse(combo) for combo in combos]))
        analyses.sort(key=operator.itemgetter(4, 3))
        return analyses

    def __get_codons(self, amino_acids):
        '''Gets codons for a list of amino acids.'''
        return [self.__aa_to_codon[amino_acid]
                for amino_acid in amino_acids]

    def __analyse(self, combo):
        '''Analyses a combo, returning nucleotides, ambigous nucleotides,
        amino acids encodes, and number of variants.'''
        transpose = map(list, zip(*combo))
        nucls = [''.join(sorted(list(set(pos)))) for pos in transpose]
        ambig_codon = ''.join([sequence_utils.NUCL_CODES[nucl]
                               for nucl in nucls])
        amino_acids = dict(Counter([self.__codon_to_aa.get(''.join(combo),
                                                           'Stop')
                                    for combo in itertools.product(*nucls)]))

        return ambig_codon, \
            '[' + ']['.join(nucls) + ']', \
            frozenset(amino_acids.items()), \
            len(amino_acids), \
            reduce(mul, [len(nucl) for nucl in nucls])


def main(args):
    '''main method.'''
    analyses = CodonSelector().optimise_codon(args[0])

    for analysis in analyses:
        print '\t'.join([str(val) for val in analysis])


if __name__ == '__main__':
    main(sys.argv[1:])
