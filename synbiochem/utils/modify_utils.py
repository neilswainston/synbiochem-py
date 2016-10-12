'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=bad-builtin
# pylint: disable=too-few-public-methods
from collections import defaultdict
from operator import mul
import itertools
import operator

from synbiochem.utils import sequence_utils
from synbiochem.utils.sequence_utils import CodonOptimiser
import Bio.Data.CodonTable as CodonTable


class CodonSelector(object):
    '''Class to optimise codon selection.'''

    def __init__(self, table_id=1):
        self.__codon_to_aa = \
            CodonTable.unambiguous_dna_by_id[table_id].forward_table
        self.__aa_to_codon = defaultdict(list)

        for codon, amino_acid in self.__codon_to_aa.iteritems():
            self.__aa_to_codon[amino_acid].append(codon)

        self.__codon_opt = {}

    def optimise_codons(self, query):
        '''Optimises codon selection.'''
        codons = self.__get_codons(set(query['aminoAcids'].upper()))
        combos = [combo for combo in itertools.product(*codons)]
        analyses = list(set([self.__analyse(combo, query['organism']['id'])
                             for combo in combos]))
        analyses.sort(key=operator.itemgetter(4, 3))
        return analyses

    def __get_codons(self, amino_acids):
        '''Gets codons for a list of amino acids.'''
        return [self.__aa_to_codon[amino_acid]
                for amino_acid in amino_acids]

    def __analyse(self, combo, tax_id):
        '''Analyses a combo, returning nucleotides, ambigous nucleotides,
        amino acids encodes, and number of variants.'''
        codon_opt = self.__get_codon_opt(tax_id)
        transpose = map(list, zip(*combo))
        nucls = [''.join(sorted(list(set(pos)))) for pos in transpose]
        ambig_codon = ''.join([sequence_utils.NUCL_CODES[nucl]
                               for nucl in nucls])
        amino_acids = defaultdict(list)

        for val in [(self.__codon_to_aa.get(''.join(combo), 'Stop'),
                     ''.join(combo),
                     codon_opt.get_codon_prob(''.join(combo)))
                    for combo in itertools.product(*nucls)]:
            amino_acids[val[0]].append(val[1:])

        return ambig_codon, \
            tuple(nucls), \
            tuple([(key, tuple(value))
                   for key, value in dict(amino_acids).iteritems()]), \
            len(amino_acids), \
            reduce(mul, [len(nucl) for nucl in nucls])

    def __get_codon_opt(self, tax_id):
        '''Gets the CodonOptimiser for the supplied taxonomy.'''
        if tax_id not in self.__codon_opt:
            self.__codon_opt[tax_id] = CodonOptimiser(tax_id)

        return self.__codon_opt[tax_id]
