'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=bad-builtin
# pylint: disable=too-few-public-methods
from collections import defaultdict
import itertools
import operator

from mock.mock import self

from synbiochem.utils import sequence_utils
from synbiochem.utils.sequence_utils import CodonOptimiser
import Bio.Data.CodonTable as CodonTable


class CodonSelector(object):
    '''Class to optimise codon selection.'''

    def __init__(self, table_id=1, wanted_pen=3, degen_pen=-1, unwanted_pen=-3,
                 stop_pen=-10):
        self.__codon_to_aa = \
            CodonTable.unambiguous_dna_by_id[table_id].forward_table
        self.__aa_to_codon = defaultdict(list)

        self.__wanted_pen = wanted_pen
        self.__degen_pen = degen_pen
        self.__unwanted_pen = unwanted_pen
        self.__stop_pen = stop_pen

        for codon, amino_acid in self.__codon_to_aa.iteritems():
            self.__aa_to_codon[amino_acid].append(codon)

        self.__codon_opt = {}

    def optimise_codons(self, query):
        '''Optimises codon selection.'''
        req_amino_acids = set(query['aminoAcids'].upper())
        codons = self.__get_codons(req_amino_acids)
        combos = [combo for combo in itertools.product(*codons)]
        analyses = [self.__analyse(combo, query['organism']['id'],
                                   req_amino_acids)
                    for combo in combos]
        analyses = list(set([codon for analysis in analyses
                             for codon in analysis]))
        analyses.sort(key=operator.itemgetter(5), reverse=True)
        return analyses

    def __get_codons(self, amino_acids):
        '''Gets codons for a list of amino acids.'''
        return [sequence_utils.CODONS[amino_acid]
                for amino_acid in amino_acids]

    def __analyse(self, combo, tax_id, req_amino_acids):
        '''Analyses a combo, returning nucleotides, ambiguous nucleotides,
        amino acids encodes, and number of variants.'''
        codon_opt = self.__get_codon_opt(tax_id)
        transpose = [sorted(list(term))
                     for term in map(set, zip(*combo))]

        nucls = [[''.join(sorted(list(set(pos))))]
                 for pos in transpose[:2]] + \
            [_optimise_pos_3(transpose[2])]

        ambig_codons = [[''.join([sequence_utils.NUCL_CODES[term]
                                  for term in cdn]),
                         cdn,
                         [''.join(c)
                          for c in itertools.product(*cdn)]]
                        for cdn in itertools.product(*nucls)]

        for ambig_codon in ambig_codons:
            num_codons = len(ambig_codon[2])
            amino_acids = defaultdict(list)

            for codon in ambig_codon[2]:
                a_a = self.__codon_to_aa.get(codon, 'Stop')
                amino_acids[a_a].append((codon,
                                         codon_opt.get_codon_prob(codon)))

            ambig_codon[2] = tuple([(key,
                                     tuple(sorted(value,
                                                  key=lambda prob: prob[1],
                                                  reverse=True)))
                                    for key, value in amino_acids.iteritems()])

            ambig_codon.append(len(amino_acids))
            ambig_codon.append(num_codons)
            ambig_codon.append(self.__score(amino_acids, req_amino_acids))

        return [tuple(ambig_codon) for ambig_codon in ambig_codons]

    def __get_codon_opt(self, tax_id):
        '''Gets the CodonOptimiser for the supplied taxonomy.'''
        if tax_id not in self.__codon_opt:
            self.__codon_opt[tax_id] = CodonOptimiser(tax_id)

        return self.__codon_opt[tax_id]

    def __score(self, amino_acids, req_amino_acids):
        scores = [value[1] * (self.__wanted_pen
                              if idx == 0 and amino_acid in req_amino_acids
                              else (self.__stop_pen if amino_acid is 'Stop'
                                    else self.__unwanted_pen))
                  for amino_acid, values in amino_acids.iteritems()
                  for idx, value in enumerate(values)]

        return sum(scores)


def _optimise_pos_3(options):
    options = list(set([tuple(sorted(set(opt)))
                        for opt in itertools.product(*options)]))
    options.sort(key=len)
    return [''.join(opt) for opt in options]
