'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=useless-object-inheritance
import re

from synbiochem.utils import seq_utils


class Mutation(object):
    '''Class to represent a mutation.'''

    def __init__(self, wt_res, pos, mut_res, typ='aa'):
        # Validate:
        alphabet = list(seq_utils.AA_CODES.values()) if typ == 'aa' \
            else seq_utils.NUCLEOTIDES
        assert wt_res in alphabet
        assert isinstance(pos, int) and pos > 0
        assert mut_res in alphabet + ['-']  # Consider deletions

        self.__wt_res = wt_res
        self.__pos = pos
        self.__mut_res = mut_res

    def get_wt_res(self):
        '''Gets wt residue.'''
        return self.__wt_res

    def get_pos(self):
        '''Gets position.'''
        return self.__pos

    def get_mut_res(self):
        '''Gets mutation residue.'''
        return self.__mut_res

    def __repr__(self):
        return self.__wt_res + str(self.__pos) + self.__mut_res

    def __cmp__(self, other):
        pos_diff = self.__pos - other.get_pos()

        if pos_diff:
            return pos_diff

        return ord(self.__mut_res) - ord(other.get_mut_res())

    def __eq__(self, other):
        return self.__pos == other.get_pos() and \
            self.__mut_res == other.get_mut_res()

    def __lt__(self, other):
        if self.__pos == other.get_pos():
            return ord(self.__mut_res) < ord(other.get_mut_res())

        return self.__pos < other.get_pos()


def parse_mut_str(mut_str):
    '''Parse mutation string.'''
    return [Mutation(mut[0], int(mut[1]), mut[2])
            for mut in [re.compile(r'(\d+)').split(mutation)
                        for mutation in mut_str.split()]]


def get_mutations(wt_seq, mut_seq):
    '''Get Mutations.'''
    assert len(wt_seq) == len(mut_seq)

    mutations = []

    for (pos, aas) in enumerate(zip(wt_seq, mut_seq)):
        if aas[0] != aas[1]:
            mutations.append(Mutation(aas[0], pos + 1, aas[1]))

    return mutations


def apply_mutations(seq, mutations):
    '''Applies mutations to sequence.'''
    seq = list(seq)
    offset = 1

    for mutation in mutations:
        if mutation.get_wt_res() != seq[mutation.get_pos() - offset]:
            err = 'Invalid mutation at position %d. ' % mutation.get_pos() + \
                'Amino acid is %s ' % seq[mutation.get_pos() - offset] + \
                'but mutation is of %s.' % mutation.get_wt_res()

            raise ValueError(err)

        if mutation.get_mut_res() == '-':
            # Deletion:
            del seq[mutation.get_pos() - offset]
            offset += 1
        else:
            # Mutation:
            seq[mutation.get_pos() - offset] = mutation.get_mut_res()

    return ''.join(seq)
