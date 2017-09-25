'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
from synbiochem.utils import seq_utils


class Mutation(object):
    '''Class to represent a mutation.'''

    def __init__(self, wt_res, pos, mut_res, typ='aa'):
        # Validate:
        alphabet = seq_utils.AA_CODES.values() if typ == 'aa' \
            else seq_utils.NUCLEOTIDES
        assert wt_res in alphabet
        assert isinstance(pos, (int, long)) and pos > 0
        assert mut_res in alphabet

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
