'''
synbiochem (c) University of Manchester 2016

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''


def get_well(idx, rows=8, columns=12):
    '''Map idx to well'''
    if idx < 0 or idx >= rows * columns:
        raise ValueError('Index %idx out of range' % idx)

    return chr(ord('A') + (idx / rows)) + str(idx % rows + 1)
