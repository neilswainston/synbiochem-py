'''
synbiochem (c) University of Manchester 2016

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''


def get_well(idx, rows=8, columns=12, col_ord=True):
    '''Map idx to well'''
    if idx < 0 or idx >= rows * columns:
        raise ValueError('Index %idx out of range' % idx)

    return _get_well_col(idx, rows) if col_ord \
        else _get_well_row(idx, columns)


def _get_well_col(idx, rows):
    '''Map idx to well, column ordered.'''
    return chr(ord('A') + (idx % rows)) + str(idx / rows + 1)


def _get_well_row(idx, columns):
    '''Map idx to well, row ordered.'''
    return chr(ord('A') + (idx / columns)) + str(idx % columns + 1)
