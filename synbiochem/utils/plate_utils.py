'''
synbiochem (c) University of Manchester 2016

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import math


def get_well_pos(pos, rows=12, columns=8):
    '''Converts a position, pos, to a well position on a microtitre plate.'''
    if pos >= rows * columns:
        raise ValueError('Position ' + str(pos) + ' too large for plate.')

    return chr(int(pos % columns) + ord('A')) + \
        str(int(math.floor(float(pos) / columns)) + 1)
