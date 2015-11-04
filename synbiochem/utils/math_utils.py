'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''


def isclose(aaa, bbb, rel_tol=1e-09, abs_tol=0.0):
    '''Compares floating point numbers.'''
    return abs(aaa - bbb) <= max(rel_tol * max(abs(aaa), abs(bbb)), abs_tol)
