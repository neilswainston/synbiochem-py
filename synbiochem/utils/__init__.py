'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import re


def pairwise(iterable):
    '''s -> (s0,s1), (s1,s2), (s2, s3), ...'''
    return [(iterable[i], iterable[i + 1]) for i in range(len(iterable) - 1)]


def sort(sortable):
    ''''Sort the given iterable in the way that humans expect.'''
    def convert(text):
        '''Convert char to int if digit.'''
        return int(text) if text.isdigit() else text

    def alphanum_key(key):
        '''Generate alphanumeric key.'''
        return [convert(c) for c in re.split('([0-9]+)', key)]

    return sorted(sortable, key=alphanum_key)
