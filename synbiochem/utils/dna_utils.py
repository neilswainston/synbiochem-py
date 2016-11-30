'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=no-member
# pylint: disable=too-many-arguments
import copy
import math
import re
import uuid

from synbiochem.utils import seq_utils


SO_CDS = 'http://purl.obolibrary.org/obo/SO_0000316'
SO_PROM = 'http://purl.obolibrary.org/obo/SO_0000167'

_RDF_NS = 'http://www.w3.org/1999/02/22-rdf-syntax-ns#'
_NS = {'ns': 'http://sbols.org/v1#',
       'rdf': _RDF_NS}


class DNA(dict):
    '''Class to represent a DNA object.'''

    def __init__(self, disp_id=None, name=None, desc=None, typ=None, seq=None,
                 start=float('NaN'), end=float('NaN'), forward=True,
                 features=None):

        dict.__init__(self)

        if math.isnan(end) and (not seq or len(seq) == 0):
            raise ValueError('Unable to determine sequence length')

        self.update({'disp_id': disp_id
                     if disp_id is not None else str(uuid.uuid4()),
                     'seq': seq,
                     'name': name,
                     'desc': desc,
                     'typ': typ,
                     'start': start if not math.isnan(start) else 1,
                     'end': end if not math.isnan(end) else len(seq),
                     'forward': forward,
                     'features': [] if features is None else features
                     })

    def set_seq(self, seq):
        '''Sets sequence.'''
        self['sequence'] = seq

        if math.isnan(self['end']):
            self['end'] = len(seq)

    def copy(self):
        '''Copies the DNA object (making a deep copy).'''
        copy_dict = copy.deepcopy(self)
        copy_dna = DNA(**copy_dict)
        return copy_dna

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __repr__(self):
        return self['name']


def concat(dnas):
    '''Concatenates a list of DNA objects into a single DNA object.'''
    concat_dna = dnas[0].copy()

    for dna in dnas[1:]:
        concat_dna = _add(concat_dna, dna)

    return concat_dna


def apply_restricts(dna, restricts, circular=False):
    '''Applies restriction site cleavage to forward and reverse strands.'''
    out_dnas = [dna]

    for restrict in restricts:
        out_dnas = _apply_restrict_to_dnas(out_dnas, restrict)

    if circular and len(out_dnas) > 1:
        return [concat([out_dnas[-1], out_dnas[0]])] + out_dnas[1:-1]
    else:
        return out_dnas


def _add(dna1, dna2):
    '''Adds two DNA objects together.'''
    # Add names, etc.
    dna1.disp_id = str(uuid.uuid4())
    dna1.name = _concat([dna1['name'], dna2['name']])
    dna1.description = _concat([dna1['desc'], dna2['desc']])

    # Add sequences:
    orig_seq_len = len(dna1['seq'])
    dna1['seq'] += dna2['seq']

    # Update features:
    for feature in dna2['features']:
        feature = feature.copy()
        feature['start'] += orig_seq_len
        feature['end'] += orig_seq_len
        dna1['features'].append(feature)

    return dna1


def _concat(strs):
    '''Concatenates non-None strings.'''
    return ' - '.join([string for string in strs if string is not None])


def _apply_restrict_to_dnas(dnas, restrict):
    '''Applies restriction site cleavage to forward and reverse strands.'''
    out_dnas = []

    for dna in dnas:
        parent_seq = dna['seq']

        for forw in _apply_restrict_to_seq(parent_seq, restrict):
            for rev in _apply_restrict_to_seq(seq_utils.get_rev_comp(forw[0]),
                                              restrict):
                rev_comp = seq_utils.get_rev_comp(rev[0])
                start = forw[1] + len(forw[0]) - rev[1] - len(rev[0])
                end = start + len(rev_comp)
                out_dnas.append(_get_concat_dna(dna, rev_comp, start, end))

    return out_dnas


def _apply_restrict_to_seq(seq, restrict):
    '''Applies restriction site cleavage to a sequence.'''
    sub_seqs = [(match.group(0), match.start())
                for match in re.finditer(restrict, seq)]
    end = sub_seqs[0][1] if len(sub_seqs) > 0 else len(seq)
    return [(seq[:end], 0)] + sub_seqs


def _get_concat_dna(parent_dna, seq, start, end):
    '''Returns a DNA object from the supplied subsequence from a parent DNA
    object.'''
    if parent_dna['seq'] == seq:
        disp_id = parent_dna['disp_id']
        frag_str = ''
    else:
        disp_id = str(uuid.uuid4())
        frag_str = ' [' + str(start) + ':' + str(end) + ']'

    dna = DNA(disp_id=disp_id,
              name=parent_dna['name'] + frag_str,
              desc=parent_dna['desc'] + frag_str,
              seq=seq)

    for feature in parent_dna['features']:
        if feature['start'] >= start and feature['end'] <= end:
            copy_feature = feature.copy()
            copy_feature['start'] -= start
            copy_feature['end'] -= start
            dna['features'].append(copy_feature)

    return dna
