'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=no-member
# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals
import copy
from itertools import product
import math
import re
import uuid

from Bio.Restriction import RestrictionBatch
from Bio.Seq import Seq


SO_PROM = 'http://purl.obolibrary.org/obo/SO_0000167'
SO_RBS = 'http://purl.obolibrary.org/obo/SO_0000139'
SO_CDS = 'http://purl.obolibrary.org/obo/SO_0000316'
SO_PART = 'http://purl.obolibrary.org/obo/SO_0000804'
SO_ASS_COMP = 'http://purl.obolibrary.org/obo/SO_0000143'
SO_PLASMID = 'http://purl.obolibrary.org/obo/SO_0000637'
SO_DESIGNED = 'http://purl.obolibrary.org/obo/SO_0000546'
SO_RANDOM = 'http://purl.obolibrary.org/obo/SO_0000449'

_RDF_NS = 'http://www.w3.org/1999/02/22-rdf-syntax-ns#'
_NS = {'ns': 'http://sbols.org/v1#',
       'rdf': _RDF_NS}


class DNA(dict):
    '''Class to represent a DNA object.'''

    def __init__(self, disp_id=None, name=None, desc=None, typ=None, seq=None,
                 start=1, end=float('NaN'), forward=True,
                 features=None, children=None, options=None,
                 links=None, parameters=None, temp_params=None):

        if seq is None:
            seq = ''

        dict.__init__(self)

        self.update({'disp_id': disp_id
                     if disp_id is not None else str(uuid.uuid4()),
                     'seq': re.sub(r'[\s]', '', seq.upper()),
                     'name': name,
                     'desc': desc,
                     'typ': typ,
                     'start': start,
                     'end': end if not math.isnan(end)
                     else start + len(re.sub(r'[\s]', '', seq.upper())) - 1,
                     'forward': forward,
                     'features': [] if features is None else features,
                     'children': [] if children is None else children,
                     'options': [] if options is None else options,
                     'links': [] if links is None else links,
                     'parameters': {} if parameters is None else parameters,
                     'temp_params': {} if temp_params is None else temp_params,
                     })

    def set_seq(self, seq):
        '''Sets sequence.'''
        self['seq'] = re.sub(r'[\s]', '', seq.upper())
        self['end'] = self['start'] + len(seq) - 1

    def copy(self):
        '''Copies the DNA object (making a deep copy).'''
        copy_dict = copy.deepcopy(self)
        copy_dna = DNA(**copy_dict)
        return copy_dna

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __repr__(self):
        return self.get('name', None)


def get_dna(dct):
    '''Factory method for constructing DNA object from dict.'''
    dna = DNA(**dct)
    dna['features'] = [get_dna(feat) for feat in dna['features']]
    dna['children'] = [get_dna(child) for child in dna['children']]
    dna['options'] = [get_dna(opt) for opt in dna['options']]
    return dna


def concat(dnas):
    '''Concatenates a list of DNA objects into a single DNA object.'''
    concat_dna = dnas[0].copy()
    concat_dna['disp_id'] = str(uuid.uuid4())

    for dna in dnas[1:]:
        concat_dna = add(concat_dna, dna)

    return concat_dna


def apply_restricts(dna, restricts, circular=False):
    '''Applies restriction site cleavage to forward and reverse strands.'''
    out_dnas = [dna]

    for restrict in restricts:
        batch = RestrictionBatch()
        batch.add(str(restrict))
        restrict = batch.get(str(restrict))
        out_dnas = _apply_restrict_to_dnas(out_dnas, restrict, circular)

    return out_dnas


def add(dna1, dna2):
    '''Adds two DNA objects together.'''
    # Add names, etc.
    dna1.disp_id = str(uuid.uuid4())
    dna1['name'] = _concat([dna1['name'], dna2['name']])
    dna1['desc'] = _concat([dna1['desc'], dna2['desc']])

    # Add sequences:
    orig_seq_len = len(dna1['seq'])
    dna1['seq'] += dna2['seq']
    dna1['end'] = len(dna1['seq'])

    # Update parameters:
    for key, value in dna2['parameters'].items():
        param = dna1['parameters'].get(key, None)

        if param is None:
            dna1['parameters'][key] = value
        elif isinstance(param, list):
            param.append(value)
            dna1['parameters'][key] = param
        else:
            dna1['parameters'][key] = [param, value]

    # Update features:
    for feature in dna2.get('features', []):
        feature = feature.copy()
        feature['start'] += orig_seq_len
        feature['end'] += orig_seq_len
        dna1['features'].append(feature)

    return dna1


def expand(dna):
    '''Expand a DNA object containing options into distinct DNA objects.'''
    features = dna.pop('features')
    dna['features'] = []

    dnas = [dna]

    for feature in features:
        if feature['options']:
            additions = feature['options']
        elif feature['seq']:
            additions = [feature]
        else:
            continue

        for addition in additions:
            addition['features'].append(addition.copy())

        dnas = [add(dna.copy(), addition)
                for dna, addition in product(dnas, additions)
                if len(addition['seq'])]

    return dnas


def _concat(strs):
    '''Concatenates non-None strings.'''
    return ' - '.join([string for string in strs
                       if string is not None and len(string)])


def _apply_restrict_to_dnas(dnas, restrict, circular):
    '''Applies restriction site cleavage to forward and reverse strands.'''
    out_dnas = []

    for dna in dnas:
        restrict.catalyse(Seq(dna['seq']), linear=not circular)
        start = 0
        seqs = []

        for result in restrict.results:
            end = result - 1
            seqs.append((dna['seq'][start:end], start, end))
            start = end

        seqs.append((dna['seq'][start:], start, len(dna['seq'])))

        if circular:
            seqs = [(seqs[-1][0] + seqs[0][0], seqs[-1][1], seqs[0][2])] + \
                seqs[1:-1]

        for seq in seqs:
            out_dnas.append(_get_concat_dna(dna, seq[0], seq[1], seq[2]))

    return out_dnas


def _apply_restrict_to_seq(seq, restrict):
    '''Applies restriction site cleavage to a sequence.'''
    sub_seqs = [(match.group(0), match.start())
                for match in re.finditer(restrict, seq)]
    end = sub_seqs[0][1] if sub_seqs else len(seq)
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
