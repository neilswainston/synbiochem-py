'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=no-member
# pylint: disable=superfluous-parens
# pylint: disable=too-many-arguments
import re
import uuid
from xml.etree import ElementTree

from synbiochem.utils.dna_utils import DNA


_RDF_NS = 'http://www.w3.org/1999/02/22-rdf-syntax-ns#'
_NS = {'ns': 'http://sbols.org/v1#',
       'rdf': _RDF_NS}


def read(filename):
    '''Parses SBOL v1.'''
    tree = ElementTree.parse(filename)
    root = tree.getroot()

    dna_comp = root.find('ns:DnaComponent', _NS)
    dna_seq = dna_comp.find('ns:dnaSequence', _NS)
    dna_seq = dna_seq.find('ns:DnaSequence', _NS)

    params = _read_dna_comp(dna_comp)
    params.update({'seq': dna_seq.find('ns:nucleotides', _NS).text})
    dna = DNA(**params)

    for annot in dna_comp.findall('ns:annotation', _NS):
        _read_annot(dna, annot)

    return dna


def write(dna, filename=None):
    '''Writes a Dna object to SBOL v1.'''
    root = ElementTree.Element('ns2:RDF', {'xmlns': 'http://sbols.org/v1#',
                                           'xmlns:ns2': _RDF_NS})

    dna_comp = _write_dna_comp(root, dna)

    dna_seq = _write(dna_comp, 'dnaSequence')
    dna_seq = _write(dna_seq, 'DnaSequence', _get_about())
    _write(dna_seq, 'nucleotides', text=dna['seq'])

    for feature in dna['features']:
        annot = _write(dna_comp, 'annotation')
        annot = _write(annot, 'SequenceAnnotation', _get_about())
        _write(annot, 'bioStart', text=str(feature['start']))
        _write(annot, 'bioEnd', text=str(feature['end']))
        _write(annot, 'strand', text='+' if feature['forward'] else '-')
        sub_comp = _write(annot, 'subComponent')
        _write_dna_comp(sub_comp, feature)

    sbol = ElementTree.tostring(root)

    with open(filename, 'wb') as outfile:
        outfile.write(sbol)

    return sbol


def _read_dna_comp(dna_comp):
    '''Read DNAComponent node.'''
    disp_id = _find_text(dna_comp, 'ns:displayId')
    name = _find_text(dna_comp, 'ns:name')
    desc = _find_text(dna_comp, 'ns:description')
    typ_node = dna_comp.find('rdf:type', _NS)
    typ = typ_node.attrib['{' + _RDF_NS + '}resource'] \
        if typ_node is not None else None

    if not re.match('^[\\w-]+$', disp_id):
        disp_id = str(uuid.uuid4())

    return {'disp_id': disp_id, 'name': name, 'desc': desc, 'typ': typ}


def _read_annot(dna, annot):
    '''Reads annotation node.'''
    seq_annot = annot.find('ns:SequenceAnnotation', _NS)
    sub_comp = seq_annot.find('ns:subComponent', _NS)
    dna_comp = sub_comp.find('ns:DnaComponent', _NS)

    params = _read_dna_comp(dna_comp)

    start = _find_text(seq_annot, 'ns:bioStart')

    if start:
        params.update({'start': int(start)})

    end = _find_text(seq_annot, 'ns:bioEnd')

    if end:
        params.update({'end': int(end)})

    forward = _find_text(seq_annot, 'ns:strand')

    if forward:
        params.update({'forward': forward == '+'})

    # Tests due to ICE eccentricities...
    try:
        feat = DNA(**params)

        pos = (feat['start'], feat['end'], feat['forward'])

        # Prevents cases where features are duplicated in the SBOL:
        if pos not in [(feature['start'], feature['end'], feature['forward'])
                       for feature in dna['features']]:
            dna['features'].append(feat)

    except ValueError:
        # Prevents cases with no end position and no sequence:
        print('Ignoring invalid feature.')


def _find_text(parent, field):
    '''Finds text from node.'''
    node = parent.find(field, _NS)
    return None if node is None else node.text


def _write_dna_comp(parent, dna):
    '''Write DNAComponent node.'''
    dna_comp = ElementTree.SubElement(parent, 'DnaComponent',
                                      _get_about(dna['disp_id']))

    if dna['typ']:
        _write(dna_comp, 'ns2:type', {'ns2:resource': dna['typ']})

    _write(dna_comp, 'displayId', text=dna['disp_id'])

    if dna['name']:
        _write(dna_comp, 'name', text=dna['name'])

    if dna['desc']:
        _write(dna_comp, 'description', text=dna['desc'])

    return dna_comp


def _write(parent, name, params=None, text=None):
    '''Writes a node.'''
    if params is None:
        params = {}

    node = ElementTree.SubElement(parent, name, params)
    node.text = text
    return node


def _get_about(uid=None):
    '''Gets about attributes.'''
    if uid is None:
        uid = str(uuid.uuid4())

    return {'ns2:about': 'https://www.synbiochem.co.uk#' + uid}
