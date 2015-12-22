'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
from matplotlib.colors import LinearSegmentedColormap
import numpy
import os
import pylab
import random
import re
import scipy.spatial

from Bio import SeqUtils
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder

import regex
import synbiochem.utils.io_utils as io_utils


_DIR = 'structure_utils'
_PDB_DIR = 'pdb'


def get_pdb_ids(max_ids=None, local_only=False):
    '''Returns PDB ids.'''
    if local_only:
        # Returns local PDB ids.
        pdb_dir = os.path.join(os.path.expanduser('~'), _DIR, _PDB_DIR)
        ids = [filename[:-4].upper()
               for _, _, files in os.walk(pdb_dir)
               for filename in files if filename.endswith('.pdb')]
    else:
        # Returns all PDB ids.
        source_url = 'http://www.uniprot.org/uniprot/?query=database:pdb' \
            + '&format=tab&columns=id,database(PDB)'
        target_filename = os.path.join(os.path.expanduser('~'), _DIR,
                                       'pdb_ids.txt')

        with open(io_utils.get_file(source_url, target_filename)) as fle:
            ids = [x for line in fle
                   for x in line.split()[1].split(';')
                   if len(x) > 0 and x != 'Cross-reference']

    return ids if max_ids is None \
        else random.sample(ids, min(len(ids), max_ids))


def get_seq_struct(pdb_ids):
    '''Returns sequence and structure.'''
    seq_struct = {}
    pdb_ids = sorted(pdb_ids)
    pdb_id = pdb_ids.pop(0)
    pdb_id_fasta = '>' + pdb_id
    in_field = False
    tokens = None
    str_data = ''

    source_url = 'http://www.rcsb.org/pdb/files/ss.txt'
    target_filename = os.path.join(os.path.expanduser('~'), _DIR,
                                   'ss.txt')

    with open(io_utils.get_file(source_url, target_filename)) as fle:
        for line in fle:
            if line.startswith(pdb_id_fasta):
                if in_field:
                    if tokens[:2] not in seq_struct:
                        seq_struct[tokens[:2]] = [None, None]

                    seq_struct[tokens[:2]][0 if tokens[2] == 'sequence'
                                           else 1] = str_data
                    str_data = ''

                tokens = tuple(re.split('>|:', line.strip())[1:])
                in_field = True

            elif in_field and line.startswith('>'):
                if tokens[:2] not in seq_struct:
                    seq_struct[tokens[:2]] = [None, None]

                seq_struct[tokens[:2]][0 if tokens[2] == 'sequence'
                                       else 1] = str_data
                str_data = ''

                if len(pdb_ids) == 0:
                    break
                else:
                    pdb_id = pdb_ids.pop(0)
                    pdb_id_fasta = '>' + pdb_id
                    in_field = False
                    tokens = None
                    str_data = ''

            elif in_field:
                str_data += line[:-1]

    seq_struct = {key: val for key, val in seq_struct.iteritems()
                  if val[0] is not None}

    return seq_struct


def get_pep_struct(outfilepath, length):
    '''Returns sequence and structure.'''
    seq = None
    pdb_id = None
    chain = None
    str_data = ''

    source_url = 'http://www.rcsb.org/pdb/files/ss.txt'
    target_filename = os.path.join(os.path.expanduser('~'), _DIR,
                                   'ss.txt')

    count = 0

    with open(io_utils.get_file(source_url, target_filename)) as infile, \
            open(outfilepath, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                tokens = tuple(re.split('>|:', line.strip())[1:])

                if tokens[2] == 'sequence':
                    if seq is not None:
                        for i in range(len(seq) - length + 1):
                            count += 1
                            outfile.write('\t'.join([seq[i:i + length],
                                                     str_data[i:i + length],
                                                     pdb_id, chain,
                                                     str(i),
                                                     str(i + length)]) + '\n')
                else:
                    seq = str_data

                str_data = ''
                pdb_id = tokens[0]
                chain = tokens[1]
            else:
                str_data += line[:-1]

    print count


def get_sequences(pdb_id, chain=None):
    '''Gets the sequences in a PDB file.'''
    return [SeqUtils.seq1(''.join([residue.get_resname()
                                   for residue in chn
                                   if 'CA' in residue.child_dict]))
            for chn in get_structure(pdb_id).get_chains()
            if chain is None or chain == chn.get_id()]


def get_structure(pdb_id):
    '''Returns a PDB structure.'''
    source_url = 'http://www.rcsb.org/pdb/files/' + pdb_id + '.pdb'
    target_filename = os.path.join(os.path.expanduser('~'), _DIR, _PDB_DIR,
                                   pdb_id + '.pdb')

    with open(io_utils.get_file(source_url, target_filename)) as pdb_file:
        parser = PDBParser(QUIET=True)
        return parser.get_structure(pdb_id, pdb_file.name)


def calc_proximities(pdb_id):
    '''Calculates residue proximities from PDB file.'''
    structure = get_structure(pdb_id)
    chains = [c for c in structure.get_chains()]

    coords = [[residue.child_dict['CA'].get_coord()
               for residue in chn
               if 'CA' in residue.child_dict]
              for chn in chains]

    return [scipy.spatial.distance.cdist(coord, coord, 'euclidean')
            if len(coord) > 0
            else None
            for coord in coords]


def get_phi_psi_data(pdb_id, chain=None):
    '''Gets phi and phi angle data.'''
    builder = PPBuilder()
    return [polypep.get_phi_psi_list()
            for model in get_structure(pdb_id)
            for chn in model
            if chain is None or chain == chn.get_id()
            for polypep in builder.build_peptides(chn)]


def plot_proximities(pdb_id):
    '''Plots proximity plot(s).'''
    all_proximities = calc_proximities(pdb_id)

    plot_format = 'png'

    for idx, proximities in enumerate(all_proximities):
        name = pdb_id + '_' + str(idx + 1)
        _plot(proximities, name + '.' + plot_format, plot_format,
              name + ' proximity plot')


def sample_seqs(sample_size, struct_patt, local_only=False):
    '''Sample sequence and structure data.'''
    seqs = []

    while len(seqs) < sample_size:
        # Get more PDB ids than necessary, as not all PDB entries will contain
        # desired pattern:
        pdb_ids = get_pdb_ids(sample_size * 10, local_only)
        seq_structs = get_seq_struct(pdb_ids)
        matches = []

        for pdb_ids, seq_struct in seq_structs.iteritems():
            for match in regex.finditer(struct_patt, seq_struct[1],
                                        overlapped=True):
                matches.append([seq_struct[0][slice(*(match.span()))],
                                seq_struct[1][slice(*(match.span()))],
                                pdb_ids,
                                match.span()])

        seqs.extend(random.sample(matches, min(len(matches),
                                               sample_size - len(seqs))))

    return seqs


def _plot(values, plot_filename, plot_format, title, max_value=None):
    '''Plots 3d matrix values.'''
    fig = pylab.figure()
    sub_plot = fig.add_subplot(111)

    cmap = LinearSegmentedColormap.from_list(name='name', colors=['g', 'w'],
                                             N=10)

    cax = sub_plot.imshow(values, interpolation='nearest', cmap=cmap)
    cax.set_clim(0.0, max_value)
    sub_plot.set_title(title)

    # Add colorbar, make sure to specify tick locations to match desired tick
    # labels
    min_val = numpy.min(values)
    max_val = numpy.max(values)
    cbar = fig.colorbar(cax, ticks=[min_val, max_val])
    cbar.set_ticks([min_val, max_val])

    pylab.savefig(plot_filename, format=plot_format)
