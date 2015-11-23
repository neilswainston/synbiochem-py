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
import tempfile
import urllib

from Bio import SeqUtils
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder

import synbiochem.utils.io_utils as io_utils

_DIR = 'structure_utils'


def get_pdb_ids(max_ids=None):
    '''Returns all PDB ids.'''
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
                # Do something:
                str_data += line[:-1]

    return seq_struct


def get_sequences(pdb_id):
    '''Gets the sequences in a PDB file.'''
    return [SeqUtils.seq1(''.join([residue.get_resname()
                                   for residue in chn
                                   if 'CA' in residue.child_dict]))
            for chn in get_structure(pdb_id).get_chains()]


def get_structure(pdb_id):
    '''Returns a PDB structure.'''
    with tempfile.TemporaryFile() as pdb_file:
        opener = urllib.URLopener()
        opener.retrieve('http://www.rcsb.org/pdb/files/' + pdb_id + '.pdb',
                        pdb_file.name)
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


def get_phi_psi_data(pdb_id):
    '''Gets phi and phi angle data.'''
    builder = PPBuilder()
    return [polypep.get_phi_psi_list()
            for polypep in builder.build_peptides(get_structure(pdb_id))]


def plot_proximities(pdb_id):
    '''Plots proximity plot(s).'''
    all_proximities = calc_proximities(pdb_id)

    plot_format = 'png'

    for idx, proximities in enumerate(all_proximities):
        name = pdb_id + '_' + str(idx + 1)
        _plot(proximities, name + '.' + plot_format, plot_format,
              name + ' proximity plot')


def sample_seqs(sample_size, struct_patt):
    '''Sample sequence and structure data.'''
    seqs = []
    patt = re.compile(struct_patt)

    while len(seqs) < sample_size:
        pdb_ids = get_pdb_ids(sample_size)
        seq_struct = get_seq_struct(pdb_ids)
        try:
            assert all([len(v[0]) == len(v[1]) for v in seq_struct.values()])
        except TypeError:
            print 'error'
        matches = set([v[0][slice(*(m.span()))]
                       for v in seq_struct.values()
                       for m in patt.finditer(v[1])])

        seqs.extend(random.sample(matches,
                                  min(len(matches), sample_size - len(seqs))))

    return seqs[:sample_size]


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
