'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=too-many-arguments
from matplotlib.colors import LinearSegmentedColormap
import collections
import itertools
import os
import random
import re
import sys

from Bio import SeqUtils
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
import numpy
import pylab

import scipy.spatial.distance as dist
import synbiochem.utils.io_utils as io_utils


_DIR = 'structure_utils'
_PDB_DIR = 'pdb'
_SEC_STRUCT = 'EBHIGTSL '


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


def get_seq_structs(pdb_ids=None):
    '''Returns sequence and structure.'''
    seq_structs = {}
    pdb_ids = sorted(pdb_ids) if pdb_ids is not None else None
    in_field = False
    tokens = None
    str_data = ''

    source_url = 'http://www.rcsb.org/pdb/files/ss.txt'
    target_filename = os.path.join(os.path.expanduser('~'), _DIR,
                                   'ss.txt')

    with open(io_utils.get_file(source_url, target_filename)) as fle:
        for line in fle:
            if line.startswith('>'):
                pdb_id = re.search('(?<=\\>)[^:]*', line).group(0)

                if pdb_ids is None or pdb_id in pdb_ids:
                    if in_field:
                        if tokens[:2] not in seq_structs:
                            seq_structs[tokens[:2]] = [None, None]

                        seq_structs[tokens[:2]][0 if tokens[2] == 'sequence'
                                                else 1] = str_data
                        str_data = ''

                    tokens = tuple(re.split('>|:', line.strip())[1:])
                    in_field = True

                elif in_field:
                    if tokens[:2] not in seq_structs:
                        seq_structs[tokens[:2]] = [None, None]

                    seq_structs[tokens[:2]][0 if tokens[2] == 'sequence'
                                            else 1] = str_data
                    str_data = ''

                    in_field = False
                    tokens = None
                    str_data = ''

            elif in_field:
                str_data += line[:-1]

    return {key: value for key, value in seq_structs.iteritems()
            if all(val is not None for val in value)}


def get_pep_structs(struct_set, length):
    '''Returns sequence and structure.'''
    for struct in struct_set:
        assert struct in _SEC_STRUCT

    return [pep for struct in struct_set for pep in _get_peps(struct, length)]


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
    all_proximities = []
    structure = get_structure(pdb_id)
    chains = [c for c in structure.get_chains()]

    for chn in chains:
        res_coords = []

        for residue in chn:
            if 'CA' in residue.child_dict:
                res_coords.append([atom.get_coord()
                                   for atom_id, atom
                                   in residue.child_dict.iteritems()
                                   if atom_id
                                   not in ['C', 'N', 'O', 'OXT']])

        residues = [residue.get_id()[1]
                    for residue in chn
                    if 'CA' in residue.child_dict]

        chn_prox = collections.defaultdict(list)

        for res1, res2 in itertools.product(zip(residues, res_coords),
                                            repeat=2):
            chn_prox[res1[0]].append(numpy.amin(dist.cdist(res1[1], res2[1],
                                                           'euclidean')))

        all_proximities.append([chn_prox.keys(), chn_prox.values()])

    return all_proximities


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
        _plot(proximities[0], proximities[1], name + '.' + plot_format,
              plot_format, name + ' proximity plot')


def _plot(residues, values, plot_filename, plot_format, title, max_value=None):
    '''Plots 3d matrix values.'''
    fig = pylab.figure()
    sub_plot = fig.add_subplot(111)

    cmap = LinearSegmentedColormap.from_list(name='name', colors=['g', 'w'],
                                             N=10)

    cax = sub_plot.imshow(values, interpolation='nearest', cmap=cmap,
                          extent=[residues[0], residues[-1],
                                  residues[-1], residues[0]])
    cax.set_clim(0.0, max_value)
    sub_plot.set_title(title)

    # Add colorbar, make sure to specify tick locations to match desired tick
    # labels
    min_val = numpy.min(values)
    max_val = numpy.max(values)
    cbar = fig.colorbar(cax, ticks=[min_val, max_val])
    cbar.set_ticks([min_val, max_val])

    pylab.savefig(plot_filename, format=plot_format)


def _get_peps(struct, length):
    '''Gets n-mers of length, matching a given structure.'''
    directory = os.path.join(os.path.expanduser('~'), _DIR, str(length))
    filename = _get_struct_filename(directory, struct)

    if not os.path.isfile(filename):
        _write_peps(length, directory)

    with open(filename) as fle:
        return [line.split()[0] for line in fle]


def _write_peps(length, directory):
    '''Writes n-mers of length, matching a given structure, to file.'''
    assert length % 2 == 1

    pep_structs = {}

    for seq_struct in get_seq_structs().values():
        seq = seq_struct[0]
        struct = seq_struct[1]

        for i in range(len(seq) - length + 1):
            pep_seq = seq[i:i + length]
            pep_struct = struct[i:i + length]

            if pep_seq not in pep_structs:
                pep_structs[pep_seq] = set([pep_struct])
            else:
                pep_structs[pep_seq].add(pep_struct)

    classified_peps = _get_classified_peps(pep_structs, length / 2)
    _write_classified_peps(classified_peps, directory)


def _get_classified_peps(pep_structs, middle_index):
    '''Gets peptides classified by middle residue secondary structure,
    e.g. H, S, B, etc.'''
    classified_peps = {key: [] for key in list(_SEC_STRUCT) + ['unknown']}

    for pep, structs in pep_structs.iteritems():
        structs = list(structs)
        clss = 'unknown' if len(structs) > 1 else structs[0][middle_index]
        classified_peps[clss].append([pep, structs])

    return classified_peps


def _write_classified_peps(classified_peps, directory):
    '''Writes classified peptides to file.'''
    if not os.path.exists(directory):
        os.makedirs(directory)

    for clss, peps in classified_peps.iteritems():
        with open(_get_struct_filename(directory, clss), 'w+') as outfile:
            for pep in peps:
                outfile.write('\t'.join([str(val) for val in pep]) + '\n')


def _get_struct_filename(directory, struct):
    '''Gets filename for a given secondary structure class.'''
    # struct = struct if struct is not ' ' else 'blank'
    return os.path.join(directory, struct + '.xls')


def main(args):
    '''main method'''
    plot_proximities(args[0])

if __name__ == '__main__':
    main(sys.argv[1:])
