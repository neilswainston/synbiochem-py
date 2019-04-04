'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=broad-except
# pylint: disable=no-member
# pylint: disable=protected-access
# pylint: disable=redefined-builtin
# pylint: disable=relative-import
# pylint: disable=superfluous-parens
# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals
# pylint: disable=ungrouped-imports
# pylint: disable=wrong-import-order
from collections import defaultdict
import itertools
import operator
import os
import random
import re
import ssl
from subprocess import call
import tempfile

from Bio import Seq, SeqIO, SeqRecord
from Bio.Blast import NCBIXML
from Bio.Data import CodonTable
from Bio.Restriction import Restriction, Restriction_Dictionary
from Bio.SeqUtils.MeltingTemp import Tm_NN
from synbiochem.biochem4j import taxonomy
from synbiochem.utils import thread_utils

from six.moves.urllib import parse
from six.moves.urllib import request


try:
    # Python 2:
    from requests.exceptions import ConnectionError
except ImportError:
    pass


NUCLEOTIDES = ['A', 'C', 'G', 'T']

AA_CODES = {'Ala': 'A',
            'Cys': 'C',
            'Asp': 'D',
            'Glu': 'E',
            'Phe': 'F',
            'Gly': 'G',
            'His': 'H',
            'Ile': 'I',
            'Lys': 'K',
            'Leu': 'L',
            'Met': 'M',
            'Asn': 'N',
            'Pro': 'P',
            'Gln': 'Q',
            'Arg': 'R',
            'Ser': 'S',
            'Thr': 'T',
            'Val': 'V',
            'Trp': 'W',
            'Tyr': 'Y',
            'Stop': '*'}

CODONS = {'A': [['G', 'C', 'ACGT']],
          'C': [['T', 'G', 'CT']],
          'D': [['G', 'A', 'CT']],
          'E': [['G', 'A', 'AG']],
          'F': [['T', 'T', 'CT']],
          'G': [['G', 'G', 'ACGT']],
          'H': [['C', 'A', 'CT']],
          'I': [['A', 'T', 'ACT']],
          'K': [['A', 'A', 'AG']],
          'L': [['C', 'T', 'ACGT'], ['T', 'T', 'AG']],
          'M': [['A', 'T', 'G']],
          'N': [['A', 'A', 'CT']],
          'P': [['C', 'C', 'ACGT']],
          'Q': [['C', 'A', 'AG']],
          'R': [['C', 'G', 'ACGT'], ['A', 'G', 'AG']],
          'S': [['T', 'C', 'ACGT'], ['A', 'G', 'CT']],
          'T': [['A', 'C', 'ACGT']],
          'V': [['G', 'T', 'ACGT']],
          'W': [['T', 'G', 'G']],
          'Y': [['T', 'A', 'CT']],
          'Stop': [['T', 'A', 'AG'], ['T', 'G', 'A']]}

NUCL_CODES = {
    'A': 'A',
    'C': 'C',
    'G': 'G',
    'T': 'T',
    'AG': 'R',
    'CT': 'Y',
    'CG': 'S',
    'AT': 'W',
    'GT': 'K',
    'AC': 'M',
    'CGT': 'B',
    'AGT': 'D',
    'ACT': 'H',
    'ACG': 'V',
    'ACGT': 'N',
}

INV_NUCL_CODES = {val: key for key, val in NUCL_CODES.items()}

NA = 'NA'
K = 'K'
TRIS = 'TRIS'
MG = 'MG'
DNTP = 'DNTP'

__DEFAULT_REAG_CONC = {NA: 0.05, K: 0, TRIS: 0, MG: 0.01, DNTP: 0}

AA_COD = defaultdict(list)

for cod, am_ac in \
        CodonTable.unambiguous_dna_by_name['Standard'].forward_table.items():
    AA_COD[am_ac].append(cod)


def get_codon_usage_organisms(expand=False, verbose=False):
    '''Gets name to taxonomy id dictionary of available codon usage tables.'''
    destination = os.path.dirname(os.path.realpath(__file__))
    filename = 'expand.txt' if expand else 'normal.txt'
    filepath = os.path.join(destination, filename)

    if not os.path.exists(filepath):
        # Download:
        if not os.path.exists(destination):
            os.makedirs(destination)

        url = 'ftp://ftp.kazusa.or.jp/pub/codon/current/species.table'
        tmp = tempfile.NamedTemporaryFile(delete=False)
        request.urlretrieve(url, tmp.name)

        # Read:
        codon_orgs = _read_codon_usage_orgs_file(tmp.name)

        # Expand:
        if expand:
            _expand_codon_usage_orgs(codon_orgs, verbose)

        # Save:
        _write_codon_usage_orgs_file(codon_orgs, filepath)

        return codon_orgs

    return _read_codon_usage_orgs_file(filepath)


class CodonOptimiser():
    '''Class to support codon optimisation.'''

    def __init__(self, taxonomy_id):
        self.__taxonomy_id = taxonomy_id
        self.__aa_to_codon_prob = self.__get_codon_usage()
        self.__codon_prob = {item[0]: item[1]
                             for lst in self.__aa_to_codon_prob.values()
                             for item in lst}

        self.__codon_to_w = {}

        for key in self.__aa_to_codon_prob:
            aa_dict = {a: b / self.__aa_to_codon_prob[key][0][1]
                       for a, b in self.__aa_to_codon_prob[key]}
            self.__codon_to_w.update(aa_dict)

    def get_codon_prob(self, codon):
        '''Gets the codon probability.'''
        return self.__codon_prob[codon]

    def get_codon_optim_seq(self, protein_seq, excl_codons=None,
                            max_repeat_nuc=float('inf'), restr_enzyms=None,
                            max_attempts=1000, tolerant=False, stepback=3):
        '''Returns a codon optimised DNA sequence.'''
        if max_repeat_nuc == float('inf') and restr_enzyms is None:
            return ''.join([self.get_random_codon(aa, excl_codons)
                            for aa in protein_seq])

        attempts = 0
        seq = ''
        i = 0
        blockage_i = -1
        inv_patterns = 0

        while attempts < max_attempts:
            amino_acid = protein_seq[i]
            new_seq = seq + self.get_random_codon(amino_acid, excl_codons)

            invalids = find_invalid(new_seq, max_repeat_nuc,
                                    restr_enzyms)

            if len(invalids) == inv_patterns or \
                    (attempts == max_attempts - 1 and tolerant):

                if i == blockage_i:
                    if attempts == max_attempts - 1:
                        inv_patterns = inv_patterns + 1

                    attempts = 0

                seq = new_seq

                if i == len(protein_seq) - 1:
                    return seq

                i += 1
            else:
                blockage_i = max(i, blockage_i)
                i = max(0, (invalids[-1] // 3) - stepback)
                seq = seq[:i * 3]
                attempts += 1

        raise ValueError('Unable to generate codon-optimised sequence.')

    def get_cai(self, dna_seq):
        '''Gets the CAI for a given DNA sequence.'''
        cai = 0

        for i in range(0, len(dna_seq), 3):
            cai += self.__codon_to_w[dna_seq[i:i + 3]]

        return cai / (len(dna_seq) / 3)

    def mutate(self, protein_seq, dna_seq, mutation_rate):
        '''Mutate a protein-encoding DNA sequence according to a
        supplied mutation rate.'''
        return ''.join([self.get_random_codon(amino_acid)
                        if random.random() < mutation_rate
                        else dna_seq[3 * i:3 * (i + 1)]
                        for i, amino_acid in enumerate(protein_seq)])

    def get_all_codons(self, amino_acid):
        '''Returns all codons for a given amino acid.'''
        return [t[0] for t in self.__aa_to_codon_prob[amino_acid]]

    def get_best_codon(self, amino_acid):
        '''Get 'best' codon for a given amino acid.'''
        return self.__aa_to_codon_prob[amino_acid][0][0]

    def get_random_codon(self, amino_acid, excl_codons=None):
        '''Returns a random codon for a given amino acid,
        based on codon probability from the codon usage table.'''
        if excl_codons is None:
            excl_codons = []

        codon_usage = [codon_usage
                       for codon_usage in self.__aa_to_codon_prob[amino_acid]
                       if codon_usage[0] not in excl_codons]

        if not codon_usage:
            raise ValueError('No codons available for ' + amino_acid +
                             ' after excluding ' + str(excl_codons))

        while True:
            rand = random.random()
            cumulative_prob = 0

            for codon, prob in iter(reversed(codon_usage)):
                cumulative_prob += prob

                if cumulative_prob > rand:
                    return codon

    def __get_codon_usage(self):
        '''Gets the codon usage table for a given taxonomy id.'''
        aa_to_codon_prob = {aa_code: {} for aa_code in AA_CODES.values()}

        url = 'http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=' \
            + self.__taxonomy_id + '&aa=1&style=GCG'

        in_codons = False

        for line in request.urlopen(url):
            line = line.decode('utf-8')

            if line == '<PRE>\n':
                in_codons = True
            elif line == '</PRE>\n':
                break
            elif in_codons:
                values = re.split('\\s+', line)
                am_acid = 'Stop' if values[0] == 'End' else values[0]

                if am_acid in AA_CODES:
                    codon_prob = aa_to_codon_prob[AA_CODES[am_acid]]
                    codon_prob[values[1]] = float(values[3])

        aa_to_codon_prob.update((x, _scale(y))
                                for x, y in aa_to_codon_prob.items())

        return aa_to_codon_prob


def find_invalid(seq, max_repeat_nuc=float('inf'), restr_enzyms=None):
    '''Finds invalid sequences.'''
    inv = []
    seq = seq.upper()

    # Invalid repeating nucleotides:
    if max_repeat_nuc != float('inf'):
        pattern = [''.join([nucl] * (max_repeat_nuc + 1))
                   for nucl in NUCLEOTIDES]
        pattern = re.compile(r'(?=(' + '|'.join(pattern) + '))')

        inv = [m.start() for m in pattern.finditer(seq)]

    # Invalid restriction sites:
    if restr_enzyms:
        for rest_enz in [_get_restr_type(name) for name in restr_enzyms]:
            inv.extend(rest_enz.search(Seq.Seq(seq)))

    return inv


def is_invalid(seq, max_repeat_nuc=float('inf'), restr_enzyms=None):
    '''Checks whether a sequence is valid.'''
    return len(find_invalid(seq, max_repeat_nuc, restr_enzyms)) > 0


def get_all_rev_trans(aa_seq):
    '''Returns all reverse translations of amino acid sequence.'''
    codons = [AA_COD[aa] for aa in aa_seq.strip()]
    return [''.join(t) for t in list(itertools.product(*codons))]


def get_random_dna(length, max_repeat_nuc=float('inf'), restr_enzyms=None):
    '''Returns a random sequence of DNA of the supplied length,
    while adhering to a maximum number of repeating nucleotides.'''
    max_attempts = 100
    attempts = 0
    len_add = 16

    seq = ''

    while True:
        attempts += 1

        if attempts > max_attempts:
            raise ValueError('Unable to optimise sequence.')

        while len(seq) < length:
            seq += _get_random_dna(len_add)

            if is_invalid(seq, max_repeat_nuc, restr_enzyms):
                seq = seq[:-len_add]

        if not is_invalid(seq, max_repeat_nuc, restr_enzyms):
            return seq[:length]

    return None


def mutate_seq(seq, mutations=1, alphabet=None):
    '''Mutates sequence.'''
    if alphabet is None:
        alphabet = NUCLEOTIDES

    seq_new = seq

    for _ in range(mutations):
        move = random.random()
        pos = int(random.random() * len(seq))
        base = random.choice(alphabet)

        # Insert:
        if move < 0.1:
            seq_new = seq_new[1:pos + 1] + base + seq_new[pos + 1:]

        # Delete:
        elif move < 0.2:
            seq_new = base + seq_new[:pos] + seq_new[pos + 1:]

        # Replace:
        else:
            seq_new = seq_new[:pos] + base + seq_new[pos + 1:]

    return seq_new


def get_melting_temp(dna1, dna2=None, reag_concs=None, strict=True):
    '''Calculates melting temperarure of DNA sequence against its
    complement, or against second DNA sequence using Nearest-Neighbour
    method.'''
    assert len(dna1) > 1

    reagent_concs = __DEFAULT_REAG_CONC

    if reag_concs is not None:
        reagent_concs.update(reag_concs)

    reagent_conc = {k: v * 1000 for k, v in reagent_concs.items()}
    dnac1 = 30

    return Tm_NN(dna1, check=True, strict=strict, c_seq=dna2, shift=0,
                 Na=reagent_conc[NA], K=reagent_conc[K],
                 Tris=reagent_conc[TRIS], Mg=reagent_conc[MG],
                 dNTPs=reagent_conc[DNTP],
                 dnac1=dnac1, dnac2=dnac1, selfcomp=dna2 is None,
                 saltcorr=7)


def get_seq_by_melt_temp(seq, target_melt_temp, forward=True,
                         terminii=None,
                         reagent_concs=None,
                         tol=0.025):
    '''Returns a subsequence close to desired melting temperature.'''
    if terminii is None:
        terminii = ['A', 'C', 'G', 'T']
    else:
        terminii = [term.upper() for term in terminii]

    best_delta_tm = float('inf')
    best_subseq = ''
    best_melt_temp = float('NaN')
    in_tol = False

    for i in range(1, len(seq)):
        subseq = seq[:(i + 1)] if forward else seq[-(i + 1):]
        melt_temp = get_melting_temp(subseq, None, reagent_concs)

        if subseq[-1 if forward else 0].upper() in terminii:
            delta_tm = abs(melt_temp - target_melt_temp)

            if delta_tm / target_melt_temp < tol:
                in_tol = True

                if delta_tm < best_delta_tm:
                    best_delta_tm = delta_tm
                    best_subseq = subseq
                    best_melt_temp = melt_temp
            elif in_tol:
                return best_subseq, best_melt_temp

    raise ValueError('Unable to get sequence of required melting temperature')


def get_rand_seq_by_melt_temp(target_melt_temp,
                              max_repeat_nuc=float('inf'),
                              restr_enzyms=None,
                              reagent_concs=None,
                              tol=0.025):
    '''Returns a random close to desired melting temperature.'''
    seq = random.choice(NUCLEOTIDES)

    while True:
        seq += random.choice(NUCLEOTIDES)

        if is_invalid(seq, max_repeat_nuc, restr_enzyms):
            seq = random.choice(NUCLEOTIDES)
            continue

        melt_temp = get_melting_temp(seq, None, reagent_concs)

        delta_tm = abs(melt_temp - target_melt_temp)

        if delta_tm / target_melt_temp < tol:
            return seq, melt_temp

    raise ValueError('Unable to get sequence of required melting temperature')


def get_uniprot_values(uniprot_ids, fields, batch_size=128, verbose=False,
                       num_threads=0):
    '''Gets dictionary of ids to values from Uniprot.'''
    values = []

    if num_threads:
        thread_pool = thread_utils.ThreadPool(num_threads)

        for i in range(0, len(uniprot_ids), batch_size):
            thread_pool.add_task(_get_uniprot_batch, uniprot_ids, i,
                                 batch_size, fields, values, verbose)

        thread_pool.wait_completion()
    else:
        for i in range(0, len(uniprot_ids), batch_size):
            _get_uniprot_batch(uniprot_ids, i, batch_size, fields, values,
                               verbose)

    return {value['Entry']: value for value in values}


def search_uniprot(query, fields, limit=128):
    '''Gets dictionary of ids to values from Uniprot.'''
    values = []

    url = 'http://www.uniprot.org/uniprot/?query=' + parse.quote(query) + \
        '&sort=score&limit=' + str(limit) + \
        '&format=tab&columns=id,' + ','.join([parse.quote(field)
                                              for field in fields])

    _parse_uniprot_data(url, values)

    return values


def do_blast(id_seqs_subjects, id_seqs_queries, program='blastn',
             dbtype='nucl', evalue=1.0, word_size=28):
    '''Performs BLAST of query sequences against subject sequences.'''
    db_filename = write_fasta(id_seqs_subjects)
    query_filename = write_fasta(id_seqs_queries)
    result_file = tempfile.NamedTemporaryFile(prefix='blast_result_',
                                              suffix='.xml',
                                              delete=False)
    log_file = tempfile.NamedTemporaryFile(prefix='makeblastdb_log',
                                           suffix='.txt',
                                           delete=False)

    call(['makeblastdb',
          '-in', db_filename,
          '-out', db_filename,
          '-dbtype', dbtype,
          '-logfile', log_file.name])

    call([program,
          '-query', query_filename,
          '-db', db_filename,
          '-out', result_file.name,
          '-evalue', str(evalue),
          '-word_size', str(word_size),
          '-outfmt', '5'])

    return NCBIXML.parse(open(result_file.name))


def do_clustal(in_data, is_fasta_file=False, result_file=None,
               guidetree_file=None):
    '''Performs Clustal Omega multiple sequence alignment.'''
    result_file = tempfile.NamedTemporaryFile(prefix='clustalo_result_',
                                              suffix='.fasta',
                                              delete=False).name \
        if result_file is None \
        else result_file

    guidetree_file = tempfile.NamedTemporaryFile(prefix='clustalo_tree_',
                                                 suffix='.dnd',
                                                 delete=False).name \
        if guidetree_file is None \
        else guidetree_file

    call(['clustalo',
          '-i', in_data if is_fasta_file else write_fasta(in_data),
          '-o', result_file,
          '--guidetree-out=' + guidetree_file,
          '--force'])

    return read_fasta(result_file)


def read_fasta(filename):
    '''Reads a fasta file.'''
    with open(filename, 'rU') as fle:
        seqs = {record.id: str(record.seq)
                for record in SeqIO.parse(fle, 'fasta')}

    return seqs


def write_fasta(id_seqs, filename=None):
    '''Writes a fasta file.'''
    if filename is None:
        temp_file = tempfile.NamedTemporaryFile(prefix='fasta_', suffix='.txt',
                                                delete=False)
        filename = temp_file.name

    records = [SeqRecord.SeqRecord(Seq.Seq(seq), str(seq_id), '', '')
               for seq_id, seq in id_seqs.items()]

    SeqIO.write(records, filename, 'fasta')

    return filename


def pcr(seq, forward_primer, reverse_primer):
    '''Apply in silico PCR.'''
    for_primer_pos = seq.find(forward_primer.upper())

    rev_primer_pos = \
        seq.find(str(Seq.Seq(reverse_primer).reverse_complement().upper()))

    if for_primer_pos > -1 and rev_primer_pos > -1:
        seq = seq[for_primer_pos:] + \
            seq[:rev_primer_pos + len(reverse_primer)]
    elif for_primer_pos > -1:
        seq = seq[for_primer_pos:]
    elif rev_primer_pos > -1:
        seq = seq[:rev_primer_pos + len(reverse_primer)]

    return seq, for_primer_pos


def _scale(codon_usage):
    '''Scale codon usage values to add to 1.'''
    codon_usage = {key: value / sum(codon_usage.values())
                   for key, value in codon_usage.items()}

    return sorted(codon_usage.items(), key=operator.itemgetter(1),
                  reverse=True)


def _get_random_dna(length):
    '''Returns a random sequence of DNA of the supplied length.'''
    return ''.join(random.choice(['A', 'C', 'G', 'T']) for _ in range(length))


def _get_uniprot_batch(uniprot_ids, i, batch_size, fields, values, verbose):
    '''Get batch of Uniprot data.'''
    if verbose:
        print('seq_utils: getting Uniprot values ' + str(i) + ' - ' +
              str(min(i + batch_size, len(uniprot_ids))) + ' / ' +
              str(len(uniprot_ids)))

    batch = uniprot_ids[i:min(i + batch_size, len(uniprot_ids))]
    query = '+or+'.join(['id:' + uniprot_id for uniprot_id in batch])
    url = 'https://www.uniprot.org/uniprot/?query=' + query + \
        '&format=tab&columns=id,' + ','.join([parse.quote(field)
                                              for field in fields])

    _parse_uniprot_data(url, values)


def _parse_uniprot_data(url, values):
    '''Parses Uniprot data.'''
    headers = None

    try:
        context = ssl._create_unverified_context()

        for line in request.urlopen(url, context=context):
            line = line.decode('utf-8')
            tokens = line.strip().split('\t')

            if headers is None:
                headers = tokens
            else:
                resp = dict(zip(headers, tokens))

                if 'Protein names' in resp:
                    regexp = re.compile(r'(?<=\()[^)]*(?=\))|^[^(][^()]*')
                    names = regexp.findall(resp.pop('Protein names'))
                    resp['Protein names'] = [nme.strip() for nme in names]

                for key in resp:
                    if key.startswith('Cross-reference'):
                        resp[key] = resp[key].split(';')

                values.append(resp)
    except Exception as err:
        print(err)


def _read_codon_usage_orgs_file(filename):
    '''Reads Codon Usage Database table of species file.'''
    codon_orgs = {}

    with open(filename, 'r') as textfile:
        next(textfile)

        for line in textfile:
            tokens = line.strip().split('\t')
            codon_orgs[tokens[0]] = tokens[1]

    return codon_orgs


def _expand_codon_usage_orgs(codon_orgs, verbose, max_errors=16):
    '''Expand Codon Usage Db table of species with children and synonyms.'''
    for tax_id in codon_orgs.values():

        if verbose:
            print('Expanding codon usage for NCBI Taxonomy id: ' + tax_id)

        errors = 0
        success = False

        while not success:
            try:
                for name in taxonomy.get_synonyms_by_id(tax_id):
                    _add_codon_usage_org(codon_orgs, name, tax_id)

                for child in taxonomy.get_children_by_id(tax_id):
                    _add_codon_usage_org(codon_orgs, child['name'], tax_id)

                    for name in child['names']:
                        _add_codon_usage_org(codon_orgs, name, tax_id)

                success = True

            except ConnectionError as err:
                errors += 1

                if errors == max_errors:
                    raise err


def _add_codon_usage_org(codon_orgs, name, tax_id):
    '''Adds name to codon_orgs.'''
    if name not in codon_orgs:
        codon_orgs[name] = tax_id


def _write_codon_usage_orgs_file(codon_orgs, filepath):
    '''Writes Codon Usage Database table of species file.'''
    with open(filepath, 'w+') as fle:
        fle.write('Name\tId\n')

        for name, tax_id in codon_orgs.items():
            fle.write(name + '\t' + tax_id + '\n')


def _get_restr_type(name):
    '''Gets RestrictionType from name.'''
    types = [
        x for _, (x, y) in Restriction_Dictionary.typedict.items()
        if name in y][0]

    enz_types = tuple(getattr(Restriction, typ)
                      for typ in types)

    return Restriction.RestrictionType(
        str(name), enz_types, Restriction_Dictionary.rest_dict[name])
