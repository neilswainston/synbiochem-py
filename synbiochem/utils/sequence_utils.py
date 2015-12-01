'''
synbiochemdev (c) University of Manchester 2015

synbiochemdev is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import csv
import itertools
import math
import operator
import os
import random
import re
import subprocess
import urllib2
import uuid


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
            'End': '*'}


NA = 'NA'
K = 'K'
TRIS = 'TRIS'
MG = 'MG'
DNTP = 'DNTP'


class CodonOptimiser(object):
    '''Class to support codon optimisation.'''

    def __init__(self, taxonomy_id):
        self.__taxonomy_id = taxonomy_id
        self.__codon_usage_table = self.__get_codon_usage_table()
        self.__codon_to_w = {}

        for key in self.__codon_usage_table:
            aa_dict = dict([(a, b / self.__codon_usage_table[key][0][1])
                            for a, b in self.__codon_usage_table[key]])
            self.__codon_to_w.update(aa_dict)

    def optimise(self, protein_seqs, max_repeat_nuc=float('inf')):
        '''Codon optimises the supplied protein sequences.'''
        optimised_seqs = []
        max_attempts = 1000
        attempts = 0

        for protein_seq in protein_seqs:
            while True:
                attempts += 1

                if attempts > max_attempts:
                    raise ValueError('Unable to optimise sequence. ' +
                                     'Greater than ' + str(max_repeat_nuc) +
                                     ' repeating nucleotides.')

                optimised_seq = self.get_codon_optimised_seq(protein_seq)

                if is_valid(optimised_seq, max_repeat_nuc):
                    optimised_seqs.append(optimised_seq)
                    break

        return optimised_seqs

    def get_codon_optimised_seq(self, protein_seq, invalid_pattern=None,
                                max_attempts=1000):
        '''Returns a codon optimised DNA sequence.
        TODO: *May* get stuck in a loop. Check max iteration. e.g. WG protein
        seq guarantees 4 consecutive nucleotides.'''
        if invalid_pattern is None:
            return ''.join([self.get_random_codon(aa)
                            for aa in protein_seq])
        else:
            attempts = 0
            seq = ''
            i = 0

            while attempts < max_attempts:
                amino_acid = protein_seq[i]
                new_seq = seq + self.get_random_codon(amino_acid)

                if count_pattern(new_seq, invalid_pattern) > 0:
                    num_problem_codons = len(max(re.findall(invalid_pattern,
                                                            new_seq),
                                                 key=len)) / 3
                    i -= num_problem_codons
                    seq = seq[:-num_problem_codons * 3]
                    attempts += 1
                else:
                    attempts = 0
                    seq = new_seq

                    if i == len(protein_seq) - 1:
                        return seq

                    i += 1

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

    def get_all_rev_trans(self, aa_seq):
        '''Returns all reverse translations of amino acid sequence.'''
        codons = [self.get_all_codons(aa) for aa in aa_seq]
        return [''.join(t) for t in list(itertools.product(*codons))]

    def get_all_codons(self, amino_acid):
        '''Returns all codons for a given amino acid.'''
        return [t[0] for t in self.__codon_usage_table[amino_acid]]

    def get_random_codon(self, amino_acid):
        '''Returns a random codon for a given amino acid,
        based on codon probability from the codon usage table.'''
        codon_usage = self.__codon_usage_table[amino_acid]
        rand = random.random()
        cumulative_prob = 0

        for codon, prob in iter(reversed(codon_usage)):
            cumulative_prob += prob

            if cumulative_prob > rand:
                return codon

    def __get_codon_usage_table(self):
        '''Gets the codon usage table for a given taxonomy id.'''
        codon_usage_table = {aa_code: {} for aa_code in AA_CODES.values()}

        url = 'http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=' \
            + self.__taxonomy_id + '&aa=1&style=GCG'

        in_codons = False

        for line in urllib2.urlopen(url):
            if line == '<PRE>\n':
                in_codons = True
            elif line == '</PRE>\n':
                break
            elif in_codons:
                values = re.split('\\s+', line)

                if values[0] in AA_CODES:
                    codon_usage = codon_usage_table[AA_CODES[values[0]]]
                    codon_usage[values[1]] = float(values[3])

        codon_usage_table.update((x, _scale(y))
                                 for x, y in codon_usage_table.items())

        return codon_usage_table


def get_minimum_free_energy(exec_dir, sequences):
    '''Returns minimum free energy of supplied DNA / RNA sequences.'''
    with open(_get_temp_filename(), 'w') as input_file, \
            open(_get_temp_filename(), 'w') as output_file:

        for i, sequence in enumerate(sequences):
            input_file.write('>Seq' + str(i) + '\n' + sequence + '\n')
            input_file.flush()

        input_file.close()

        proc = subprocess.Popen(os.path.join(exec_dir, 'RNAfold'),
                                stdin=open(input_file.name),
                                stdout=output_file)

        proc.wait()

        _cleanup(os.getcwd(), 'Seq\\d+_ss.ps')

        mfes = []
        pattern = re.compile(r'[+-]?\d+\.\d+')

        with open(output_file.name) as out_file:
            for line in out_file.readlines():
                src = pattern.search(line)

                if src:
                    mfes.append(float(src.group()))

        _cleanup(os.getcwd(), input_file.name)
        _cleanup(os.getcwd(), output_file.name)

        return mfes


def get_random_dna(length, max_repeat_nuc=float('inf')):
    '''Returns a random sequence of DNA of the supplied length,
    while adhering to a maximum number of repeating nucleotides.'''
    max_attempts = 1000
    attempts = 0

    while True:
        attempts += 1

        if attempts > max_attempts:
            raise ValueError('Unable to optimise sequence. ' +
                             'Greater than ' + str(max_repeat_nuc) +
                             ' repeating nucleotides.')

        random_dna = _get_random_dna(length)

        if is_valid(random_dna, max_repeat_nuc):
            return random_dna

    return None


def get_melting_temp(dna1, dna2=None, reagent_concs=None):
    '''Calculates melting temperarure of DNA sequence against its
    complement, or against second DNA sequence.

    See von Ahsen N, Wittwer CT, Schutz E: Oligonucleotide melting
    temperatures under PCR conditions: deoxynucleotide Triphosphate and
    Dimethyl sulfoxide concentrations with comparison to alternative
    empirical formulas. Clin Chem 2001, 47:1956-1961.'''

    base_melt_temp = _get_base_melting_temp(reagent_concs)

    if dna2 is None:
        dna_gc_content = len(re.findall('G|C', dna1.upper()))
        return _get_melting_temp(dna_gc_content, 0, len(dna1), base_melt_temp)
    else:
        # Not implemented...
        return -1


def is_valid(dna_seq, max_repeat_nuc):
    '''Checks whether a DNA sequence is valid, in terms of a supplied maximum
    number of repeating nucleotides.'''
    nuc_count = 0
    prev_nuc = ''

    for nuc in dna_seq:
        if prev_nuc == nuc:
            nuc_count += 1
        else:
            prev_nuc = nuc
            nuc_count = 1

        if nuc_count > max_repeat_nuc:
            return False

    return True


def get_sequences(protein_ids):
    '''Returns sequences from protein ids, which may be either Uniprot ids,
    or a protein sequence itself.'''
    uniprot_id_pattern = \
        '[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}'

    sequences = {}

    for idx, protein_id in enumerate(protein_ids):
        if re.match(uniprot_id_pattern, protein_id):
            sequences.update(get_uniprot_seqs([protein_id]))
        else:
            sequences[str(idx)] = protein_id

    return sequences


def get_uniprot_seqs(uniprot_ids):
    '''Gets dictionary of ids to sequences from Uniprot.'''
    query = '+or+'.join(['id:' + uniprot_id for uniprot_id in uniprot_ids])
    url = 'http://www.uniprot.org/uniprot/?query=' + query + \
        '&format=tab&columns=id,sequence'
    return dict(list(csv.reader(urllib2.urlopen(url), delimiter='\t'))[1:])


def count_pattern(strings, pattern):
    '''Counts pattern in string of list of strings.'''
    if isinstance(strings, str) or isinstance(strings, unicode):
        return len(re.findall(pattern, strings))
    elif strings is None:
        return 0
    else:
        return [count_pattern(s, pattern) for s in strings]


def _scale(codon_usage):
    '''Scale codon usage values to add to 1.'''
    codon_usage = dict([(key, value / sum(codon_usage.values()))
                        for key, value in codon_usage.items()])

    return sorted(codon_usage.items(), key=operator.itemgetter(1),
                  reverse=True)


def _get_random_dna(length):
    '''Returns a random sequence of DNA of the supplied length.'''
    return ''.join(random.choice(['A', 'C', 'G', 'T']) for _ in range(length))


def _get_base_melting_temp(reagent_concs=None):
    '''Gets the base melting temperature based on reagent concentrations.'''
    na_equivalence_parameter = 3.79

    if reagent_concs is None:
        reagent_concs = {NA: 0.05, K: 0, TRIS: 0, MG: 0.002, DNTP: 1e-7}

    na_equiv = reagent_concs[NA] + reagent_concs[K] \
        + reagent_concs[TRIS] / 2 \
        + na_equivalence_parameter \
        * math.sqrt(reagent_concs[MG] - reagent_concs[DNTP])

    # See Wetmur J: DNA probes: applications of the principles of nucleic
    # acid hybridization. Crit Rev Biochem Mol Biol 1991, 26:227-259.
    return 81.5 + 16.6 * math.log10(na_equiv / (1.0 + 0.7 * na_equiv))


def _get_melting_temp(dna_gc_content, mismatches, length, base_melting_temp):
    '''Calculates melting temperature.'''
    gc_increment = 41
    mismatch_decrement = 100
    fixed_decrement = 500
    return base_melting_temp + \
        ((gc_increment * dna_gc_content) - fixed_decrement -
         (mismatch_decrement * mismatches)) / float(length)


def _cleanup(drctry, pattern):
    '''Deletes files in directory matching pattern.'''
    for filename in os.listdir(drctry):
        if re.search(pattern, filename):
            os.remove(os.path.join(drctry, filename))


def _get_temp_filename():
    '''Gets temporary file name.'''
    return str(uuid.uuid1()) + '.txt'
