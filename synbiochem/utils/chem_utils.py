'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import itertools
import re

import synbiochem.utils.math_utils as math_utils


__ELEMENTAL_MASSES = {
    'H': 1.007825032,
    'He': 3.01602932,
    'Li': 6.015122887,
    'Be': 9.012183065,
    'B': 10.01293695,
    'C': 12,
    'N': 14.003074,
    'O': 15.99491462,
    'F': 18.99840316,
    'Ne': 19.99244018,
    'Na': 22.98976928,
    'Mg': 23.9850417,
    'Al': 26.98153853,
    'Si': 27.97692653,
    'P': 30.973762,
    'S': 31.97207117,
    'Cl': 34.96885268,
    'Ar': 35.96754511,
    'K': 38.96370649,
    'Ca': 39.96259086,
    'Sc': 44.95590828,
    'Ti': 45.95262772,
    'Cr': 49.94604183,
    'V': 49.94715601,
    'Fe': 53.93960899,
    'Mn': 54.93804391,
    'Ni': 57.93534241,
    'Co': 58.93319429,
    'Cu': 62.92959772,
    'Zn': 63.92914201,
    'Ga': 68.9255735,
    'Ge': 69.92424875,
    'Se': 73.92247593,
    'As': 74.92159457,
    'Kr': 77.92036494,
    'Br': 78.9183376,
    'Sr': 83.9134191,
    'Rb': 84.91178974,
    'Y': 88.9058403,
    'Zr': 89.9046977,
    'Mo': 91.90680796,
    'Nb': 92.906373,
    'Ru': 95.90759025,
    'Tc': 96.9063667,
    'Pd': 101.9056022,
    'Rh': 102.905498,
    'Cd': 105.9064599,
    'Ag': 106.9050916,
    'Sn': 111.9048239,
    'In': 112.9040618,
    'Te': 119.9040593,
    'Sb': 120.903812,
    'Xe': 123.905892,
    'I': 126.9044719,
    'Ba': 129.9063207,
    'Cs': 132.905452,
    'Ce': 135.9071292,
    'La': 137.9071149,
    'Pr': 140.9076576,
    'Nd': 141.907729,
    'Sm': 143.9120065,
    'Pm': 144.9127559,
    'Eu': 150.9198578,
    'Gd': 151.9197995,
    'Dy': 155.9242847,
    'Tb': 158.9253547,
    'Er': 161.9287884,
    'Ho': 164.9303288,
    'Yb': 167.9338896,
    'Tm': 168.9342179,
    'Hf': 173.9400461,
    'Lu': 174.9407752,
    'W': 179.9467108,
    'Ta': 179.9474648,
    'Os': 183.9524885,
    'Re': 184.9529545,
    'Pt': 189.9599297,
    'Ir': 190.9605893,
    'Hg': 195.9658326,
    'Au': 196.9665688,
    'Tl': 202.9723446,
    'Pb': 203.973044,
    'Bi': 208.9803991,
    'Po': 208.9824308,
    'At': 209.9871479,
    'Rn': 210.9906011,
    'Ra': 223.0185023,
    'Fr': 223.019736,
    'Ac': 227.0277523,
    'Th': 230.0331341,
    'Pa': 231.0358842,
    'U': 233.0396355,
    'Np': 236.04657,
    'Pu': 238.0495601,
    'Am': 241.0568293,
    'Cm': 243.0613893,
    'Bk': 247.0703073,
    'Cf': 249.0748539,
    'Es': 252.08298,
    'Fm': 257.0951061,
    'Md': 258.0984315,
    'No': 259.10103,
    'Lr': 262.10961,
    'Rf': 267.12179,
    'Db': 268.12567,
    'Hs': 270.13429,
    'Sg': 271.13393,
    'Bh': 272.13826,
    'Mt': 276.15159,
    'Rg': 280.16514,
    'Ds': 281.16451,
    'Uut': 284.17873,
    'Cn': 285.17712,
    'Uup': 288.19274,
    'Fl': 289.19042,
    'Uus': 292.20746,
    'Lv': 293.20449,
    'Uuo': 294.21392,
}


def get_molecular_mass(formula):
    '''Calculate and return molecular mass from chemical formula.'''
    return sum([__ELEMENTAL_MASSES[element] * count
                if element in __ELEMENTAL_MASSES
                else float('NaN')
                for element, count in get_elem_comp(formula).iteritems()])


def get_elem_comp(formula):
    '''Gets elemental composition as a dict from formula.'''
    elem_comp = {}

    for term in re.findall('[A-Z]{1}[0-9]*[a-z]{0,1}[0-9]*', formula):
        element = re.search('[A-z]*', term).group(0)
        result = re.search('[0-9]+', term)
        elem_comp[element] = int(result.group(0)) if result else 1

    return elem_comp


def parse_equation(equation, separator='='):
    '''Parses chemical equation strings.'''
    equation_terms = [re.split('\\s+\\+\\s+', equation_side)
                      for equation_side in
                      re.split('\\s*' + separator + '\\s*', equation)]

    # Add reactants and products:
    return _get_reaction_participants(equation_terms[0], -1) + \
        _get_reaction_participants(equation_terms[1], 1)


def balance(reaction_def, optional_comp=None, max_stoich=8.0):
    '''Applies linear programming to balance reaction.'''
    if optional_comp is None:
        optional_comp = [('H', 1, 'CHEBI:24636'), ('H2O', 0, 'CHEBI:15377')]

    all_formulae = [[x[0] for x in reaction_def if x[2] <= 0],
                    [x[0] for x in reaction_def if x[2] > 0]]
    all_formulae.extend([[x[0] for x in optional_comp] for _ in range(2)])

    all_charges = [[x[1] for x in reaction_def if x[2] <= 0],
                   [x[1] for x in reaction_def if x[2] > 0]]
    all_charges.extend([[x[1] for x in optional_comp] for _ in range(2)])

    all_ids = [[x[3] for x in reaction_def if x[2] <= 0],
               [x[3] for x in reaction_def if x[2] > 0]]
    all_ids.extend([[x[2] for x in optional_comp] for _ in range(2)])

    all_elem_comp = [[_get_elem_comp(formula, idx)
                      for formula in formulae]
                     for idx, formulae in enumerate(all_formulae)]

    (balanced, stoichs) = _optimise(_get_elem_matrix(all_elem_comp,
                                                     all_charges),
                                    max_stoich, len(optional_comp))

    balanced_def = sorted(_get_reaction_def(stoichs, all_formulae,
                                            all_charges, all_ids)) \
        if balanced else reaction_def

    return balanced, \
        balanced and _compare_reaction_defs(reaction_def, balanced_def), \
        balanced_def


def _get_reaction_participants(equation_term, stoich_factor):
    '''Adds reaction participants to a list of participants.'''
    if len(equation_term) == 1 and len(equation_term[0]) == 0:
        return []

    all_terms = [participant.split() for participant in equation_term]
    return [(terms[0], stoich_factor)
            if len(terms) == 1
            else (terms[1], stoich_factor * float(terms[0]))
            for terms in all_terms]


def _get_elem_comp(formula, idx):
    '''Returns elemental compositions, multiplying stoichiometry by -1
    if reactants.'''
    elem_comp = get_elem_comp(formula)
    elem_comp.update((x, y * (-1 if idx % 2 == 0 else 1))
                     for x, y in elem_comp.items())
    return elem_comp


def _get_elem_matrix(all_elem_comp, all_charges):
    '''Gets the elemental (and charge) matrix representing rows of elements
    and columns of compounds.'''
    a_matrix = []
    elements = [elem_comp.keys() for elem_comps in all_elem_comp
                for elem_comp in elem_comps]

    for element in set([item for sublist in elements for item in sublist]):
        a_matrix.append([elem_comp[element] if element in elem_comp else 0
                         for elem_comps in all_elem_comp
                         for elem_comp in elem_comps])

    a_matrix.append(list(itertools.chain(*([charge * (-1 if idx % 2 == 0
                                                      else 1)
                                            for charge in charges]
                                           for idx, charges
                                           in enumerate(all_charges)))))

    return a_matrix


def _optimise(a_matrix, max_stoich, num_opt_formulae):
    '''Optimise linear program and return result.'''
    bounds = [(1.0, max_stoich)] * (len(a_matrix[0]) -
                                    num_opt_formulae * 2) + \
        [(0.0, max_stoich)] * num_opt_formulae * 2

    return math_utils.linprog([1] * len(a_matrix[0]),
                              a_matrix,
                              [0.0] * (len(a_matrix)),
                              bounds)


def _get_reaction_def(stoichs, all_formulae, all_charges, all_ids):
    '''Formats the input into (formula, charge, stoichiometry) tuples.'''
    return [(a[0], b, a[1] * c, d)
            for a, b, c, d in zip([(x, -1 if idx % 2 == 0 else 1)
                                   for idx, formulae in enumerate(all_formulae)
                                   for x in formulae],
                                  [y for charges in all_charges
                                   for y in charges],
                                  _simplify_stoichs(stoichs),
                                  [z for ids in all_ids for z in ids])
            if c > 1e-8]


def _compare_reaction_defs(def1, def2):
    '''Returns True/False depending on whether reaction definitions are
    equal.'''
    return def1 is not None and def2 is not None and \
        len(def1) == len(def2) and \
        all([x[0] == y[0] and x[1] == y[1] and math_utils.isclose(x[2], y[2])
             for x, y in zip(sorted(def1), sorted(def2))])


def _simplify_stoichs(stoichs):
    '''Attempts to simplify stoichs of 1.00000001 to 1.0.'''
    return [float(round(stoich)) if math_utils.isclose(stoich, round(stoich))
            else stoich for stoich in stoichs]
