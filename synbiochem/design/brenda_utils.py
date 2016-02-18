'''
SYNBIOCHEM-DB (c) University of Manchester 2015

SYNBIOCHEM-DB is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import hashlib
import sys
from SOAPpy import WSDL


class BrendaReader(object):
    '''Client of BRENDA webservices.'''

    __WSDL = 'http://www.brenda-enzymes.info/soap/brenda.wsdl'

    def __init__(self, username, password, wsdl=__WSDL):
        self.__username = username
        self.__password = hashlib.sha256(password).hexdigest()
        self.__client = WSDL.Proxy(wsdl)

    def get_km_values(self, ec_number, organism=None):
        '''Returns Km values.'''
        parameters = self.__get_parameters(ec_number, organism)
        results = _parse_results(self.__client.getKmValue(parameters))
        ligands = self.__get_ligands(parameters)
        _map_ligand_ids(results, ligands)
        return results

    def get_turnover_number(self, ec_number, organism=None):
        '''Returns kcat values.'''
        parameters = self.__get_parameters(ec_number, organism)
        return _parse_results(self.__client.getTurnoverNumber(parameters))

    def __get_ligands(self, params):
        '''Returns ligand values.'''
        return {result['ligandStructureId']: result['ligand']
                for result in _parse_results(self.__client.getLigands(params))}

    def __get_parameters(self, ec_number, organism=None):
        '''Generates query string.'''
        return self.__username + ',' + self.__password + \
            ',ecNumber*' + ec_number + \
            ('#organism*' + organism if organism is not None else '')


def _parse_results(result_str):
    '''Parses BRENDA-formatted results string into something sensible.'''
    return [{val[0]: val[1] for val in [item.split("*")
                                        for item in line.split("#")]
             if len(val) == 2}
            for line in result_str.split('!')]


def _map_ligand_ids(results, ligands):
    '''Maps ligand id to ligand name.'''
    for result in results:
        ligand_id = result.pop('ligandStructureId')

        if len(ligand_id) > 0 and ligand_id in ligands:
            result['ligand'] = ligands[ligand_id]


def main(argv):
    '''main method'''
    brenda_reader = BrendaReader(argv[0], argv[1])
    km_values = brenda_reader.get_km_values(argv[2], argv[3])
    print km_values


if __name__ == '__main__':
    main(sys.argv[1:])
