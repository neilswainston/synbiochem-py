'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import csv
import traceback
import urllib2
import synbiochem.utils.chem_utils as chem_utils

_METANETX_URL = 'http://metanetx.org/cgi-bin/mnxget/mnxref/'


class MnxRefReader(object):
    '''Class to read MnxRef data from the chem_prop.tsv, the chem_xref.tsv and
    reac_prop.tsv files.'''

    def __init__(self, source=_METANETX_URL):
        self.__source = source
        self.__chem_data = {}
        self.__reac_data = {}

    def get_chem_data(self):
        '''Gets chemical data.'''
        if len(self.__chem_data) == 0:
            self.__read_chem_prop()
            self.__read_chem_xref()

        return self.__chem_data

    def get_reac_data(self):
        '''Gets reaction data.'''
        if len(self.__reac_data) == 0:
            self.__read_reac_prop()

        return self.__reac_data

    def __read_chem_prop(self):
        '''Read chemical properties and create Nodes.'''
        chem_prop_keys = ['id', 'name', 'formula', 'charge', 'mass', 'inchi',
                          'smiles', 'source']

        for values in self.__read_data('chem_prop.tsv'):
            if not values[0].startswith('#'):
                props = dict(zip(chem_prop_keys, values))
                self.__chem_data[values[0]] = props

    def __read_chem_xref(self):
        '''Read chemical xrefs and update Nodes.'''
        chem_xref_keys = ['XREF', 'MNX_ID', 'Evidence', 'Description']

        for values in self.__read_data('chem_xref.tsv'):
            if not values[0].startswith('#'):
                xrefs = dict(zip(chem_xref_keys, values))
                xref = xrefs['XREF'].split(':')

                if xrefs['MNX_ID'] in self.__chem_data:
                    chem = self.__chem_data[xrefs['MNX_ID']]
                    chem[xref[0]] = xref[1]

    def __read_reac_prop(self):
        '''Read reaction properties and create Nodes.'''
        equation_key = 'equation'
        reac_prop_keys = ['id', equation_key, 'description', 'balance', 'ec',
                          'Source']

        for values in self.__read_data('reac_prop.tsv'):
            if not values[0].startswith('#'):
                props = dict(zip(reac_prop_keys, values))
                self.__reac_data[values[0]] = props

                try:
                    participants = chem_utils.parse_equation(
                        props[equation_key])

                    for participant in participants:
                        if participant[0] not in self.__chem_data:
                            self.__add_chem(participant[0])
                except ValueError:
                    print traceback.print_exc()

    def __add_chem(self, chem_id):
        '''Adds a chemical with given id.'''
        props = {'id': chem_id}
        self.__chem_data[chem_id] = props
        return props

    def __read_data(self, filename):
        '''Downloads and reads tab-limited files into lists of lists of
        strings.'''
        return list(csv.reader(urllib2.urlopen(self.__source + filename),
                               delimiter='\t'))


def main():
    '''main method'''
    reader = MnxRefReader()
    print len(reader.get_chem_data())
    print len(reader.get_reac_data())


if __name__ == "__main__":
    main()
