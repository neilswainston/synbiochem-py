'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=no-member
# pylint: disable=redundant-unittest-assert
# pylint: disable=too-many-public-methods
import getpass
import os
import unittest

from synbiochem.utils import ice_utils, sbol_utils
from synbiochem.utils.ice_utils import ICEClient, ICEEntry


class TestICEEntry(unittest.TestCase):
    '''Test class for ICEEntry.'''
    @classmethod
    def setUpClass(cls):
        try:
            std_input = raw_input
        except NameError:
            std_input = input

        ice_url = std_input('ICE url: ')
        username = std_input('ICE username: ')
        password = getpass.getpass(prompt='ICE password: ')
        cls.__ice_client = ICEClient(ice_url, username, password)

    def test_get_ice_number(self):
        '''Tests get_ice_number method.'''
        ice_entry1 = ICEEntry(typ='PLASMID')
        self.__ice_client.set_ice_entry(ice_entry1)
        self.assertNotEqual(ice_entry1.get_ice_number(), None)

        ice_entry2 = ICEEntry(metadata={'type': 'PLASMID'})
        self.__ice_client.set_ice_entry(ice_entry2)
        self.assertNotEqual(ice_entry2.get_ice_number(), None)

    def test_get_ice_number_none(self):
        '''Tests get_ice_number method.'''
        ice_entry1 = ICEEntry(typ='PLASMID')
        self.assertEqual(ice_entry1.get_ice_number(), None)
        ice_entry2 = ICEEntry(metadata={'type': 'PLASMID'})
        self.assertEqual(ice_entry2.get_ice_number(), None)

    def test_get_record_id(self):
        '''Tests get_record_id method.'''
        ice_entry1 = ICEEntry(typ='PLASMID')
        self.__ice_client.set_ice_entry(ice_entry1)
        self.assertNotEqual(ice_entry1.get_record_id(), None)

        ice_entry2 = ICEEntry(metadata={'type': 'PLASMID'})
        self.__ice_client.set_ice_entry(ice_entry2)
        self.assertNotEqual(ice_entry2.get_record_id(), None)

    def test_get_record_id_none(self):
        '''Tests get_record_id method.'''
        ice_entry1 = ICEEntry(typ='PLASMID')
        self.assertEqual(ice_entry1.get_record_id(), None)
        ice_entry2 = ICEEntry(metadata={'type': 'PLASMID'})
        self.assertEqual(ice_entry2.get_record_id(), None)

    def test_get_type(self):
        '''Tests get_type method.'''
        ice_entry1 = ICEEntry(typ='PLASMID')
        self.assertEqual(ice_entry1.get_type(), 'PLASMID')
        ice_entry2 = ICEEntry(metadata={'type': 'PLASMID'})
        self.assertEqual(ice_entry2.get_type(), 'PLASMID')

    def test_get_metadata(self):
        '''Tests get_metadata method.'''
        ice_entry1 = ICEEntry(typ='PLASMID')
        self.assertEqual(ice_entry1.get_metadata()['type'], 'PLASMID')
        ice_entry2 = ICEEntry(metadata={'type': 'PLASMID'})
        self.assertEqual(ice_entry2.get_metadata()['type'], 'PLASMID')

    def test_get_sbol_doc(self):
        '''Tests get_sbol_doc method.'''
        dna = _read('sbol.xml')
        ice_entry = ICEEntry(typ='PLASMID', dna=dna)
        self.assertEqual(ice_entry.get_dna(), dna)

    def test_set_value(self):
        '''Tests set_value method.'''
        ice_entry1 = ICEEntry(typ='PLASMID')
        ice_entry1.set_value('creator', 'God')
        self.assertEqual(ice_entry1.get_metadata()['type'], 'PLASMID')
        self.assertEqual(ice_entry1.get_metadata()['creator'], 'God')
        self.__ice_client.set_ice_entry(ice_entry1)

        ice_entry2 = self.__ice_client.get_ice_entry(
            ice_entry1.get_ice_number())
        ice_entry2.set_value('creator', 'Aitor Karanka')
        self.__ice_client.set_ice_entry(ice_entry2)

        ice_entry3 = self.__ice_client.get_ice_entry(
            ice_entry1.get_ice_number())
        self.assertEqual(ice_entry3.get_metadata()['creator'], 'Aitor Karanka')

        ice_entry4 = ICEEntry(metadata={'type': 'PLASMID'})
        ice_entry4.set_value('creator', 'God')
        self.assertEqual(ice_entry4.get_metadata()['type'], 'PLASMID')
        self.assertEqual(ice_entry4.get_metadata()['creator'], 'God')

    def test_set_values(self):
        '''Tests set_values method.'''
        ice_entry1 = ICEEntry(typ='PLASMID')
        ice_entry1.set_values({'creator': 'God', 'name': 'Test'})
        self.assertEqual(ice_entry1.get_metadata()['type'], 'PLASMID')
        self.assertEqual(ice_entry1.get_metadata()['creator'], 'God')
        self.assertEqual(ice_entry1.get_metadata()['name'], 'Test')
        ice_entry2 = ICEEntry(metadata={'type': 'PLASMID'})
        ice_entry2.set_values({'creator': 'God', 'name': 'Test'})
        self.assertEqual(ice_entry2.get_metadata()['type'], 'PLASMID')
        self.assertEqual(ice_entry2.get_metadata()['creator'], 'God')
        self.assertEqual(ice_entry2.get_metadata()['name'], 'Test')

    def test_set_dna(self):
        '''Tests set_dna method.'''
        dna1 = _read('sbol.xml')
        dna2 = _read('sbol2.xml')

        ice_entry = ICEEntry(typ='PLASMID', dna=dna1)
        self.__ice_client.set_ice_entry(ice_entry)

        ice_entry.set_dna(dna2)
        self.__ice_client.set_ice_entry(ice_entry)

        self.assertEqual(ice_entry.get_seq(), dna2['seq'])

    def test_get_seq(self):
        '''Tests get_seq method.'''
        dna = _read('sbol.xml')

        ice_entry = ICEEntry(typ='PLASMID', dna=dna)
        self.__ice_client.set_ice_entry(ice_entry)
        ice_entry = self.__ice_client.get_ice_entry(ice_entry.get_ice_id())

        self.assertEqual(ice_entry.get_seq(), dna['seq'])

    def test_set_get_parameter(self):
        '''Tests set/get_parameter method.'''

        ice_entry = ICEEntry(typ='PLASMID')
        ice_entry.set_parameter('Cheese', 'Brie')
        self.__ice_client.set_ice_entry(ice_entry)
        ice_entry = self.__ice_client.get_ice_entry(ice_entry.get_ice_id())

        self.assertEqual(ice_entry.get_parameter('Cheese'), 'Brie')


class TestICEClient(unittest.TestCase):
    '''Test class for ICEClient.'''
    @classmethod
    def setUpClass(cls):
        try:
            std_input = raw_input
        except NameError:
            std_input = input

        ice_url = std_input('ICE url: ')
        username = std_input('ICE username: ')
        password = getpass.getpass(prompt='ICE password: ')
        cls.__ice_client = ICEClient(ice_url, username, password)

    def test_get_ice_entry(self):
        '''Tests get_ice_entry method.'''
        dna = _read('sbol.xml')

        ice_entry_in = ICEEntry(typ='PART', dna=dna)
        self.__ice_client.set_ice_entry(ice_entry_in)

        ice_entry_out = self.__ice_client.get_ice_entry(
            ice_entry_in.get_ice_number())
        self.assertEqual(ice_entry_out.get_seq(), dna['seq'])

    def test_set_ice_entry(self):
        '''Tests set_ice_entry method.'''
        dna_in = _read('sbol.xml')

        ice_entry_in = ICEEntry(typ='PLASMID', dna=dna_in)
        self.__ice_client.set_ice_entry(ice_entry_in)

        ice_entry_out = self.__ice_client.get_ice_entry(
            ice_entry_in.get_ice_number())

        self.assertNotEqual(ice_entry_out, None)

    def test_do_blast(self):
        '''Tests do_blast method.'''
        result = self.__ice_client.do_blast('tcgagaattcaaaagatctgagataggtaga' +
                                            'agctagacgagaaaccttcccaatcttatca' +
                                            'ttacgaaaggacgtccctatgagcctgatta')
        self.assertTrue(result['resultCount'] > 0)

    def test_do_blast_2(self):
        '''Tests do_blast method.'''
        result = self.__ice_client.do_blast('aggcaaattcagtgaggctgacttctcatct' +
                                            'taaatagttcccttcacgatagccgcctga')
        self.assertTrue(result['resultCount'] > 0)

    def test_get_ice_entries_by_seq(self):
        '''Tests get_ice_entries_by_seq method.'''
        dna = _read('sbol.xml')
        # dna.set_seq(sequence_utils.get_random_dna(4096))

        ice_entry = ICEEntry(typ='PLASMID', dna=dna)
        self.__ice_client.set_ice_entry(ice_entry)

        self.__ice_client.reconnect()
        result = self.__ice_client.get_ice_entries_by_seq(dna['seq'])
        self.assertTrue(len(result) > 0)

    def test_get_ice_entries_by_seq_2(self):
        '''Tests get_ice_entries_by_seq method.'''
        seq = 'aggcaaattcagtgaggctgacttctcatcttaaatagttcccttcacgatagccgcctga'
        result = self.__ice_client.get_ice_entries_by_seq(seq)
        self.assertTrue(len(result) > 0)

    def test_add_link(self):
        '''Tests add_link method.'''
        dna1 = _read('sbol.xml')
        dna2 = _read('sbol2.xml')

        ice_entry1 = ICEEntry(typ='PLASMID', dna=dna1)
        self.__ice_client.set_ice_entry(ice_entry1)

        ice_entry2 = ICEEntry(typ='PLASMID', dna=dna2)
        self.__ice_client.set_ice_entry(ice_entry2)

        self.__ice_client.add_link(
            ice_entry1.get_ice_number(), ice_entry2.get_ice_number())

        # "Refresh" (update metadata)
        ice_entry1 = self.__ice_client.get_ice_entry(ice_entry1.get_ice_id())
        ice_entry2 = self.__ice_client.get_ice_entry(ice_entry2.get_ice_id())

        self.assertTrue(ice_entry1.get_ice_number() in
                        [par['id']
                         for par in ice_entry2.get_metadata()['parents']])

        self.assertTrue(ice_entry2.get_ice_number() in
                        [link['id']
                         for link in ice_entry1.get_metadata()['linkedParts']])

    def test_search_groups(self):
        '''Tests get_group_id and search_groups method.'''
        groups = self.__ice_client.get_groups()

        for name in groups:
            groups = self.__ice_client.search_groups(name[:-1])
            self.assertTrue(name in [grp['label'] for grp in groups])

    def test_add_permission(self):
        '''Tests add_permission method.'''
        dna = _read('sbol.xml')
        ice_entry = ICEEntry(typ='PLASMID', dna=dna)
        self.__ice_client.set_ice_entry(ice_entry)

        groups = self.__ice_client.get_groups()

        for group_num in groups.values():
            self.__ice_client.add_permission(ice_entry.get_ice_id(), group_num)

        # If test progresses to here, it has succeeded:
        self.assertTrue(True)

    def test_search(self):
        '''Tests search method.'''
        resp = self.__ice_client.search('PLASMID')

        # If test progresses to here, it has succeeded:
        self.assertTrue(resp['resultCount'] > 0)

    def test_search_name(self):
        '''Tests advanced search by name method.'''
        resp = self.__ice_client.search_name('SBC_DE15_PL01', 'PLASMID')

        # If test progresses to here, it has succeeded:
        self.assertTrue(resp)

    def test_search_design(self):
        '''Tests advanced search by name method.'''
        resp = self.__ice_client.search_design(15)

        # If test progresses to here, it has succeeded:
        self.assertTrue(len(resp) == 46)

    def test_advanced_search(self):
        '''Tests search method.'''
        typ = 'PLASMID'
        resp = self.__ice_client.advanced_search('*', typ)

        self.assertTrue(len(resp['results']) == 5)

        for result in resp['results']:
            self.assertTrue(result['entryInfo']['type'] == typ)

    def test_get_genbank(self):
        '''Tests get_genbank method.'''
        resp = self.__ice_client.get_genbank(6592)
        self.assertTrue(resp)


class Test(unittest.TestCase):
    '''Test class for ice_utils.'''

    def test_get_ice_number(self):
        '''Tests get_ice_number method.'''
        self.assertEqual(ice_utils.get_ice_number('TEST00063', 'TEST'), '63')
        self.assertEqual(ice_utils.get_ice_number('TEST63', 'TEST'), '63')
        self.assertEqual(ice_utils.get_ice_number('63', 'TEST'), '63')
        self.assertEqual(ice_utils.get_ice_number(63, 'TEST'), '63')


def _read(filename):
    '''Reads sbol file.'''
    directory = os.path.dirname(os.path.realpath(__file__))
    return sbol_utils.read(os.path.join(directory, filename))


if __name__ == "__main__":
    unittest.main()
