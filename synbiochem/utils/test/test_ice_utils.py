'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=no-member
# pylint: disable=too-many-public-methods
import getpass
import unittest

from synbiochem.utils.ice_utils import ICEClient, ICEEntry
import synbiochem.utils.dna_utils as dna_utils
import synbiochem.utils.ice_utils as ice_utils


class TestICEEntry(unittest.TestCase):
    '''Test class for ICEEntry.'''
    @classmethod
    def setUpClass(cls):
        ice_url = raw_input('ICE url: ')
        username = raw_input('ICE username: ')
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
        dna = dna_utils.read('sbol.xml')
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
        dna1 = dna_utils.read('sbol.xml')
        dna2 = dna_utils.read('sbol2.xml')

        ice_entry = ICEEntry(typ='PLASMID', dna=dna1)
        self.__ice_client.set_ice_entry(ice_entry)

        ice_entry.set_dna(dna2)
        self.__ice_client.set_ice_entry(ice_entry)

        self.assertEqual(ice_entry.get_dna().sequence, dna2.sequence)


class TestICEClient(unittest.TestCase):
    '''Test class for ICEClient.'''
    @classmethod
    def setUpClass(cls):
        ice_url = raw_input('ICE url: ')
        username = raw_input('ICE username: ')
        password = getpass.getpass(prompt='ICE password: ')
        cls.__ice_client = ICEClient(ice_url, username, password)

    def test_get_ice_entry(self):
        '''Tests get_ice_entry method.'''
        dna = dna_utils.read('sbol.xml')

        ice_entry_in = ICEEntry(typ='PART', dna=dna)
        self.__ice_client.set_ice_entry(ice_entry_in)

        ice_entry_out = self.__ice_client.get_ice_entry(
            ice_entry_in.get_ice_number())
        self.assertEqual(ice_entry_out.get_dna().sequence, dna.sequence)

    def test_set_ice_entry(self):
        '''Tests set_ice_entry method.'''
        dna_in = dna_utils.read('sbol.xml')

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

    def test_get_ice_entries_by_seq(self):
        '''Tests get_ice_entries_by_seq method.'''
        dna = dna_utils.read('sbol.xml')
        # dna.set_seq(sequence_utils.get_random_dna(4096))

        ice_entry = ICEEntry(typ='PLASMID', dna=dna)
        self.__ice_client.set_ice_entry(ice_entry)

        self.__ice_client.reconnect()
        result = self.__ice_client.get_ice_entries_by_seq(dna.sequence)
        self.assertTrue(len(result) > 0)


class Test(unittest.TestCase):
    '''Test class for ice_utils.'''

    def test_get_ice_number(self):
        '''Tests get_ice_number method.'''
        self.assertEqual(ice_utils.get_ice_number('TEST00063', 'TEST'), '63')
        self.assertEqual(ice_utils.get_ice_number('TEST63', 'TEST'), '63')
        self.assertEqual(ice_utils.get_ice_number('63', 'TEST'), '63')
        self.assertEqual(ice_utils.get_ice_number(63, 'TEST'), '63')


if __name__ == "__main__":
    unittest.main()
