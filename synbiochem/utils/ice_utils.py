'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  pablocarbonell / neilswainston / alanwilliams
'''
# pylint: disable=too-many-arguments
# pylint: disable=too-many-instance-attributes
import json
import tempfile

from synbiochem.utils import net_utils as net_utils, dna_utils


_DEFAULT_ID_PREFIX = 'SBC'
_SESSION_KEY = 'X-ICE-Authentication-SessionId'


class ICEEntry(object):
    '''Class to represent an ICE entry.'''

    def __init__(self, dna=None, typ=None, metadata=None):

        assert typ is not None or metadata is not None

        self.__dna = dna
        self.__dna_updated = dna is not None

        if metadata is None:
            self.__metadata = {'type': typ}
        else:
            if 'type' not in metadata:
                metadata['type'] = typ

            self.__metadata = metadata

    def get_ice_number(self):
        '''Gets the ICE number.'''
        return self.__metadata['id'] if 'id' in self.__metadata else None

    def get_ice_id(self):
        '''Gets the ICE id.'''
        ice_number = self.get_ice_number()
        return get_ice_id(ice_number) if ice_number is not None else None

    def get_record_id(self):
        '''Gets the ICE record id.'''
        return self.__metadata['recordId'] if 'recordId' in self.__metadata \
            else None

    def get_type(self):
        '''Gets the ICE type.'''
        return self.__metadata['type']

    def get_name(self):
        '''Gets the ICE name.'''
        return self.__metadata['name'] if 'name' in self.__metadata else ''

    def get_metadata(self):
        '''Gets the metadata.'''
        return self.__metadata

    def get_dna(self):
        '''Gets the DNA object.'''
        return self.__dna

    def get_dna_updated(self):
        '''Gets the DNA object updated flag.'''
        return self.__dna_updated

    def set_values(self, new_metadata):
        '''Sets multiple metadata values.'''
        self.__metadata.update(new_metadata)

    def set_value(self, key, value):
        '''Sets a metadata value.'''
        self.__metadata[key] = value

    def set_dna(self, dna):
        '''Sets the DNA object.'''
        self.__dna_updated = self.__dna is not None \
            or dna is not None
        self.__dna = dna

    def unset_dna_updated(self):
        '''Sets the DNA object updated flag.'''
        self.__dna_updated = False

    def __repr__(self):
        return str(self.__metadata) + \
            ('\n' + str(self.__dna)
             if self.__dna is not None
             else '')


class ICEClient(object):
    '''Class representing an ICE client.'''

    def __init__(self, url, username, psswrd, id_prefix=_DEFAULT_ID_PREFIX):
        self.__url = url[:-1] if url[-1] == '/' else url
        self.__username = username
        self.__psswrd = psswrd
        self.__id_prefix = id_prefix

        self.__headers = {'Accept': 'application/json',
                          'Content-Type': 'application/json'}

        self.__sid, self.__user, self.__email = self.reconnect()
        self.__headers[_SESSION_KEY] = self.__sid

    def reconnect(self):
        '''Reconnects to ICE server.'''
        try:
            resp = self.__get_access_token(
                '/accesstoken', self.__username, self.__psswrd)
        except net_utils.NetworkError:
            resp = self.__get_access_token(
                '/accesstokens', self.__username, self.__psswrd)

        return resp['sessionId'], \
            resp['firstName'] + ' ' + resp['lastName'], \
            resp['email']

    def get_ice_entry(self, ice_id):
        '''Gets an ICEEntry object from the ICE database.'''
        metadata = self.__get_meta_data(ice_id)
        dna = self.__get_dna(ice_id) if metadata['hasSequence'] \
            else None

        return ICEEntry(dna, metadata=metadata)

    def get_dna(self, ice_id):
        '''Gets a DNA object from the ICE database.'''
        ice_entry = self.get_ice_entry(ice_id)
        return ice_entry.get_dna()

    def set_ice_entry(self, ice_entry):
        '''Saves an ICEEntry object in the ICE database.'''
        if ice_entry.get_ice_number() is None:
            response = self.__create_entry(ice_entry)
        else:
            response = self.__update_entry(ice_entry.get_ice_number(),
                                           self.__form_metadata(ice_entry))

        metadata = self.__get_meta_data(self.__get_ice_id(response['id']))
        ice_entry.set_values(metadata)

        if ice_entry.get_dna_updated():

            if 'hasSequence' in metadata and metadata['hasSequence']:
                self.__delete_seq(ice_entry.get_ice_number())

            dna = ice_entry.get_dna()

            if dna is not None:
                self.__upload_dna(ice_entry.get_record_id(),
                                  ice_entry.get_type(),
                                  dna)

                ice_entry.unset_dna_updated()

        metadata = self.__get_meta_data(self.__get_ice_id(response['id']))
        ice_entry.set_values(metadata)
        return response['id']

    def do_blast(self, seq):
        '''Performs BLAST search against database.'''
        data = {'blastQuery': {'blastProgram': 'BLAST_N',
                               'sequence': seq.lower()}}
        return _read_resp(net_utils.post(self.__url + '/rest/search',
                                         json.dumps(data), self.__headers))

    def get_ice_entries_by_seq(self, seq):
        '''Returns entries matching supplied sequence.'''
        entries = []
        response = self.do_blast(seq)

        if 'results' in response:
            for result in response['results']:
                if result['queryLength'] == int(result['alignment']):
                    entry = self.get_ice_entry(result['entryInfo']['id'])

                    if entry.get_dna().sequence == seq:
                        entries.append(entry)

        return entries

    def add_permission(self, ice_id, group_number, read=True):
        '''Adds user permissions to a given ICE entry.'''
        url = self.__url + '/parts/' + self.__get_ice_number(ice_id) + \
            '/permissions'
        data = {'type': 'READ_ENTRY' if read else 'WRITE_ENTRY',
                'article': 'GROUP',
                'articleId': group_number}

        return _read_resp(net_utils.post(url, json.dumps(data),
                                         self.__headers))

    def __get_access_token(self, service, username, psswrd):
        '''Gets access token response.'''
        return _read_resp(net_utils.post(self.__url + '/rest' + service,
                                         json.dumps({'email': username,
                                                     'password': psswrd}),
                                         self.__headers))

    def __form_metadata(self, ice_entry):
        '''Forms metadata dictionary.'''
        metadata = ice_entry.get_metadata()

        if 'creator' not in metadata:
            metadata['creator'] = self.__user

        if 'creatorEmail' not in metadata:
            metadata['creatorEmail'] = self.__email

        return metadata

    def __upload_seq_file(self, record_id, typ, filename):
        '''Uploads a sequence file (not necessarily SBOL).'''
        return _read_resp(net_utils.post_file(self.__url +
                                              '/rest/file/sequence',
                                              {'file': open(filename, 'r'),
                                               'entryType': typ,
                                               'entryRecordId': record_id},
                                              {_SESSION_KEY: self.__sid}))

    def __delete_seq(self, ice_number):
        '''Deletes the sequence associated with supplied record id.'''
        net_utils.delete(self.__url + '/rest/parts/' + str(ice_number) +
                         '/sequence',
                         headers={_SESSION_KEY: self.__sid})

    def __get_meta_data(self, ice_id):
        '''Returns an ICE entry metadata.'''
        return _read_resp(net_utils.get(
            self.__url + '/rest/parts/' + self.__get_ice_number(ice_id),
            self.__headers))

    def __get_dna(self, ice_id):
        '''Gets the sequence ICE entry.'''
        url = self.__url + '/rest/file/' + self.__get_ice_number(ice_id) + \
            '/sequence/sbol1?sid=' + self.__sid
        temp_file = tempfile.NamedTemporaryFile(delete=False)

        with open(temp_file.name, 'w') as text_file:
            text_file.write(net_utils.get(url))

        return dna_utils.read(temp_file.name)

    def __create_entry(self, ice_entry):
        '''Creates a new ICE entry in the database.'''
        url = self.__url + '/rest/parts'
        return _read_resp(
            net_utils.post(url, json.dumps(self.__form_metadata(ice_entry)),
                           self.__headers))

    def __update_entry(self, ice_id, metadata):
        '''Updates an ICE entry in the database.'''
        ice_number = self.__get_ice_number(ice_id)
        url = self.__url + '/rest/parts/' + str(ice_number)

        return _read_resp(net_utils.put(url, json.dumps(metadata),
                                        self.__headers))

    def __upload_dna(self, record_id, typ, dna):
        '''Uploads an SBOLDocument to ICE database.'''
        sbol_file = tempfile.NamedTemporaryFile(suffix='.xml')
        dna_utils.write(dna, sbol_file.name)
        return self.__upload_seq_file(record_id, typ, sbol_file.name)

    def __get_ice_number(self, ice_identifier):
        '''Maps ICE number to ICE id, i.e. from SBC000123 to 123,
        or if a number is supplied, returns the number.'''
        return get_ice_number(ice_identifier, self.__id_prefix)

    def __get_ice_id(self, ice_number):
        '''Maps ICE id to ICE number, i.e. from 123 to SBC000123.'''
        return get_ice_id(ice_number, self.__id_prefix)


def get_ice_number(ice_identifier, id_prefix=_DEFAULT_ID_PREFIX):
    '''Maps ICE number to ICE id, i.e. from SBC000123 to 123,
    or if a number is supplied, returns the number.'''
    try:
        ice_number = int(ice_identifier.replace(id_prefix, ''))
    except AttributeError:
        # "Ask forgiveness, not permission" and assume ice_identifier is
        # the ice_number:
        ice_number = ice_identifier

    return str(ice_number)


def get_ice_id(ice_number, id_prefix=_DEFAULT_ID_PREFIX):
    '''Maps ICE id to ICE number, i.e. from 123 to SBC000123.'''
    return id_prefix + format(ice_number, '06')


def _read_resp(response):
    '''Parses a string response into json.'''
    return json.loads(str(response))
