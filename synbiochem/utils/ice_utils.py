'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  pablocarbonell / neilswainston / alanwilliams
'''
# pylint: disable=too-many-arguments
import json
import tempfile

import sbol

from synbiochem.utils import net_utils as net_utils, sbol_utils


_DEFAULT_ID_PREFIX = 'SBC'
_SESSION_KEY = 'X-ICE-Authentication-SessionId'


class ICEEntry(object):
    '''Class to represent an ICE entry.'''

    def __init__(self, sbol_doc=None, typ=None, metadata=None):

        assert typ is not None or metadata is not None

        self.__sbol_doc = sbol_doc
        self.__sbol_doc_updated = sbol_doc is not None

        if metadata is None:
            self.__metadata = {'type': typ}
        else:
            self.__metadata = metadata

    def get_ice_number(self):
        '''Gets the ICE number.'''
        return self.__metadata['id'] if 'id' in self.__metadata else None

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

    def get_sbol_doc(self):
        '''Gets the SBOL Document.'''
        return self.__sbol_doc

    def get_sbol_doc_updated(self):
        '''Gets the SBOL Document updated flag.'''
        return self.__sbol_doc_updated

    def set_values(self, new_metadata):
        '''Sets multiple metadata values.'''
        self.__metadata.update(new_metadata)

    def set_value(self, key, value):
        '''Sets a metadata value.'''
        self.__metadata[key] = value

    def set_sbol_doc(self, sbol_doc):
        '''Sets the SBOL Document.'''
        self.__sbol_doc_updated = self.__sbol_doc is not None \
            or sbol_doc is not None
        self.__sbol_doc = sbol_doc

    def unset_sbol_doc_updated(self):
        '''Sets the SBOL Document updated flag.'''
        self.__sbol_doc_updated = False

    def __repr__(self):
        return str(self.__metadata) + \
            ('\n' + str(self.__sbol_doc)
             if self.__sbol_doc is not None
             else '')


class ICEClient(object):
    '''Class representing an ICE client.'''

    def __init__(self, url, username, psswrd, id_prefix=_DEFAULT_ID_PREFIX):
        self.__url = url + ('' if url[-1] == '/' else '/') + 'rest'
        self.__headers = {'Accept': 'application/json',
                          'Content-Type': 'application/json'}

        try:
            resp = self.__get_access_token(
                '/accesstoken', username, psswrd)
        except net_utils.NetworkError:
            resp = self.__get_access_token(
                '/accesstokens', username, psswrd)

        self.__sid = resp['sessionId']
        self.__headers[_SESSION_KEY] = self.__sid
        self.__id_prefix = id_prefix

    def get_ice_entry(self, ice_id):
        '''Gets an ICEEntry object from the ICE database.'''
        metadata = self.__get_meta_data(ice_id)
        sbol_doc = self.__get_sbol_doc(ice_id) if metadata['hasSequence'] \
            else None

        return ICEEntry(sbol_doc, metadata=metadata)

    def set_ice_entry(self, ice_entry):
        '''Saves an ICEEntry object in the ICE database.'''
        if ice_entry.get_ice_number() is None:
            response = self.__create_entry(ice_entry)
        else:
            response = self.__update_entry(ice_entry.get_ice_number(),
                                           ice_entry.get_metadata())

        metadata = self.__get_meta_data(self.__get_ice_id(response['id']))
        ice_entry.set_values(metadata)

        if ice_entry.get_sbol_doc_updated():

            if 'hasSequence' in metadata and metadata['hasSequence']:
                self.__delete_seq(ice_entry.get_ice_number())

            sbol_doc = ice_entry.get_sbol_doc()

            if sbol_doc is not None:
                self.__upload_sbol(ice_entry.get_record_id(),
                                   ice_entry.get_type(),
                                   sbol_doc)

                ice_entry.unset_sbol_doc_updated()

        metadata = self.__get_meta_data(self.__get_ice_id(response['id']))
        ice_entry.set_values(metadata)

    def do_blast(self, seq):
        '''Performs BLAST search against database.'''
        data = {'blastQuery': {'blastProgram': 'BLAST_N',
                               'sequence': seq.lower()}}
        return _read_resp(net_utils.post(self.__url + '/search',
                                         json.dumps(data), self.__headers))

    def get_ice_entries_by_seq(self, seq):
        '''Returns entries matching supplied sequence.'''
        entries = []
        response = self.do_blast(seq)

        if 'results' in response:
            for result in response['results']:
                if '100%' in result['alignment']:
                    entry = self.get_ice_entry(result['entryInfo']['id'])

                    if sbol_utils.get_seq(entry.get_sbol_doc()) == seq:
                        entries.append(entry)

        return entries

    def __get_access_token(self, service, username, psswrd):
        '''Gets access token response.'''
        return _read_resp(net_utils.post(self.__url + service,
                                         json.dumps({'email': username,
                                                     'password': psswrd}),
                                         self.__headers))

    def __upload_seq_file(self, record_id, typ, filename):
        '''Uploads a sequence file (not necessarily SBOL).'''
        return _read_resp(net_utils.post_file(self.__url +
                                              '/file/sequence',
                                              {'file': open(filename, 'r'),
                                               'entryType': typ,
                                               'entryRecordId': record_id},
                                              {_SESSION_KEY: self.__sid}))

    def __delete_seq(self, ice_number):
        '''Deletes the sequence associated with supplied record id.'''
        net_utils.delete(self.__url + '/parts/' + str(ice_number) +
                         '/sequence',
                         headers={_SESSION_KEY: self.__sid})

    def __get_meta_data(self, ice_id):
        '''Returns an ICE entry metadata.'''
        return _read_resp(net_utils.get(
            self.__url + '/parts/' + self.__get_ice_number(ice_id),
            self.__headers))

    def __get_sbol_doc(self, ice_id):
        '''Gets the sequence ICE entry.'''
        url = self.__url + '/file/' + self.__get_ice_number(ice_id) + \
            '/sequence/sbol?sid=' + self.__sid
        temp_file = tempfile.NamedTemporaryFile(delete=False)

        with open(temp_file.name, 'w') as text_file:
            text_file.write(net_utils.get(url))

        document = sbol.Document()
        document.read(temp_file.name)
        return document

    def __create_entry(self, ice_entry):
        '''Creates a new ICE entry in the database.'''
        url = self.__url + '/parts'
        return _read_resp(
            net_utils.post(url, json.dumps(ice_entry.get_metadata()),
                           self.__headers))

    def __update_entry(self, ice_id, metadata):
        '''Updates an ICE entry in the database.'''
        ice_number = self.__get_ice_number(ice_id)
        url = self.__url + '/parts/' + str(ice_number)
        return _read_resp(net_utils.put(url, json.dumps(metadata),
                                        self.__headers))

    def __upload_sbol(self, record_id, typ, sbol_doc):
        '''Uploads an SBOLDocument to ICE database.'''
        sbol_file = tempfile.NamedTemporaryFile()
        sbol_doc.write(sbol_file.name)

        return self.__upload_seq_file(record_id, typ, sbol_file.name)

    def __get_ice_number(self, ice_identifier):
        '''Maps ICE number to ICE id, i.e. from SBC000123 to 123,
        or if a number is supplied, returns the number.'''
        return get_ice_number(ice_identifier, self.__id_prefix)

    def __get_ice_id(self, ice_number):
        '''Maps ICE id to ICE number, i.e. from 123 to SBC000123.'''
        return self.__id_prefix + format(ice_number, '06')


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


def _read_resp(response):
    '''Parses a string response into json.'''
    return json.loads(str(response))
