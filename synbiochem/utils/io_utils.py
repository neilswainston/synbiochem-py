'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import gzip
import os
import tempfile
import zipfile
import requests


def get_file(source_url, target_filename):
    '''Downloads file from url and saves to target.'''

    if not os.path.isfile(target_filename):
        if not os.path.exists(os.path.dirname(target_filename)):
            os.makedirs(os.path.dirname(target_filename))

        resp = requests.get(source_url, allow_redirects=True)

        with open(target_filename, 'w') as target_file:
            target_file.write_bytes(resp.content)

    destination = os.path.dirname(target_filename)

    if target_filename.endswith('.zip'):
        zfile = zipfile.ZipFile(target_filename, 'r')
        target_filename = os.path.join(destination, zfile.namelist()[0])
        zfile.extractall(destination)
    elif target_filename.endswith('.gz'):
        unzipped_filepath = target_filename[:-len('.gz')]

        if os.path.exists(unzipped_filepath):
            target_filename = unzipped_filepath
        else:
            input_file = gzip.open(target_filename, 'rb')
            target_filename = os.path.join(destination,
                                           input_file.name[:-len('.gz')])
            output_file = open(target_filename, 'wb')

            for line in input_file:
                output_file.write(line)

            input_file.close()
            output_file.close()

    return target_filename


def get_filename(filename):
    '''Returns a filename, generating a temp file if necessary.'''
    if filename is None:
        fle = tempfile.NamedTemporaryFile(delete=False)
        return fle.name

    return filename


def get_filenames(filepaths, max_files=int(1e16)):
    '''Get filename.'''
    all_filenames = []

    for filepath in filepaths:
        all_filenames.extend(_get_filenames(filepath))

    return all_filenames[:max_files]


def _get_filenames(filepath):
    '''Get filename.'''
    filenames = []

    if os.path.isdir(filepath):
        for filename in os.listdir(os.path.abspath(filepath)):
            filenames.append(os.path.join(filepath, filename))

        return filenames

    # else:
    return [filepath]
