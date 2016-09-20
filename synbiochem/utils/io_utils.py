'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import gzip
import os
import tempfile
import urllib
import zipfile


def get_file(source_url, target_filename):
    '''Downloads file from url and saves to target.'''

    if not os.path.isfile(target_filename):
        if not os.path.exists(os.path.dirname(target_filename)):
            os.makedirs(os.path.dirname(target_filename))

        urllib.urlretrieve(source_url, target_filename)

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
