'''
synbiochem (c) University of Manchester 2017

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=too-many-arguments
import os
import shutil
import subprocess


def create_db(neo4j_root_loc, nodes_files, rels_files, db_name='graph.db',
              ignore_duplicates=False, delimiter=',', array_delimiter=';'):
    '''Creates the database from csv files.

    Import format is:

    neo4j-admin import --nodes NODES_FILE --relationships RELS_FILE

    '''
    # Stop database:
    bin_loc = os.path.join(neo4j_root_loc, 'bin')
    subprocess.call([os.path.join(bin_loc, 'neo4j'), 'stop'])

    # Remove database directory if it already exists:
    db_loc = os.path.join(neo4j_root_loc, 'data', 'databases', db_name)

    if os.path.exists(db_loc):
        shutil.rmtree(db_loc)

    # Import database:
    admin_loc = os.path.join(bin_loc, 'neo4j-admin')

    params = [admin_loc,
              'import',
              '--database', db_name]

    for node_file in nodes_files:
        params.extend(['--nodes', node_file])

    for rels_file in rels_files:
        params.extend(['--relationships', rels_file])

    if ignore_duplicates:
        params.append('--ignore-duplicate-nodes')

    params.extend(['--delimiter', delimiter])
    params.extend(['--array-delimiter', array_delimiter])

    subprocess.call(params)
