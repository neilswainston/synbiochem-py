'''
synbiochem (c) University of Manchester 2017

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=too-many-arguments
import os
import re
import shutil
import subprocess

import numpy as np
import pandas as pd


def create_db(neo4j_root_loc, nodes_files, rels_files, db_name='graph.db',
              ignore_duplicates=False, delimiter=';', array_delimiter='|'):
    '''Creates the database from csv files.

    Import format is:

    neo4j-admin import --nodes NODES_FILE --relationships RELS_FILE

    '''
    # Apply types to import files:
    for filename in nodes_files + rels_files:
        _type_file(filename, delimiter)

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

    # not sure why, but this bit is required by ENI model
    # seems like one more bug regarding reading csv
    params.extend(['--multiline-fields', 'True'])

    params.extend(['--delimiter', delimiter])
    params.extend(['--array-delimiter', array_delimiter])

    subprocess.call(params)


def _type_file(filename, sep):
    '''Update column names in file based on dtype.'''
    df = pd.read_csv(filename, sep=sep, low_memory=False)
    _type_cols(df)
    df.to_csv(filename, sep, index=False)


def _type_cols(df):
    '''Update column names based on dtype.'''
    # Id columns regexps:
    id_str = r'(:ID|:START_ID|:END_ID)'
    id_regex = re.compile(id_str)

    # Ensure absence of colons in non-id column names:
    type_cols = {col: col.replace(':', '_')
                 for col in df
                 if col != ':TYPE' and col != ':LABEL' and ':' in col
                 and not id_regex.search(col)}

    # Convert appropriate cols to boolean:
    for col in df:
        if True in pd.unique(df[col]) and df[col].dtype is np.dtype('O'):
            type_cols.update({col: col + ':boolean'})
            df[col] = df[col].astype(str).str.lower()

    # Add types for int and float columns *except* those of ids,
    # which are always assumed to be string:
    col_grouped_by_type = df.columns.to_series().groupby(df.dtypes).groups
    col_grouped_by_type.pop(np.dtype('O'))

    for type_name, cols in col_grouped_by_type.iteritems():
        type_name = ''.join([i for i in str(type_name).lower()
                             if not i.isdigit()])

        type_cols.update({col: col + ':' + type_name
                          for col in cols
                          if not id_regex.search(col)})

    df.rename(columns=type_cols, inplace=True)
