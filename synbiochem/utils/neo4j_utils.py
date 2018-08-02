'''
synbiochem (c) University of Manchester 2017

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=too-many-arguments
import functools
import os
import shutil
import subprocess

import numpy as np
import pandas as pd


def create_db(neo4j_root_loc, nodes_files, rels_files,
              delimiter, array_delimiter,
              db_name='graph.db',
              ignore_duplicates=False):
    '''Creates the database from csv files.

    Import format is:

    neo4j-admin import --nodes NODES_FILE --relationships RELS_FILE

    '''
    # Apply types to import files:
    for filename in nodes_files + rels_files:
        df = pd.read_csv(filename, sep=delimiter,
                         low_memory=False, dtype=object)
        df = type_df(df, array_delimiter)
        df.to_csv(filename, delimiter, encoding='utf-8', index=False)

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


def type_df(df, array_delimiter):
    '''Format dataframe and update column names in file based on dtype.'''
    # Id / special case columns regexps:
    _format_list_delim = functools.partial(_format_list,
                                           array_delimiter=array_delimiter)
    new_df = pd.DataFrame()

    for col in df:
        if ':' in col:
            new_df[col] = df[col]
        else:
            # Check if list and reformat if so:
            res = df[col].apply(_format_list_delim)
            dtype = _get_type(res[0].unique())

            if dtype:
                new_df[col + ':' + dtype + '[]'] = res[1]
            else:
                # Convert appropriate cols to int / float / boolean:
                try:
                    dtype = _get_type(str(df[col].dtypes))
                    new_df[col.replace(':', '_') + ':' + dtype] = res[1]
                except (TypeError, ValueError):
                    # Ignore errors if data is not numerical.
                    pass

    return new_df


def _format_list(cell, array_delimiter):
    '''Format cell if of type list.'''
    if not isinstance(cell, list):
        return pd.Series([None, cell])

    return pd.Series([pd.Series(cell).dtype,
                      array_delimiter.join(map(str, cell))])


def _get_type(dtype):
    '''Convert dtype to neo4j import type.'''
    if isinstance(dtype, np.ndarray):
        if len(dtype) == 1:
            return _get_type(str(dtype[0]))

        return None

    if dtype == 'bool':
        return 'boolean'
    if dtype == 'object':
        return 'string'
    if dtype in ['nan', 'None']:
        return None

    return ''.join([c for c in dtype if not c.isdigit()])
