#! /usr/bin/python3

import argparse
import os
import re
import sqlite3 as sql
import sys
import xml.etree.cElementTree as et
import traceback

import lib.initialize as initialize
import lib.sqlite_interface as misc
import lib.meta as meta


# ==================
# EXPORTED FUNCTIONS
# ==================

def parse(parent, *args, **kwargs):
    parser = parent.add_parser(
        'blast',
        help="Read BLAST XML report into SQL database",
        parents=args)
    parser.set_defaults(func=parse_blast_xml)

def parse_blast_xml(args, cur):
    try:
        for f in args.input:
            con = et.iterparse(f, events=('end', 'start'))
            _parse_blast_xml(args, cur, con)
    except AttributeError:
        con = et.iterparse(args.input, events=('end', 'start'))
        _parse_blast_xml(args, cur, con)

def _parse_blast_xml(args, cur, con):
    # Initialize tables as necessary
    if(not misc.table_exists('hsp', cur)):
        initialize.init_hsp(cur, verbose=False)
    if(not misc.table_exists('iteration', cur)):
        initialize.init_iteration(cur, verbose=False)
    if(not misc.table_exists('database', cur)):
        initialize.init_database(cur, verbose=False)

    db, qseqid, sseqid, nhit, nhsp = ("", "", "", 0, 0)
    for event, elem in con:
        if event == "end":
            if elem.tag in Hsp.fields:
                print("hsp: {} {}".format(elem.tag, elem.text))
            elif elem.tag in Iteration.fields:
                print("iteration: {} {}".format(elem.tag, elem.text))
            elif elem.tag in Database.fields:
                print("database: {} {}".format(elem.tag, elem.text))

class Table:
    pass

class Hsp(Table):
    fields = {
        'Hit_num',
        'Hit_len',
        'Hsp_num',
        'Hsp_bit_score',
        'Hsp_score',
        'Hsp_evalue',
        'Hsp_query_from',
        'Hsp_query_to',
        'Hsp_hit_from',
        'Hsp_hit_to',
        'Hsp_query_frame',
        'Hsp_hit_frame',
        'Hsp_identity',
        'Hsp_positive',
        'Hsp_align_len',
        'Hsp_gaps'
    }

class Iteration(Table):
    fields = {
        'Iteration_iter_num',
        'Iteration_query_ID',
        'Iteration_query_def',
        'Iteration_query_len',
        'Iteration_message',
        'Statistics_db_num',
        'Statistics_db_len',
        'Statistics_hsp_len',
        'Statistics_eff_space',
        'Statistics_kappa',
        'Statistics_lambda',
        'Statistics_entropy'
    }

class Database(Table):
    fields = {
        'BlastOutput_program',
        'BlastOutput_version',
        'Parameters_matrix',
        'Parameters_expect',
        'Parameters_gap_open',
        'Parameters_gap_extend',
        'Parameters_filter'
    }



        # if(event == 'start'): continue
        # if(elem.tag == 'Hsp'):
        #     bdat.add_partial_row()
        #     bdat.clear_hsp()
        # elif(elem.tag == 'Hit'):
        #     bdat.clear_hit()
        # elif(elem.tag == 'Iteration'):
        #     if(not bdat.has_hits()):
        #         bdat.add_partial_row()
        #     bdat.clear_iter()
        #     elem.clear()
        # elif('BlastOutput_db' in elem.tag):
        #     base = os.path.basename(elem.text)
        #     if(not misc.entry_exists('blastdatabase', 'database', base, cur)):
        #         misc.insert({'database': base}, 'blastdatabase', cur)
        #     bdat.add(elem.tag, base)
        # else:
        #     bdat.add(elem.tag, elem.text)

# class Blastdat:
#     def __init__(self, cur, args):
#         self.cur = cur
#         self.args = args
#         self.dat = {'root':{}, 'iter':{}, 'stat':{}, 'hit':{}, 'hsp':{}}
#         self.dat['root']['collection'] = args.collection
#         self.dat['root']['db_desc'] = args.db_desc
#         self.iter_dicts = []
#         self.row_by_col = {}
#
#     def write_rows_to_sqldb(self):
#         for col in self.row_by_col.keys():
#             misc.insertmany(col, self.row_by_col[col], 'BlastReport',
#                             self.cur, replace=True)
#
#     def has_hits(self):
#         try:
#             if('No hits found' in dat['iter']['Iteration_message']):
#                 return False
#         except:
#             pass
#         return True
#
#     def add_partial_row(self):
#         table = {}
#         for key in self.dat.keys():
#             for subkey in self.dat[key].keys():
#                 table[subkey] = self.dat[key][subkey]
#         self.iter_dicts.append(table)
#
#     def _add_rows(self):
#         if(not self.iter_dicts):
#             self.add_partial_row()
#         else:
#             for d in self.iter_dicts:
#                 if(int(d['Hit_num']) <= self.args.max_hits):
#                     for key, val in self.dat['stat'].items():
#                         d[key] = val
#                     col = tuple(sorted(d.keys()))
#                     row = tuple(map(d.get, col))
#                     try:
#                         self.row_by_col[col].append(row)
#                     except:
#                         self.row_by_col[col] = [row]
#         self.iter_dicts = []
#
#     def clear_iter(self):
#         '''
#         Adds all data from current iteration to the database and frees the
#         iteration and its children hits and hsps from memory
#         '''
#         self._add_rows()
#         self.dat['iter'] = {}
#         self.dat['stat'] = {}
#         self.clear_hit()
#
#     def clear_hit(self):
#         '''
#         Clears the current hit and all children hsps from memory
#         '''
#         self.dat['hit'] = {}
#         self.clear_hsp()
#
#     def clear_hsp(self):
#         '''
#         Clears hsp from memmory
#         '''
#         self.dat['hsp'] = {}
#
#     def add(self, tag, text):
#         '''
#         Input: One tag and its text (possibly None)
#         '''
#         tag = re.sub('-', '_', tag)
#         if(text is None or text.isspace()): pass
#         elif('Hsp_' in tag):
#             if(tag in ('Hsp_qseq', 'Hsp_hseq', 'Hsp_midline') and self.args.small):
#                 pass
#             else:
#                 self.dat['hsp'][tag] = text
#         elif('Hit_' in tag):
#             self.dat['hit'][tag] = text
#         elif('Iteration_' in tag):
#             if(tag == 'Iteration_query_def'):
#                 self.dat['iter']['query_seqid'] = re.sub('(\S+).*', '\\1', text)
#             self.dat['iter'][tag] = text
#         elif('Statistics_' in tag):
#             self.dat['stat'][tag] = text
#         elif('BlastOutput_' in tag or 'Parameters_' in tag):
#             if('reference' in tag or 'query' in tag):
#                 pass
#             else:
#                 self.dat['root'][tag] = text
