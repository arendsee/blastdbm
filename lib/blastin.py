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
    parser.add_argument(
        '-c', '--collection',
        metavar="COL",
        help="blast collection")
    parser.add_argument(
        '-m', '--db_desc',
        metavar="DESC",
        help="BLAST database description")
    parser.add_argument(
        '-s', '--small',
        help="Reduce database size by not writing alignment sequences",
        action=('store_true'), default=False)
    parser.add_argument(
        '-x', '--max-hits',
        metavar="INT",
        help='Maximum number of hits to store (default 500)',
        type=int,
        default=500
    )
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
    if(not misc.table_exists('blastreport', cur)):
        initialize.init_blastreport(cur, verbose=False)
    if(not misc.table_exists('blastdatabase', cur)):
        initialize.init_blastdatabase(cur, verbose=False)

    bdat = Blastdat(cur, args)
    for event, elem in con:
        if(event == 'start'): continue
        if(elem.tag == 'Hsp'):
            bdat.add_partial_row()
            bdat.clear_hsp()
        elif(elem.tag == 'Hit'):
            bdat.clear_hit()
        elif(elem.tag == 'Iteration'):
            if(not bdat.has_hits()):
                bdat.add_partial_row()
            bdat.clear_iter()
            elem.clear()
        elif('BlastOutput_db' in elem.tag):
            base = os.path.basename(elem.text)
            if(not misc.entry_exists('blastdatabase', 'database', base, cur)):
                misc.insert({'database': base}, 'blastdatabase', cur)
            bdat.add(elem.tag, base)
        else:
            bdat.add(elem.tag, elem.text)
    bdat.write_rows_to_sqldb()
    meta.update_dbinfo(cur, verbose=True)
    meta.update_mrca(cur, verbose=True)

def _parse_fasta_header(header):
    dic = {}
    try:
        for match in re.finditer('([^|]+)\|([^|]+)', header):
            for tag in ('locus', 'gi', 'taxon', 'gb', 'gene'):
                if(match.group(1) == tag and match.group(2) != None):
                    dic['Query_' + tag] = match.group(2)
        return(dic)
    except:
        print("Cannot parse header {}".format(header), file=sys.stderr)
        return({})


class Blastdat:
    def __init__(self, cur, args):
        self.cur = cur
        self.args = args
        self.dat = {'root':{}, 'iter':{}, 'stat':{}, 'hit':{}, 'hsp':{}}
        self.dat['root']['collection'] = args.collection
        self.dat['root']['db_desc'] = args.db_desc
        self.iter_dicts = []
        self.row_by_col = {}

    def write_rows_to_sqldb(self):
        for col in self.row_by_col.keys():
            misc.insertmany(col, self.row_by_col[col], 'BlastReport',
                            self.cur, replace=True)

    def has_hits(self):
        try:
            if('No hits found' in dat['iter']['Iteration_message']):
                return False
        except:
            pass
        return True

    def add_partial_row(self):
        table = {}
        for key in self.dat.keys():
            for subkey in self.dat[key].keys():
                table[subkey] = self.dat[key][subkey]
        self.iter_dicts.append(table)

    def _add_rows(self):
        if(not self.iter_dicts):
            self.add_partial_row()
        else:
            for d in self.iter_dicts:
                if(int(d['Hit_num']) <= self.args.max_hits):
                    for key, val in self.dat['stat'].items():
                        d[key] = val
                    col = tuple(sorted(d.keys()))
                    row = tuple(map(d.get, col))
                    try:
                        self.row_by_col[col].append(row)
                    except:
                        self.row_by_col[col] = [row]
        self.iter_dicts = []

    def clear_iter(self):
        '''
        Adds all data from current iteration to the database and frees the
        iteration and its children hits and hsps from memory
        '''
        self._add_rows()
        self.dat['iter'] = {}
        self.dat['stat'] = {}
        self.clear_hit()

    def clear_hit(self):
        '''
        Clears the current hit and all children hsps from memory
        '''
        self.dat['hit'] = {}
        self.clear_hsp()

    def clear_hsp(self):
        '''
        Clears hsp from memmory
        '''
        self.dat['hsp'] = {}

    def add(self, tag, text):
        '''
        Input: One tag and its text (possibly None)
        '''
        tag = re.sub('-', '_', tag)
        if(text is None or text.isspace()): pass
        elif('Hsp_' in tag):
            if(tag in ('Hsp_qseq', 'Hsp_hseq', 'Hsp_midline') and self.args.small):
                pass
            else:
                self.dat['hsp'][tag] = text
        elif('Hit_' in tag):
            self.dat['hit'][tag] = text
        elif('Iteration_' in tag):
            if(tag == 'Iteration_query_def'):
                self.dat['iter']['query_seqid'] = re.sub('(\S+).*', '\\1', text)
            self.dat['iter'][tag] = text
        elif('Statistics_' in tag):
            self.dat['stat'][tag] = text
        elif('BlastOutput_' in tag or 'Parameters_' in tag):
            if('reference' in tag or 'query' in tag):
                pass
            else:
                self.dat['root'][tag] = text
