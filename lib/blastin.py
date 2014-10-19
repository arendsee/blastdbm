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

    database = Database(cur)
    iteration = Iteration(cur)
    hit = Hit()
    hsp = Hsp(cur)

    for event, elem in con:
        if event == "end":
            if elem.tag in hsp.fields:
                hsp.add(elem.tag, elem.text)
            elif elem.tag in hit.fields:
                hit.add(elem.tag, elem.text)
            elif elem.tag in iteration.fields:
                iteration.add(elem.tag, elem.text)
            elif elem.tag in database.fields:
                database.add(elem.tag, elem.text)
        elif event == "start":
            if elem.tag == "Hsp":
                hsp.new_entry(hit=hit, iteration=iteration)
            elif elem.tag == "Hit":
                hit.new_entry()
            elif elem.tag == "Iteration":
                iteration.new_entry(database=database)
            elif elem.tag == "BlastOutput":
                database.new_entry()


class Table:
    def __init__(self, cur, table, fields, commit_size=10000):
        self.cur = cur
        self.entry = dict()
        self.table = table
        self.commit_size = commit_size
        self.fields = fields
        self.ordered_fields = sorted(fields)
        self.colnames = [re.sub('-', '_', s) for s  in self.ordered_fields]
        self.entries = []

    def __del__(self):
        self._enter()
        self.commit()

    def add(self, key, val):
        try:
            self.entry[key] = val
        except KeyError:
            pass

    def _enter(self):
        if not self.entry:
            return
        try:
            self.entries.append([self.entry[k] for k in self.ordered_fields])
            if len(self.entries) == self.commit_size:
                self.commit()
        except KeyError:
            print("error in {}".format(self.table))
            print("\tFields: %s" % self.ordered_fields)
            print("\tFound: %s" % sorted(self.entry.keys()))

    def commit(self):
        self._enter()
        misc.insertmany(col=self.colnames,
                        rows=self.entries,
                        table=self.table,
                        cur=self.cur,
                        replace=True)
        self.entries = []

    def new_entry(self):
        self._enter()
        self.entry = dict()

class Hit():
    def __init__(self):
        self.entry = dict()
        self.fields = {
            'Hit_def',
            'sseqid',
            'Hit_num',
            'Hit_len'
        }

    def add(self, key, val):
        if key == 'Hit_def':
            self.entry['sseqid'] = re.match('\S+', val).group(0)
        self.entry[key] = val

    def new_entry(self):
        self.entry = dict()

class Hsp(Table):
    def __init__(self, cur):
        fields = {
            'db',
            'qseqid',
            'qlen',
            'hlen',
            'sseqid',
            'hit_num',
            'Hsp_num',
            'Hsp_bit-score',
            'Hsp_score',
            'Hsp_evalue',
            'Hsp_query-from',
            'Hsp_query-to',
            'Hsp_hit-from',
            'Hsp_hit-to',
            'Hsp_query-frame',
            'Hsp_hit-frame',
            'Hsp_identity',
            'Hsp_positive',
            'Hsp_align-len',
            'Hsp_gaps'
        }
        super().__init__(cur=cur, table='Hsp', fields=fields)
        self.colnames = [re.sub('Hsp_num', 'hsp_num', s) for s in self.colnames]
        self.colnames = [re.sub('Hsp_', '', s) for s in self.colnames]


    def new_entry(self, hit, iteration):
        if self.entry:
            self._enter()
        self.entry = {
            'db':iteration.entry['db'],
            'qseqid':iteration.entry['qseqid'],
            'qlen':iteration.entry['Iteration_query-len'],
            'hlen':hit.entry['Hit_len'],
            'sseqid':hit.entry['sseqid'],
            'hit_num':hit.entry['Hit_num']
        }

class Iteration(Table):
    def __init__(self, cur):
        fields = {
            'db',
            'qseqid',
            'Iteration_iter-num',
            'Iteration_query-ID',
            'Iteration_query-def',
            'Iteration_query-len',
            'Iteration_message',
            'Statistics_db-num',
            'Statistics_db-len',
            'Statistics_hsp-len',
            'Statistics_eff-space',
            'Statistics_kappa',
            'Statistics_lambda',
            'Statistics_entropy'
        }
        super().__init__(cur=cur, table='Iteration', fields=fields)

    def add(self, key, val):
        if key == 'Iteration_query-def':
            self.entry['qseqid'] = re.match('\S+', val).group(0)
        self.entry[key] = val

    def new_entry(self, database):
        if self.entry:
            self._enter()
        self.entry = {
            'db':database.entry['BlastOutput_db'],
            'Iteration_message':''
        }

class Database(Table):

    def __init__(self, cur):
        fields = {
            'BlastOutput_db',
            'BlastOutput_program',
            'BlastOutput_version',
            'Parameters_matrix',
            'Parameters_expect',
            'Parameters_gap-open',
            'Parameters_gap-extend',
            'Parameters_filter'
        }
        super().__init__(cur=cur, table='Database', fields=fields)
        self.colnames = sorted([re.sub('BlastOutput_db', 'db', s) for s in self.colnames])

    def add(self, key, val):
        if key == 'BlastOutput_db':
            self.entry[key] = re.sub('\..*$', '', os.path.basename(val))
        else:
            self.entry[key] = val
