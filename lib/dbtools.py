#! /usr/bin/python3

import argparse
import csv
import sys

import lib.meta as meta
import lib.initialize as init
import lib.sqlite_interface as misc

def parse(parent, *args, **kwargs):
    # Add blast database info to sql database
    db_parser = parent.add_parser(
        'update',
        parents=args,
        help="Update blast database info")
    db_parser.add_argument(
        '-t', '--taxids',
        help="Space delimited list of taxids to add to MRCA table",
        type=int,
        nargs='+')
    db_parser.add_argument(
        '-d', '--deep',
        help="Search blastreport for database names (usually unnecessary)",
        action='store_true',
        default=False)
    db_parser.add_argument(
        '--destroy',
        help="Destroy current blastdatabase file (does not lose blast results)",
        action='store_true',
        default=False)
    db_parser.set_defaults(func=update)

    dump_parser = parent.add_parser(
        'dump',
        parents=args,
        help="Dump tables in various ways (currently only --mrca is implemented)")
    dump_parser.add_argument(
        '--mrca',
        help="Dumps MRCA data relative to given taxon id")
    dump_parser.set_defaults(func=dump_dbinfo)

def update(args, cur, con=None):
    meta.update_dbinfo(cur, deep=args.deep, destroy=args.destroy)
    meta.update_mrca(cur, sync=True, taxids=args.taxids)
    meta.update_besthits(cur, con)

def dump_dbinfo(args, cur):
    if args.mrca:
        dfields = ('database', 'species', 'taxid')
        dtable = 'blastdatabase'
        res1 = misc.get_fields(dfields, dtable, cur)
        mfields = ('mrca', 'phylostratum')
        cmd = 'select {} from mrca where taxid_1 = {} and taxid_2 = {}'.format(
                ', '.join(mfields), args.mrca, '{}')
        c = csv.writer(sys.stdout)
        c.writerow(('Database', 'Species', 'Taxid', 'MRCA_Taxid', 'Stratum', 'MRCA_Name'))
        for line in res1:
            try:
                mrca = misc.fetch(cmd.format(line[2]), cur)[0]
            except IndexError:
                sys.exit("Invalid taxonid (must be a leaf taxon)")
            mrca_cmd = "select sciname from taxid2name where taxid = {}"
            mrca_name = misc.fetch(mrca_cmd.format(mrca[0]), cur)[0][0]
            c.writerow(line + mrca + (mrca_name, ))
