#! /usr/bin/python3

import argparse
import sys
import csv
import orphanlib.sqlite_interface as misc
import re
import numpy
import pandas

# =========
# CONSTANTS
# =========

MISSING_DATA = 'MISSING'

# ==================
# EXPORTED FUNCTIONS
# ==================

def parse(parent, *args, **kwargs):
    """ Handles arguments for the 'query' subcommand
    """

    parser = parent.add_parser(
        'query',
        help="Retrieve data from database")

    # Parent parser for sub-commands that take an identifier such as
    # locus, gi or gb
    _identifiers = argparse.ArgumentParser(add_help=False)
    _identifiers.add_argument(
        '-e', '--identifier',
        help="Identifier choice",
        choices=["gb", "gi", "locus"],
        default="locus")
    _identifiers.add_argument(
        '-f', '--from_file',
        help="Read identifiers from file")
    _identifiers.add_argument(
        '-j', '--from_list',
        help="Read from here",
        nargs="+")
    _identifiers.add_argument(
        '-c', '--col',
        help="Draw from this blast collection")
    _identifiers.add_argument(
        '-a', '--all',
        help="Search entire database",
        action='store_true', default=False)

    # Parent criterion class
    _criterion = argparse.ArgumentParser(add_help=False)
    _criterion.add_argument(
        '-r', '--criterion',
        help="Criterion by which clade specifivity is determined",
        nargs=2,
        default=['score', 100])

    # Subparser for all sub-commands
    sub = parser.add_subparsers(
        help="add_subparsers STUB",
        dest="query_function")

    # Submit a raw SQL command
    raw = sub.add_parser(
        'raw',
        help="Submit raw SQLite command",
        parents=(args + (kwargs['csv'], )),
        description="raw description stub")
    raw.add_argument(
        'sqlcmd',
        help="raw sql command (this function should be removed)")

    # Retrieve a score matrix
    mat = sub.add_parser(
        'mat',
        help="Fetch a matrix",
        parents=(args + (_identifiers, )))
    mat.add_argument(
        '-l', '--filling',
        help="Fill matrix with this hsp value",
        default="hsp_bit_score")
    mat.add_argument(
        '-o', '--output',
        help="Output csv file",
        default=sys.stdout)

    # Retrieve the clade specificity of a single locus
    spec = sub.add_parser(
        'spec',
        help="Get the clade specificity of a single sequence",
        parents=(args + (_identifiers, _criterion, )))

    # Retrieve phylostrata info for a list of loci
    phylo = sub.add_parser(
        'phylo',
        help="Get phylostrata info",
        parents=(args + (_identifiers, _criterion, kwargs['csv'])))
    phylo.add_argument(
        '-p', '--pathway',
        help="Pathway or network name",
        default="Nameless")
    phylo.add_argument(
        '-t', '--outfmt',
        help="Output format",
        default='csv')

    maxattr = sub.add_parser(
        'maxattr',
        help="Get attributes of maximum scoring database/query pairs",
        parents=(args + (kwargs['csv'], )))
    maxattr.add_argument(
        '--fields',
        nargs='+')
    maxattr.add_argument(
        '--no-basic-fields',
        help="Don't use default set of generally useful fields",
        dest='basic_fields',
        action='store_false',
        default=True)
    maxattr.add_argument(
        '--condition',
        help="Condition upon which to keep")

    parser.set_defaults(func=_dispatch)


# =================
# UTILITY FUNCTIONS
# =================

def _dispatch(args, cur):
    call = {'raw': _fetch_and_print,
            'mat': _get_mat,
            'phylo': _phylo,
            'spec': _spec,
            'maxattr': _get_maxattr}
    call[args.query_function](args, cur)

def _phylo(args, cur):
    qdat = _get_query_data(args, cur)
    if(args.outfmt == 'json'):
        _phylo_json(qdat, args.pathway)
    elif(args.outfmt == 'csv'):
        _phylo_csv(qdat, args.delimiter)
    else:
        print("Unrecognized format '{}', dying".format(args.outfmt), file=sys.stderr)

def _fetch_and_print(args, cur):
    rows = misc.fetch(args.sqlcmd, cur)
    for row in rows:
        print((args.delimiter).join(map(str, row)))

def _get_mat(args, cur):
    args.identifier = 'query_{}'.format(args.identifier)
    queries = _get_identifiers(args, cur)

    # Every element (locus or gb) must be a str (not a tuple, as can easily
    # happen)
    for x in queries:
        assert isinstance(x, str), "Found {} expected string".format(str(type(x)))

    dbs = misc.get_db(cur)

    mat = pandas.DataFrame(index=queries, columns=dbs)

    if(args.filling == 'hsp_evalue'):
        act = 'min'
    else:
        act = 'max'
    d = {'d':'blastoutput_db',
         'a':act,
         'f':args.filling,
         't':'blastreport',
         'i':args.identifier,
         'c':"and collection = '{}'".format(args.col) if args.col else "",
         'v':'{}'}
    cmd = "select {d}, {a}({f}) from {t} where {i} = '{v}' {c} group by {d}".format(**d)

    for q in queries:
        vals = misc.fetch(cmd.format(q), cur)
        # A return value of none implies an iteration with no hits, this
        # corresponds to a high evalue or a 0 score
        if(args.filling == 'hsp_evalue'):
            vals = [(x[0], 99) if x[1] is None else x for x in vals]
        else:
            vals = [(x[0], 0) if x[1] is None else x for x in vals]
        for pair in vals:
            mat.at[q, pair[0]] = pair[1]

    # Drop empty columns (databases outside collection)
    mat = mat.dropna(axis=1, how='all')

    # Write output csv file
    mat.to_csv(args.output)

def _phylo_json(qdat, pathway):
    # Write output in JSON format
    out = '{\n"Pathway":"' + pathway + '",\n'
    out += '"Protein":\n\t[\n\t\t{\n'
    gb_dat = []
    for q in qdat:
        d = []
        for field, value in q.items():
            d.append('\t\t\t"{}":"{}"'.format(field, value))
        gb_dat.append(',\n'.join(d) + "\n\t\t}")
    out += ",\n\t\t{\n".join(gb_dat) + "\n"
    out += "\t]\n}"
    # Remove sinlge quotes
    out = re.sub("'", "", out)
    print(out)

def _phylo_csv(qdat, delimiter):
    fields = qdat[0].keys()
    print('{}'.format(delimiter).join(fields))
    for q in qdat:
        print('{}'.format(delimiter).join([str(q[y]) for y in fields]))

def _read_single_column(filename):
    with open(filename) as f:
        dat = f.readlines()
    dat = [re.sub("^\s*['\"]?\s*|\s*['\"]?\s*,?\s*$", "", x) for x in dat]
    return(dat)

def _get_identifiers(args, cur):
    q = []
    if(args.from_list):
        q += args.from_list
    if(args.from_file):
        q += _read_single_column(args.from_file)
    if(args.col):
        q += misc.get_fields(args.identifier, 'blastreport', cur,
                             ident='collection', value=args.col, is_distinct=True)
    return(q)

def _get_query_data(args, cur):
    ids = _get_identifiers(args, cur)
    qdat = []
    ps = _spec(args.identifier, ids, args.criterion, cur)
    for i in ids:
        q = misc.get_query_info(args.identifier, i, cur)[0]
        if(args.criterion):
            q['phylostratum'] = ps[i]
        qdat.append(q)
    return(qdat)

def _get_maxattr(args, cur):
    fields = set()
    if(args.fields):
        fields = args.fields
    if(args.basic_fields):
        fields.update(['iteration_query_len', 'hit_len',
                       'hit_def',
                       'hsp_query_from', 'hsp_query_to',
                       'hsp_hit_from', 'hsp_hit_to',
                       'hsp_positive', 'hsp_bit_score',
                       'mrca.phylostratum', 'mrca.mrca',
                       'taxid2name.sciname'])
    raw_results = misc.get_maxattr(fields, cur, condition=args.condition)
    writer = csv.writer(sys.stdout, delimiter=args.delimiter,
                        quoting=csv.QUOTE_MINIMAL)
    writer.writerow(['database','query'] + list(fields))
    for row in raw_results:
        writer.writerow(row[:2] + row[3:])

def _spec(identifier, values, criterion, cur):
    out = {}
    for val in values:
        out[val] = misc.spec(identifier, val, criterion, cur)
    return(out)
