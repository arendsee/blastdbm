#! /usr/bin/python3

import sys
import argparse
import orphanlib.sqlite_interface  as misc
import orphanlib.blastin           as blastin
import orphanlib.initialize        as initialize
import orphanlib.meta              as meta
import orphanlib.query             as query
import orphanlib.dbtools           as tools

# Top parser
parser = argparse.ArgumentParser(
    description="Manages io for blast result database (FIX SUBCOMMAND HELP)",
    epilog='For help on any sub-command, enter: %(prog)s [sub-command] -h')

# Input parent parser
_input = argparse.ArgumentParser(add_help=False)
_input.add_argument(
    '-i', '--input',
    help="Input file (from stdin)",
    default=sys.stdin)

_sqldb = argparse.ArgumentParser(add_help=False)
_sqldb.add_argument(
    '-q', '--sqldb',
    help="SQL database name")

_csv = argparse.ArgumentParser(add_help=False)
_csv.add_argument(
    '--delimiter',
    help="CSV file delimiter (default=',')",
    default=',')

sub = parser.add_subparsers(
    help='sub-command help')

# orphanlib.blast parser
blastin.parse(sub, _input, _sqldb)

# orphanlib.query parser
query.parse(sub, _sqldb, **{'csv':_csv})

# orphanlib.dbtools parser
tools.parse(sub, _sqldb)

# Parse arguments
args = parser.parse_args()

if(not getattr(args, 'sqldb')):
    sys.stderr.write("Please provide an SQL filename (-q <filename>)\n")
    sys.exit()

with misc.open_db(args.sqldb) as con:
    # Pass arguments to proper sub-command
    cur = con.cursor()
    cur.execute('pragma shrink_memory')

    # This is kind of a dirty hack, but I need two cursors, one for select on
    # for insert in update_besthits. So I need to send in the con not just the
    # cur
    if(args.func.__name__ == 'update'):
        args.func(args, cur, con)
    else:
        args.func(args, cur)
