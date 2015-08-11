#!/usr/bin/env python3

import sys
import argparse
import lib.sqlite_interface  as misc
import lib.blastin           as blastin
import lib.initialize        as initialize
import lib.meta              as meta
import lib.query             as query
import lib.dbtools           as tools

__version__ = "0.1.1"

def parser(argv=None):
    # Top parser
    parser = argparse.ArgumentParser(
        description="Manages io for blast result database (FIX SUBCOMMAND HELP)",
        prog="fagin",
        epilog='For help on any sub-command, enter: %(prog)s [sub-command] -h')

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {}'.format(__version__)
    )

    # Input parent parser
    _input = argparse.ArgumentParser(add_help=False)
    _input.add_argument(
        '-i', '--input',
        help="Input file (from stdin)",
        type=argparse.FileType('r'),
        nargs="*",
        metavar="FILE"
    )

    _sqldb = argparse.ArgumentParser(add_help=False)
    _sqldb.add_argument(
        '-q', '--sqldb',
        help="SQL database name",
        metavar="FILE"
    )

    _csv = argparse.ArgumentParser(add_help=False)
    _csv.add_argument(
        '--delimiter',
        help="CSV file delimiter (default=',')",
        metavar="DEL",
        default=','
    )

    sub = parser.add_subparsers(
        help='sub-command help'
    )

    # lib.blast parser
    blastin.parse(sub, _input, _sqldb)

    # lib.query parser
    query.parse(sub, _sqldb, **{'csv':_csv})

    # lib.dbtools parser
    tools.parse(sub, _sqldb)

    # Parse arguments
    args = parser.parse_args()

    if(not getattr(args, 'sqldb')):
        sys.stderr.write("Please provide an SQL filename (-q <filename>)\n")
        sys.exit()

    return(args)


if __name__ == '__main__':
    args = parser()
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
