#! /usr/bin/python3

import sqlite3 as sql
import sys
import argparse
import traceback

# ==============
# BASE FUNCTIONS
# ==============

def insert(dic, table, cur, replace=False):
    row = tuple(dic.values())
    col = tuple(dic.keys())
    rep = "OR REPLACE" if replace else ""
    cmd = ' '.join(("INSERT", rep, "INTO", table, "(",
                      ', '.join(col),
                      ") VALUES (",
                      ', '.join("?"*len(col)),
                      ")"))
    try:
        # ONLY values can be substituted for '?', NOT table or column names
        cur.execute(cmd, row)
    except Exception as e:
        _sql_err(e, cmd)

def insertmany(col, rows, table, cur, replace=False):
    rep = 'OR REPLACE' if replace else ''
    cmd = ' '.join(("INSERT", rep, "INTO", table, "(",
                      ', '.join(col),
                      ") VALUES (",
                      ', '.join("?"*len(col)),
                      ")"))
    try:
        # ONLY values can be substituted for '?', NOT table or column names
        cur.executemany(cmd, rows)
    except Exception as e:
        print(rows)
        _sql_err(e, cmd)

def update(dic, table, cond, cur):
    sets = []
    if(isinstance(cond[1], str)):
        condition = "where {} = '{}'".format(*cond)
    else:
        condition = "where {} in ({})".format(cond[0], _quote(', '.join(cond[1])))
    for col, val in dic.items():
        val = _quote(val)
        sets.append("{} = {}".format(col, val))
    cmd = ' '.join(("UPDATE {} SET".format(table),
                   ', '.join(sets), ' ', condition))
    try:
        cur.execute(cmd)
    except Exception as e:
        _sql_err(e, cmd)

def fetch(cmd, cur):
    try:
        cur.execute(cmd)
        result = cur.fetchall()
    except Exception as e:
        _sql_err(e, cmd)
    return(result)


# =====================
# SQL COMMAND FUNCTIONS
# =====================

def open_db(filename):
    try:
        con = sql.connect(filename)
        return(con)
    except Exception as e:
        print("Error opening {}: {}".format(filename, e), file=sys.stderr)
        sys.exit(1)

def get_query_info(ident, value, cur):
    value = _quote(value)
    ident = _ident2field(ident)
    fields = (('query_locus','locus'),
              ('query_taxon','taxon'),
              ('query_gi','gi'),
              ('query_gb','gb'),
              ('iteration_query_len','length'))
    cmd = "select distinct {} from blastreport where {} = {}".format(
            ','.join([x[0] for x in fields]), ident, value)
    results = fetch(cmd, cur)
    out = [dict(zip([f[1] for f in fields], q)) for q in results]
    return(out)

def get_fields(fields, table, cur, ident=None, value=None, is_distinct=False):
    if(value):
        value = _quote(value)
        if(isinstance(value, str)):
            condition = "where {} = {}".format(ident, value)
        else:
            condition = ' '.join(("where {} in (", ''.join(value), ")")).format(ident)
    else:
        condition = ''

    result = _prepare_simple_select(fields, table, condition, cur,
                                    is_distinct=is_distinct)
    return(result)

def get_db(cur):
    result = _prepare_simple_select('database', 'blastdatabase', '', cur,
                                    is_distinct=True)
    return(result)

def spec(ident, value, criterion, cur):
    value = _quote(value)
    if(criterion[0] == 'evalue'):
        condition = 'hsp_evalue <=' + str(criterion[1])
    else:
        condition = 'hsp_bit_score >=' + str(criterion[1])

    column = _ident2field(ident)
    cmd = """
    select min(phylostratum) from mrca where
        taxid_1 =
        (
            select distinct query_taxon from blastreport
                where {0} = {1}
        )
    and
        taxid_2 in
        (
            select distinct taxid from blastdatabase where database in
            (
                select distinct blastoutput_db from blastreport
                    where {0} = {1} and {2}
            )
        )
    ;""".format(column, value, condition)

    out = fetch(cmd, cur)[0][0]

    # If no result, check to ensure the input identifier exists
    # If input exists, then the protein is specific to the input taxa,
    # so find the highest possible phylostratum
    if(out is None):
        cmd = ''.join(("select distinct ", column,
                       " from blastreport where ",
                       column, " = ", value))
        check = fetch(cmd, cur)
        if(0 == len(check)):
            print("Column {} does not contain value {}".format(column, value),
                  " dying painfully...",
                  file=sys.stderr)
            sys.exit(1)
        else:
            cmd = ' '.join((
                "select max(phylostratum) from mrca where taxid_1 = ",
                "(select distinct query_taxon from blastreport",
                "where {0} = {1})")).format(column, value)
            species_ps = fetch(cmd, cur)[0][0]
            out = species_ps
    return(out)

def table_exists(table, cur):
    cmd = """SELECT name
             FROM sqlite_master
             WHERE type='table'
             AND name='{}' COLLATE NOCASE""".format(table)
    result = fetch(cmd, cur)
    table_exists = True if len(result) > 0 else False
    return(table_exists)

def entry_exists(table, field, value, cur, condition=None):
    value = _quote(value)
    cmd = "select {0} from {1} where {0} = {2}".format(field, table, value)
    result = fetch(cmd, cur)
    entry_exists = True if len(result) > 0 else False
    return(entry_exists)

def get_maxattr(fields, cur, condition=None):
    # Retrieve all input fields where score is maximum for database/query pair
    condition = "where {}".format(condition) if condition else ''
    cmd = ''' \
            select blastdatabase.species, query_locus, max(hsp_bit_score), {}
            from (
                    ((blastreport inner join blastdatabase on
                        blastreport.blastoutput_db=blastdatabase.database)
                    inner join mrca on
                        blastreport.query_taxon=mrca.taxid_1 and
                        blastdatabase.taxid=mrca.taxid_2)
                    inner join taxid2name on
                        mrca.mrca=taxid2name.taxid
                 )
            {}
            group by blastoutput_db, query_locus
            order by query_locus, mrca.phylostratum, mrca.mrca
        '''.format(', '.join(fields), condition)
    return(fetch(cmd, cur))


# =================
# UTILITY FUNCTIONS
# =================

def _prepare_simple_select(fields, table, condition, cur, is_distinct=False):
    fields = fields if isinstance(fields, (tuple, list, set)) else (fields,)
    field_str = ', '.join(fields)
    dis_str = ' distinct ' if is_distinct else ''

    cmd = "select {} {} from {} {}".format(dis_str, field_str, table, condition)

    result = fetch(cmd, cur)

    result = _unnest(result)

    if(is_distinct):
        result = set(result)

    return(result)

def _ident2field(ident):
    if(ident in ('gb', 'gi', 'locus')):
        ident = "query_{}".format(ident)
    return(ident)

def _unnest(dat):
    out = []
    if(isinstance(dat, (tuple, list))):
        for el in dat:
            if(isinstance(el, (tuple, list)) and len(el) > 1):
                return(dat)
            out.append(el[0])
    else:
        return(dat)
    return(out)

def _quote(x):
    if(not isinstance(x, (list, tuple))):
        if(x is None):
            return 'NULL'
        else:
            return(''.join(("'", x, "'")))
    else:
        return([_quote(y) for y in x])

def _sql_err(e, cmd):
    print("Error on command\n{}".format(cmd))
    print(e, file=sys.stderr)
    traceback.print_exc(file=sys.stderr)
    sys.exit(1)


