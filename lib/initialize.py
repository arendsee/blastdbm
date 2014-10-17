#! /usr/bin/python3

import lib.sqlite_interface as misc
import traceback
import sys

# ==================
# EXPORTED FUNCTIONS
# ==================

def init_hsp(cur, verbose=False):
    VAL = """
        db     TEXT NOT NULL COLLATE NOCASE,
        qseqid TEXT NOT NULL COLLATE NOCASE,
        hseqid TEXT NOT NULL COLLATE NOCASE,

        Hit_num       INTEGER DEFAULT 0,
        Hit_len       INTEGER,
        Hsp_num         INTEGER DEFAULT 0,
        Hsp_bit_score   REAL DEFAULT 0,
        Hsp_score       REAL DEFAULT 0,
        Hsp_evalue      REAL,
        Hsp_query_from  INTEGER,
        Hsp_query_to    INTEGER,
        Hsp_hit_from    INTEGER,
        Hsp_hit_to      INTEGER,
        Hsp_query_frame INTEGER,
        Hsp_hit_frame   INTEGER,
        Hsp_identity    INTEGER DEFAULT 0,
        Hsp_positive    INTEGER DEFAULT 0,
        Hsp_align_len   INTEGER DEFAULT 0,
        Hsp_gaps        INTEGER DEFAULT 0,

        PRIMARY KEY(db, qseqid, Hit_num, Hsp_num)
    """

    diagnostics = """
        CHECK(Hit_num >= 0),
        CHECK(Hit_len >= 0),
        CHECK(Hsp_num >= 0),
        CHECK(Hsp_bit_score >= 0),
        CHECK(Hsp_score >= 0),
        CHECK(Hsp_evalue >= 0),
        CHECK(Hsp_identity >= 0),
        CHECK(Hsp_positive >= 0),
        CHECK(Hsp_align_len >= 0),
        CHECK(Hsp_gaps >= 0),
        CHECK(Hsp_identity >= 0),
        CHECK(Hsp_positive >= 0),
        CHECK(Hsp_align_len >= 0),
        CHECK(Hsp_gaps >= 0)
    """

    cmds = (
        "DROP TABLE IF EXISTS hsp",
        "CREATE TABLE hsp(" + VAL + ")"
    )
    create_table(cur, cmds)

def init_iteration(cur, verbose=False):
    VAL = """
        db     TEXT COLLATE NOCASE,
        qseqid TEXT COLLATE NOCASE,

        Iteration_iter_num   INTEGER NOT NULL,
        Iteration_query_ID   TEXT NOT NULL COLLATE NOCASE,
        Iteration_query_def  TEXT COLLATE NOCASE,
        Iteration_query_len  INTEGER NOT NULL,
        Iteration_message    TEXT COLLATE NOCASE,

        Statistics_db_num    INTEGER,
        Statistics_db_len    INTEGER,
        Statistics_hsp_len   INTEGER,
        Statistics_eff_space REAL,
        Statistics_kappa     REAL,
        Statistics_lambda    REAL,
        Statistics_entropy   REAL,

        PRIMARY KEY(db, qseqid)
    """

    diagnostics = """
        CHECK(Iteration_iter_num >= 0),
        CHECK(Iteration_query_len >= 0),
        CHECK(Statistics_db_num >= 0),
        CHECK(Statistics_db_len >= 0),
        CHECK(Statistics_hsp_len >= 0),
        CHECK(Statistics_eff_space >= 0),
    """

    cmds = (
        "DROP TABLE IF EXISTS iteration",
        "CREATE TABLE iteration(" + VAL + ")")
    create_table(cur, cmds)

def init_database(cur, verbose=False):
    VAL = """
        db                    TEXT PRIMARY KEY,
        BlastOutput_program   TEXT NOT NULL COLLATE NOCASE,
        BlastOutput_version   TEXT NOT NULL COLLATE NOCASE,
        Parameters_matrix     TEXT NOT NULL COLLATE NOCASE,
        Parameters_expect     REAL NOT NULL,
        Parameters_gap_open   INTEGER NOT NULL,
        Parameters_gap_extend INTEGER NOT NULL,
        Parameters_filter     TEXT NOT NULL COLLATE NOCASE
    """

    diagnostics = """
        CHECK(Parameters_expect >= 0),
        CHECK(Parameters_gap_open >= 0),
        CHECK(Parameters_gap_extend >= 0),
    """

    cmds = (
        "DROP TABLE IF EXISTS database",
        "CREATE TABLE database(" + VAL + ")")
    create_table(cur, cmds)

# =================
# UTILITY FUNCTIONS
# =================

def create_table(cur, cmds):
    for cmd in cmds:
        try:
            cur.execute(cmd)
        except Exception as e:
            print("Error on command\n{}".format(cmd))
            print(e, file=sys.stderr)
            traceback.print_exc(file=sys.stderr)
