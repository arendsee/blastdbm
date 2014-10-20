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
        sseqid TEXT NOT NULL COLLATE NOCASE,
        qlen INTEGER DEFAULT 0,

        hit_num     INTEGER DEFAULT 0,
        hlen        INTEGER,
        hsp_num     INTEGER DEFAULT 0,
        bit_score   REAL DEFAULT 0,
        score       REAL DEFAULT 0,
        evalue      REAL,
        query_from  INTEGER,
        query_to    INTEGER,
        hit_from    INTEGER,
        hit_to      INTEGER,
        query_frame INTEGER,
        hit_frame   INTEGER,
        identity    INTEGER DEFAULT 0,
        positive    INTEGER DEFAULT 0,
        align_len   INTEGER DEFAULT 0,
        gaps        INTEGER DEFAULT 0,

        PRIMARY KEY(db, qseqid, hit_num, hsp_num)
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
        "CREATE TABLE hsp(" + VAL + ")",
        "CREATE INDEX hsp_qseqid ON hsp(qseqid)",
        "CREATE INDEX hsp_db ON hsp(db)",
        "CREATE INDEX hsp_sseqid ON hsp(sseqid)"
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
        "CREATE TABLE iteration(" + VAL + ")",
        "CREATE INDEX iteration_qseqid on iteration(qseqid)"
    )
    create_table(cur, cmds)

def init_database(cur, verbose=False):
    VAL = """
        db                    TEXT,
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
