#! /usr/bin/python3

import lib.sqlite_interface as misc
import traceback
import sys

# ==================
# EXPORTED FUNCTIONS
# ==================

def init_blastreport(cur, verbose=False):
    BLAST_VAL = """
        -- User Input
        collection TEXT COLLATE NOCASE,
        db_desc    TEXT COLLATE NOCASE,

        BlastOutput_program   TEXT NOT NULL COLLATE NOCASE,
        BlastOutput_version   TEXT NOT NULL COLLATE NOCASE,
        BlastOutput_db        TEXT NOT NULL COLLATE NOCASE,
        Parameters_matrix     TEXT NOT NULL COLLATE NOCASE,
        Parameters_expect     REAL NOT NULL,
        Parameters_gap_open   INTEGER NOT NULL,
        Parameters_gap_extend INTEGER NOT NULL,
        Parameters_filter     TEXT NOT NULL COLLATE NOCASE,

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

        -- Query data, parsed, if possible, from the fasta header
        Query_locus TEXT COLLATE NOCASE,
        Query_gi    TEXT NOT NULL,
        Query_gb    TEXT NOT NULL COLLATE NOCASE,
        Query_gene  TEXT COLLATE NOCASE,
        Query_taxon TEXT,

        -- Don't change this default ... please
        -- Hit_num and Hsp_num are two elements in the PK,
        -- ensuring they have default values gaurds against
        -- weird errors
        Hit_num       INTEGER DEFAULT 0,
        Hit_id        TEXT COLLATE NOCASE,
        Hit_def       TEXT COLLATE NOCASE,
        Hit_accession TEXT COLLATE NOCASE,
        Hit_len       INTEGER,

        -- Don't change Hsp_num default
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
        Hsp_qseq        TEXT COLLATE NOCASE,
        Hsp_hseq        TEXT COLLATE NOCASE,
        Hsp_midline     TEXT COLLATE NOCASE,

        CHECK(Parameters_expect >= 0),
        CHECK(Parameters_gap_open >= 0),
        CHECK(Parameters_gap_extend >= 0),
        CHECK(Iteration_iter_num >= 0),
        CHECK(Iteration_query_len >= 0),
        CHECK(Statistics_db_num >= 0),
        CHECK(Statistics_db_len >= 0),
        CHECK(Statistics_hsp_len >= 0),
        CHECK(Statistics_eff_space >= 0),
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

        PRIMARY KEY(blastoutput_db, query_gb, hit_num, hsp_num)
        """

    cmds = (
        "DROP TABLE IF EXISTS BlastReport",
        "CREATE TABLE BlastReport(" + BLAST_VAL + ")",
        "CREATE INDEX query_gene_idx ON BlastReport (query_gene)",
        "CREATE INDEX query_locus_idx ON BlastReport (query_locus)",
        "CREATE INDEX query_gi_idx ON BlastReport (query_gi)",
        "CREATE INDEX query_gb_idx ON BlastReport (query_gb)",
        "CREATE INDEX Iteration_iter_num_idx ON BlastReport (Iteration_iter_num)")
    create_table(cur, cmds)

def init_blastdatabase(cur, verbose=False):
    DATABASE_VAL = """
        database TEXT PRIMARY KEY,
        taxid    INTEGER CHECK(taxid >= 0),
        species  TEXT COLLATE NOCASE,
        alphabet TEXT COLLATE NOCASE,
        source   TEXT COLLATE NOCASE,
        comment  TEXT COLLATE NOCASE
    """

    cmds = (
        "DROP TABLE IF EXISTS BlastDatabase",
        "CREATE TABLE BlastDatabase(" + DATABASE_VAL + ")")
    create_table(cur, cmds)

def init_mrca(cur, verbose=False):
    MRCA_VAL = """
        taxid_1 INTEGER NOT NULL CHECK(taxid_1 >= 0),
        taxid_2 INTEGER NOT NULL CHECK(taxid_2 >= 0),
        mrca    INTEGER NOT NULL CHECK(mrca >= 0),
        phylostratum INTEGER NOT NULL CHECK(phylostratum >= 0),
        PRIMARY KEY(taxid_1, taxid_2)
    """

    cmds = (
        "DROP TABLE IF EXISTS MRCA",
        "CREATE TABLE MRCA(" + MRCA_VAL + ")")
    create_table(cur, cmds)

def init_taxid2name(cur, verbose=False):
    TAX2NAME_VAL = """
        taxid   INTEGER NOT NULL CHECK(taxid >= 0),
        sciname TEXT NOT NULL COLLATE NOCASE,
        PRIMARY KEY(taxid, sciname)
    """


    cmds = (
        "DROP TABLE IF EXISTS Taxid2Name",
        "CREATE TABLE Taxid2Name(" + TAX2NAME_VAL + ")")
    create_table(cur, cmds)

def init_besthits(cur, verbose=False):
    BESTHITS_VAL = \
    """
    database NOT NULL,

    qgene,
    qgb NOT NULL,
    qgi NOT NULL,
    qlocus,
    qtaxon INTEGER,
    qlen   INTEGER NOT NULL,

    -- Summed HSP values
    snhsp,

    shit   INTEGER DEFAULT 0,
    shlen  INTEGER DEFAULT 0,
    salen  INTEGER DEFAULT 0,
    sgaps  INTEGER DEFAULT 0,
    sident INTEGER DEFAULT 0,
    spos   INTEGER DEFAULT 0,
    sscore INTEGER DEFAULT 0,

    -- Max single HSP
    mnhsp,

    mevalue NUMERIC DEFAULT 999,

    mhit   INTEGER DEFAULT 0,
    mhlen  INTEGER DEFAULT 0,
    malen  INTEGER DEFAULT 0,
    mgaps  INTEGER DEFAULT 0,
    mident INTEGER DEFAULT 0,
    mpos   INTEGER DEFAULT 0,
    mscore INTEGER DEFAULT 0,

    -- Maximum non-overlapping path (neither query nor hit is allowed to
    -- overlap)
    pnhsp,

    phit   INTEGER DEFAULT 0,
    phlen  INTEGER DEFAULT 0,
    palen  INTEGER DEFAULT 0,
    pgaps  INTEGER DEFAULT 0,
    pident INTEGER DEFAULT 0,
    ppos   INTEGER DEFAULT 0,
    pscore INTEGER DEFAULT 0,

    CHECK (sident <= spos),
    CHECK (pident <= ppos),
    CHECK (mident <= mpos),

    CHECK (salen  >= 0),
    CHECK (sgaps  >= 0),
    CHECK (sident >= 0),
    CHECK (spos   >= 0),
    CHECK (sscore >= 0),

    CHECK (palen  >= 0),
    CHECK (pgaps  >= 0),
    CHECK (pident >= 0),
    CHECK (ppos   >= 0),
    CHECK (pscore >= 0),

    CHECK (malen  >= 0),
    CHECK (mgaps  >= 0),
    CHECK (mident >= 0),
    CHECK (mpos   >= 0),
    CHECK (mscore >= 0)

    PRIMARY KEY(database, qgb)
    """

    cmds = (
        "DROP TABLE IF EXISTS BestHits",
        "CREATE TABLE BestHits(" + BESTHITS_VAL + ")",
        "CREATE INDEX qgene_idx  ON BestHits (qgene )",
        "CREATE INDEX qlocus_idx ON BestHits (qlocus)",
        "CREATE INDEX qgi_idx    ON BestHits (qgi   )"
        )
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
