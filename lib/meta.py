#! /usr/bin/python3

import re
import sys
from collections import defaultdict
import traceback

import lib.sqlite_interface as misc
import lib.entrez_interface as entrez
import lib.initialize as initialize
from   lib.lineage import Lineage, MRCA

# ==================
# EXPORTED FUNCTIONS
# ==================

def update_dbinfo(cur, deep=False, destroy=False, verbose=False):
    if(not misc.table_exists('blastdatabase', cur) or destroy):
        initialize.init_blastdatabase(cur, verbose)
        deep = True
    if(deep):
        print("Retrieving database names from blastreport (this may take awhile)")
        db = misc.get_fields('Blastoutput_db', 'blastreport', cur, is_distinct=True)
        for d in db:
            misc.insert({'database':d}, 'blastdatabase', cur)


    dbfiles = misc.get_fields('database', 'blastdatabase', cur, is_distinct=True)

    for f in dbfiles:
        ori = misc.get_fields(('species','taxid'), 'blastdatabase',
                               cur, ident='database', value=f)[0]

        if(ori[1] is not None): continue

        if(ori[0] is None):
            name = re.sub("\..*", "", f)
            name = re.sub("_", " ", name)
        else:
            name = ori[0]

        taxid = entrez.sciname2taxid(name)

        if(taxid is None):
            print("No taxid found for '{}' (taxon parsed as '{}')".format(f, name))
            name, taxid = _void_taxid(name)

        misc.update({'species': name, 'taxid': taxid}, 'blastdatabase',
                    ('database', f), cur)

def update_mrca(cur, sync=True, taxids=None, verbose=False):
    if(not misc.table_exists('mrca', cur)):
        initialize.init_mrca(cur, verbose)
    if(not misc.table_exists('Taxid2Name', cur)):
        initialize.init_taxid2name(cur, verbose)

    taxid_in = set()
    if(taxids):
        taxid_in.update(taxids)
    if(sync):
        db_taxids = misc.get_fields('taxid', 'blastdatabase', cur, is_distinct=True)
        taxid_in.update(db_taxids)

    # Retrieve taxonomy xml file from entrez
    lin = entrez.taxid2lineage(taxid_in)

    # Find mrca
    mrca = {} # A dict that holds mrca for pairs of taxids
    for k1 in lin:
        for k2 in lin:
            mrca[(k1, k2)] = MRCA(k1, k2)

    # Update MRCA SQL database |taxid1|taxid2|mrca|phylostratum|
    for key in mrca:
        dic = {'taxid_1':key[0].taxid,
               'taxid_2':key[1].taxid,
               'mrca':mrca[key].taxid,
               'phylostratum':mrca[key].level}
        misc.insert(dic, 'MRCA', cur, replace=True)

    # Update Taxid2Name SQL database |taxid|name|
    def t2n_insert(t, s):
        dic = {'sciname': s, 'taxid': t}
        misc.insert(dic, 'Taxid2Name', cur, replace=True)

    for l in lin:
        t2n_insert(l.taxid, l.sciname)
        for pair in l.lineage:
            t2n_insert(*pair)

def update_besthits(cur, con, verbose=False):
    if(not misc.table_exists('besthits', cur)):
        initialize.init_besthits(cur, verbose)
    colmap = (
        ('blastoutput_db', 'database'),
        ('query_seqid', 'qseqid'),
        ('query_gb', 'qgb'),
        ('Iteration_query_len', 'qlen'),
        ('query_gi', 'qgi'),
        ('query_locus', 'qlocus'),
        ('query_taxon', 'qtaxon'),
        ('query_gene', 'qgene'),
        ('hit_num', 'hit'),
        ('hit_len', 'hlen'),
        ('hsp_evalue', 'mevalue'),
        ('hsp_num', 'hsp'),
        ('hsp_query_from', 'qfrom'),
        ('hsp_query_to', 'qto'),
        ('hsp_hit_from', 'hfrom'),
        ('hsp_hit_to', 'hto'),
        ('hsp_align_len', 'alen'),
        ('hsp_gaps', 'gaps'),
        ('hsp_positive', 'pos'),
        ('hsp_identity', 'ident'),
        ('hsp_bit_score', 'score')
    )
    selection = ', '.join(["{} as {}".format(x,y) for x,y in colmap])
    cmd = """ \
        select
            {}
        from
            blastreport
        order by
            database,
            qseqid
        """.format(selection)
    cur.execute(cmd)

    scol = [x[1] for x in colmap]

    # There is a slight problem with the way I've been doing things Having my
    # only tie to the database be a cursor object does not allow me to read and
    # write simultaneously (I need to create a new cursor).
    writecur = con.cursor()
    scorers = (PathScorer(), SumScorer(), MaxScorer())
    for col, rows in block_of_rows_generator(cur, scol, scorers):
        misc.insertmany(col, rows, 'besthits', writecur, replace=True)


# =================
# UTILITY FUNCTIONS
# =================

def _set_taxid(name=None):
    qstr = "Please enter taxid (e.g. 3702): "
    taxid = input(qstr)
    try:
        n = entrez.taxid2sciname(taxid)
        if(taxid is None):
            print("Entrez does not recognize this taxon id")
            _set_taxid(name)
        qstr = "You mean this database is comprised solely of '{}' (y/n)? "
        r = input(qstr.format(n)).lower()
        if('n' in r):
            qstr = "Would you like to leave the taxid null (y/n)? "
            r = input(qstr).lower()
            if('y' in r):
                return(name, None)
            else:
                return(_set_taxid())
        qstr = "Would you like to set the taxon name to '{}': "
        r = input(qstr.format(n)).lower()
        if('y' in r):
            name = n
    except:
        print("Error...")
        return(_set_taxid())
    return(name, taxid)

def _void_taxid(name):
    taxid = None
    print("Entrez does not recognize taxon {}".format(name))
    r = input("Do you want to reset the name and try again (y/n)? ").lower()
    if('y' in r):
        name = input("Then give me a new taxon name: ")
        taxid = entrez.sciname2taxid(name)
        if(taxid is None):
            name, taxid = _void_taxid(name)
    else:
        r = input("Would you like to manually assign a taxon id (y/n)? ").lower()
        if('y' in r):
            name, taxid = _set_taxid()
        else:
            print("The taxid and name will be left NULL in the SQL database")
            name, taxid = (None, None)
    return((name, taxid))


# ===============
# BESTHIT CLASSES
# ===============

def block_of_rows_generator(cur, col, scorers, blocksize=1000):
    """\
    Generates tuples containing a tuple of column names and a tuple of rows
    """
    def to_col_and_row(rows_):
        dbcol = sorted(rows_[0].keys())
        dbrows = [tuple(r[k] for k in dbcol) for r in rows_]
        return(tuple((dbcol, dbrows)))

    rows = []
    for row in row_generator(cur, col, scorers):
        if(len(rows) > blocksize):
            yield to_col_and_row(rows)
            rows = []
        rows.append(row)
    yield to_col_and_row(rows)

def row_generator(cur, col, scorers):
    """\
    Generates dicts containing all sql database columns as keys and values as
    values
    """
    for pair in pair_generator(cur, col):
        out = {'database': pair.db, 'qseqid': pair.qseqid}
        out['mevalue'] = pair.get_evalue()
        for k,v in pair.qdat.items():
            out[k] = v
        for scorer in scorers:
            hit, s = pair.bestscore(scorer)
            prefix = scorer.name[0]
            out[prefix + 'hit'] = hit.hit
            out[prefix + 'hlen'] = hit.hlen
            s.add_prefixed_dat(out, prefix)
        yield out

def pair_generator(cur, col):
    """ \
    Yields a stream of query/database pairs. Each yielded value is comprised of
    a list of hits.
    """
    pair = DQPair()
    for hit in hit_generator(cur, col):
        if(not pair.add(hit)):
            yield pair
            pair = DQPair(hit)
    yield pair

def hit_generator(cur, col):
    """ \
    Yields a stream of hits, each including one or more hsp rows (some of which
    may be filled with null or 0 values for query/database pairs with no
    matches).
    """
    hit = Hit()
    for hsp in hsp_generator(cur, col):
        # If addition of the new hsp fails, yield the current hit and add the
        # new hsp to a new hit
        if(not hit.add(hsp)):
            yield hit
            hit = Hit(hsp)
    yield hit

def hsp_generator(cur, col, chunksize=1000):
    """ \
    Single rows of data from the sql database. $chunksize elements are accessed
    at a given time to increase efficiency.
    """
    while(True):
        rows = cur.fetchmany(chunksize)
        if(not rows):
            break
        for row in rows:
            yield Hsp(row, col)


class _Scorer:
    def __init__(self, name):
        self.name = name

    def score(self, hit):
        if(hit.len() == 1):
            hsp = hit.col[0]
            score = Score(hsp)
        else:
            score = self._calculate_score(hit.col)
        return(score)

    # TODO this should be an abstract method
    def _calculate_score(hit):
        return

class PathScorer(_Scorer):
    def __init__(self, cutoff=40):
        super().__init__('path')
        self.cutoff = cutoff

    def _calculate_score(self, hit):
        """\
        The _global_val function is very slow for large numbers of hsps, so the
        approximating function, _avg_val, will be called when size is above the
        given cutoff
        """
        if(len(hit) <= self.cutoff):
            out = self._global_val(hit)
        else:
            out = self._avg_val(hit)
        return(out)

    def _avg_val(self, hit):
        """ \
        An approximated score is produced by multiplying the total alignment length
        by the average score per aligned character Other values (e.g. positive,
        identical, gaps) are dealt with similarly
        """
        # This sorting step is essential to the algorithm
        # that calculates the non-overlapping alignment length
        hit = sorted(hit, key=lambda x: x.qfrom)

        # Number of query characters that are aligned against at least
        # one HSPs
        alen = 0

        # Beginning of query alignment
        qfrom = hit[0].qfrom

        # End of query alignment
        qto = hit[0].qto

        # The true alignment length is found be combining all overlapping
        # segments and then summing their lengths
        sum_fields = defaultdict(int)
        for i in range(1, len(hit)):
            hsp = hit[i]
            sum_fields['gaps'] += hsp.gaps
            sum_fields['ident'] += hsp.ident
            sum_fields['pos'] += hsp.pos
            sum_fields['score'] += hsp.score
            if(hit[i].qfrom <= qto):
                qto = max(hsp.qto, qto)
            else:
                alen += qto - qfrom + 1
                qfrom, qto = hsp.qfrom, hsp.qto
        alen += qto - qfrom + 1

        tlen = sum([x.alen for x in hit])

        avg_fields = {k:((v / tlen) * alen) for k,v in sum_fields.items()}
        avg_fields['alen'] = alen
        avg_fields['nhsp'] = len(hit)
        score = Score.fromdict(avg_fields)

        return(score)

    def _bestpath(self, v, end=-1, s=0, p=[]):
        """\
        Search all possible paths through the directed acyclic graph of
        non-overlapping subsequences (HSPs) with weights equal to the bitscores. I
        will use a very crude bruteforce algorithm. The time could be improved with
        a bit of branch pruning. It also might be worthwhile to implement this is C
        and optimize the hell out of it.  More importantly, this algorithm ensures
        a path with no overlaps in query sequence, however it still allows overlaps
        in hits. This should be remedied.
        """
        scores = []
        eol = True
        for i in range(len(v)):
            hsp = v[i]
            nv = v[0:i] + v[(i+1):]
            if(hsp.qfrom > end and not Hsp.contains_overlap(hsp, p)):
                eol = False
                scores.append(self._bestpath(
                    nv,
                    hsp.qto,
                    s + hsp.score,
                    p + [hsp]))
        if(eol):
            return((s,p))
        else:
            return(max(scores, key=lambda x: x[0]))

    def _global_val(self, hit):
        """\
        If there is only one HSP, simply return the single HSPs relevant values. In
        biological cases, the number of HSPs very often is 0 or 1.
        """
        path = self._bestpath(hit)
        out = {}
        out['gaps'] = sum(hsp.gaps for hsp in path[1])
        out['ident'] = sum(hsp.ident for hsp in path[1])
        out['pos'] = sum(hsp.pos for hsp in path[1])
        out['score'] = sum(hsp.score for hsp in path[1])
        out['alen'] = sum(hsp.alen for hsp in path[1])
        out['nhsp'] = len(path[1])
        score = Score.fromdict(out)
        return(score)

class SumScorer(_Scorer):
    def __init__(self):
        super().__init__('sum')

    def _calculate_score(self, hit):
        """ \
        Straight sum of all hsps summable values
        """
        out = {}
        out['gaps'] = sum(hsp.gaps for hsp in hit)
        out['ident'] = sum(hsp.ident for hsp in hit)
        out['pos'] = sum(hsp.pos for hsp in hit)
        out['score'] = sum(hsp.score for hsp in hit)
        out['alen'] = sum(hsp.alen for hsp in hit)
        out['nhsp'] = len(hit)
        score = Score.fromdict(out)
        return(score)

class MaxScorer(_Scorer):
    def __init__(self):
        super().__init__('max')

    def _calculate_score(self, hit):
        """ \
        Straight max of all hsps summable values
        """
        bestscore = -1
        maxhsp = None
        for hsp in hit:
            if(hsp.score > bestscore):
                maxhsp = hsp
                bestscore = hsp.score
        score = Score(maxhsp)
        return(score)

class _Collection:
    def __init__(self, init=None):
        self.col = []
        self.db = None
        self.qseqid = None
        self.qdat = {}
        if(init):
            self._first_item(init)

    def _first_item(self, item):
        self.col = [item]
        self.db = item.db
        self.qseqid = item.qseqid
        self.qdat = item.qdat

    def next(self):
        for c in self.col:
            yield c

    def add(self, item):
        if(not self.col):
            self._first_item(item)
            return(True)
        elif(self.col and type(self.col[0]).same_group(self.col[0], item)):
            self.col.append(item)
            return(True)
        else:
            return(False)

    def len(self):
        return(len(self.col))

class DQPair(_Collection):
    def bestscore(self, scorer, by='score'):
        """\
        Returns the pair's best hit's first hsp and score selected based on the
        criterion by (score by default)
        """
        best = (self.col[0], Score())
        for hit in self.col:
            score = hit.score(scorer)
            if(score.dat[by] >= best[1].dat[by]):
                best = (hit, score)
        return(best)

    def get_evalue(self):
        evalue = 999
        for hit in self.col:
            hevalue = min([hsp.mevalue for hsp in hit.col])
            if(hevalue is not None):
                if(hevalue < evalue):
                    evalue = hevalue
        return(evalue)

class Hit(_Collection):
    """\
    Contains a list of Hsp objects with the same database, query, and hit
    numbers
    """
    def __init__(self, hsp=None):
        self.hit = None
        self.hlen = 0
        super().__init__(hsp)

    def _first_item(self, hsp):
        super()._first_item(hsp)
        self.hit = hsp.hit
        self.hlen = hsp.hlen

    def score(self, scorer):
        return(scorer.score(self))

    @staticmethod
    def same_group(a, b):
        if(a.col[0].db == b.col[0].db and a.col[0].qseqid == b.col[0].qseqid):
            return True
        else:
            return False

class Hsp:
    """\
    A class that holds HSP data, a single row from the sql database
    """
    def __init__(self, row, col):
        # A dict with column name keys and row name values
        try:
            data = {col[i]:row[i] for i in range(len(row))}
            self.db = data['database']
            self.qseqid = data['qseqid']
            self.qdat = {x:data[x] for x in ('qgi', 'qgb', 'qgene', 'qlocus', 'qtaxon', 'qlen')}
            self.hit = data['hit']
            self.hlen = data['hlen']
            self.mevalue = data['mevalue']
            self.hsp = data['hsp']
            self.qto = data['qto']
            self.qfrom = data['qfrom']
            self.hto = data['hto']
            self.hfrom = data['hfrom']
            self.alen = data['alen']
            self.score = data['score']
            self.ident = data['ident']
            self.pos = data['pos']
            self.gaps = data['gaps']
        except IndexError as e:
            traceback.print_exc(file=sys.stderr)
            print(e, file=sys.stderr)
            print("HSP row values and col must be of equal length", file=sys.stderr)
            sys.exit(1)
        except KeyError as e:
            traceback.print_exc(file=sys.stderr)
            print(e, file=sys.stderr)
            print("Input columns must match column names as received from the SQL database",
                  file=sys.stderr)
            sys.exit(1)

    @staticmethod
    def same_group(a, b):
        if(a.db == b.db and a.qseqid == b.qseqid and a.hit == b.hit):
            return True
        else:
            return False

    @staticmethod
    def is_overlapping(a, b):
        """\
        Takes two Hsps and returns true if they overlap
        """
        # If a boundaries are within b
        if((a.hfrom >= b.hfrom and a.hfrom <= b.hto) or
           (a.hto >= b.hfrom and a.hto <= b.hto)):
            return True
        # Else if b boundaries are within a
        elif((b.hfrom >= a.hfrom and b.hfrom <= a.hto) or
             (b.hto >= a.hfrom and b.hto <= a.hto)):
            return True
        else:
            return False

    @staticmethod
    def contains_overlap(a, hsp_list):
        """\
        Takes a single Hsp and searches through a list of Hsps. If any overlap,
        returns true.
        """
        for b in hsp_list:
            if(Hsp.is_overlapping(a, b)):
                return True
        return False

class Score:
    def __init__(self, hsp=None):
        if(hsp):
            self.dat = {
                'alen': hsp.alen,
                'gaps': hsp.gaps,
                'ident': hsp.ident,
                'pos': hsp.pos,
                'score': hsp.score
            }
            if(hsp.score == 0):
                self.dat['nhsp'] = 0
            else:
                self.dat['nhsp'] = 1
        else:
            self.dat = {
                'nhsp': 0,
                'alen': 0,
                'gaps': 0,
                'ident': 0,
                'pos': 0,
                'score': 0
            }

    def add_prefixed_dat(self, dat, prefix=''):
        for col, val in self.dat.items():
            dat[prefix + col] = val

    @classmethod
    def fromdict(cls, d):
        """ Alterative constructor """
        score = Score()
        for key, val in d.items():
            score.dat[key] = val
        return(score)

    @staticmethod
    def max(scores, by='score'):
        best = Score()
        for score in scores:
            if(best.dat[by] < score.dat[by]):
                best = score
        return(best)
