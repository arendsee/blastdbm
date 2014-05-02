#! /usr/bin/python3

import urllib.request
import urllib.parse
import re
import xml.etree.ElementTree as et
import sys
import time
from orphanlib.lineage import Lineage

def BASE_URL():
    return('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/{}.fcgi')


# ==================
# EXPORTED FUNCTIONS
# ==================

def sciname2taxid(name):
    val = {'db': 'taxonomy', 'term': name}
    xml = _query('esearch', val)
    taxid = _simple_search('Id', xml)
    return(taxid)

def taxid2lineage(taxids):
    if(not isinstance(taxids, (list, tuple))):
        taxids = (taxids, )
    val = {'db' : 'taxonomy', 'id' : ','.join(map(str, taxids)),
           'rettype' : 'xml', 'retmode' : 'text'}
    xml = _query('efetch', val)
    root = et.fromstring(xml)
    out = []
    for taxon in root.findall('./Taxon'):
        taxid = taxon.find('TaxId').text
        name = taxon.find('ScientificName').text
        lineage = []
        for ancestor in taxon.findall('./LineageEx/'):
            ataxid = ancestor.find('TaxId').text
            aname = ancestor.find('ScientificName').text
            lineage.append((ataxid, aname))
        lin_obj = Lineage(taxid, name, lineage)
        out.append(lin_obj)
    return(out)

def taxid2sciname(taxid):
    xml = _query('efetch', {'db': 'taxonomy', 'id': taxid})
    sciname = _simple_search('ScientificName', xml)
    return(sciname)


# =================
# UTILITY FUNCTIONS
# =================

def _query(cmd, val):
    url = BASE_URL().format(cmd)
    arg = urllib.parse.urlencode(val)
    arg = arg.encode('ascii')
    req = urllib.request.Request(url, arg)
    try:
        # In order to comply with ENTREZ policy
        time.sleep(1)
        f = urllib.request.urlopen(req)
    except:
        print("Failed to retrieve data from entrez", file=sys.stderr)
        sys.exit(1)
    body = f.read().decode('ascii')
    return(body)

def _simple_search(tag, xml):
    """ Returns the string bound by the first appearance of a tag, this may
        be an XML
    """

    xml = re.sub("\n", " ", xml)
    m = re.search(r"<{0}>\s*(.*?)\s*<\/{0}>".format(tag), xml)
    try:
        text = m.group(1)
    except:
        text = None
    return(text)
