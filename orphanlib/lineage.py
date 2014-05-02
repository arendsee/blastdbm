#! /usr/bin/python3

import sys

class Lineage:
    def __init__(self, taxid, sciname, lineage):
        try:
            self.taxid = int(taxid)
        except ValueError:
            print("Taxid must be an integer", file=sys.stderr)
            sys.exit(1)
        self.sciname = sciname
        self.lineage = lineage

    def print(self):
        print("taxid: {} sciname: {} lineage: {}".format(
                self.taxid, self.sciname, self.lineage))

class MRCA:
    def __init__(self, lin1, lin2):
        self.taxid    = None
        self.sciname  = None
        self.level    = None

        a, b = lin1.lineage, lin2.lineage
        if(lin1.taxid == 1 or lin2.taxid == 1):
            self.taxid   = 1
            self.sciname = 'root'
            self.level   = 0
        elif(lin1.taxid == lin2.taxid):
            self.taxid   = lin1.taxid
            self.sciname = lin1.sciname
            # +1 for species
            self.level   = len(a) + 1
        else:
            maxlen = min(len(a), len(b))
            for i in range(maxlen):
                if(a[i] != b[i]):
                    self.taxid, self.sciname = a[i - 1]
                    self.level = i
                    break
                elif(i == (maxlen - 1)):
                    self.taxid, self.sciname = a[i]
                    self.level = i + 1

    def print(self):
        print("taxid: {} sciname: {} level: {}".format(
                self.taxid, self.sciname, self.level))

