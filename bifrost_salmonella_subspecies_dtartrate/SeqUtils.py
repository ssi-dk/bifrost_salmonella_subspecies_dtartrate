#!/usr/bin/env python3

# Utils.py Utility library for various useful classes and functions.
#    Copyright (C) 2013  Kristoffer Kiil

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.


import gzip
import os.path
import re
import collections
import io
import sys

tr = bytes.maketrans(b"atgcATGC", b"tacgTACG")
ascii = 'ascii'
codontable = dict()
codontable[11] = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
    }


class FqStream:
    def __init__(self, fname=None, fh=None, outQualBase=33, inQualBase=None):
        self.fh = fh
        if fname:
            try:
                self.fh = gzip.open(fname)
                self.fh.peek(1)
            except OSError:
                self.fh = open(fname, 'rb')
        if not inQualBase is None:
            self.inQualBase = inQualBase
        else:
            self.inQualBase = None
        self.outQualBase = outQualBase
        self.buffer = dict()

    def __next__(self):
        return Fq(self)

    def __iter__(self):
        return self

    def readline(self):
        line = self.fh.readline()
        if not line:
            raise StopIteration
        else:
            return line

    def readlines(self, n):
        res = list()
        while(len(res) < n):
            r = self.fh.readline()
            if not r:
                raise StopIteration
            elif(len(r) > 1):
                res.append(r.strip(b'\n'))
        return res


class QualError(Exception):
    def __init__(self, errstr="QualityError"):
        Exception.__init__(self)


class Fq:
    def __init__(self, parent):
        self.id = b''
        self.seq = b''
        self.comment = b''
        self.qual = []
        try:
            self.id = parent.readline().lstrip(b'@').strip(b'\n')
            s = []
            i = 0
            while True:
                line = parent.readline().strip(b'\n')
                if line.startswith(b'+'):
                    self.comment = line
                    break
                i += 1
                s.append(line)
            self.seq = b"".join(s)
            for j in range(i):
                self.qual.extend(parent.readline().strip(b'\n'))
        except AttributeError:
            pass
        self.parent = parent
        self.setQual()

    def __len__(self):
        return len(self.seq)

    def __str__(self):
        try:
            buffer = self.parent.buffer['write']
        except KeyError:
            self.parent.buffer['write'] = list()
            buffer = self.parent.buffer['write']
        buffer.append(self)
        
        iQB = self.parent.inQualBase
        oQB = self.parent.outQualBase
        if not (iQB is None):
            s = b"\n".join(
                [b"\n".join(
                    (b'@' + fq.id, fq.bseq(), fq.comment,  bytes(
                        [c-iQB+oQB for c in fq.qual]))) for fq in buffer]
            ).decode('ascii')
            buffer.clear()
            return s
        else:
            raise QualError()

    def sseq(self, seq=None):
        if seq:
            self.seq = seq.encode(ascii)
        else:
            return self.seq.decode()

    def bseq(self, seq=None):
        if seq:
            self.seq = seq
        else:
            return self.seq

    def squal(self, qual=None):
        if qual:
            self.qual = list(bytes(qual))
        else:
            iQB = self.parent.inQualBase
            oQB = self.parent.outQualBase
            return bytes([c-iQB+oQB for c in self.qual]).decode(ascii)

    def bqual(self, qual=None):
        if qual:
            self.qual = list(qual)
        else:
            iQB = self.parent.inQualBase
            oQB = self.parent.outQualBase
            return bytes([c-iQB+oQB for c in self.qual])

    def iqual(self, qual=None):
        if qual:
            self.qual = qual
        else:
            return self.qual

    def sid(self, id=None):
        if id:
            self.id = id.encode(ascii)
        else:
            return self.id.decode(ascii)

    def bid(self, id=None):
        if id:
            self.id = id
        else:
            return self.id

    def setQual(self):
        if self.parent.inQualBase is None:
            self.parent.inQualBase = self.detectQuality()

    def write(self, fh=sys.stdout.buffer):
        try:
            buffer = self.parent.buffer['write']
        except KeyError:
            self.parent.buffer['write'] = list()
            buffer = self.parent.buffer['write']
        buffer.append(self)
        iQB = self.parent.inQualBase
        oQB = self.parent.outQualBase
        if iQB:
            while len(buffer) > 0:
                fq = buffer.pop(0)
                fh.write(
                    b"@"+b"\n".join(
                        (fq.id, fq.seq, fq.comment,
                         bytes([c-iQB+oQB for c in fq.qual])))+b"\n")
        else:
            raise QualError()

    def detectQuality(self):
        for q in self.iqual():
            if q < 59:
                # q=59 is minimum possible quality with base 64 (weird but true)
                print("Detected quality scores base 33 ({})".format(chr(q)),file=sys.stderr)
                return 33
            elif q>74:
                # q=74 is maximum possible quality with base 33
                print("Detected quality scores base 64 ({})".format(chr(q)),file=sys.stderr)
                return 64
        # In case the base quality could not be determined return None
        return None

    def trim(self, ltrm=50, rtrm=None):
        if rtrm is None:
            rtrm = ltrm
        if rtrm == 0:
            self.bseq(self.bseq()[ltrm:])
            self.iqual(self.iqual()[ltrm:])
        else:
            self.bseq(self.bseq()[ltrm:-rtrm])
            self.iqual(self.iqual()[ltrm:-rtrm])
        return self

    def ltrim(self, length=0):
        self.bseq(self.bseq()[length:])
        self.iqual(self.iqual()[length:])
        return self

    def qtrim(self, qual=25, window=10):
        try:
            buffer = self.parent.buffer['qtrim']
        except KeyError:
            self.parent.buffer['qtrim'] = list()
            buffer = self.parent.buffer['qtrim']
        buffer.append(self)
        iQB = self.parent.inQualBase
        if iQB:
            while len(buffer) > 0:
                fq = buffer.pop(0)
                l = len(fq)
                ltrim = 0
                t = (qual+iQB)*window
                m = int(l*0.25)
                s = sum(fq.iqual()[m-window:m])
                for i in range(m, 0+(window-1), -1):
                    if s < t:
                        ltrim = i-int(window/2)
                        break
                    s += fq.iqual()[i-window]-fq.iqual()[i]
                rtrim = l
                s = sum(fq.iqual()[m:m+window])
                for i in range(m, l-(window)):
                    if s < t:
                        rtrim = i+int(window/2)
                        break
                    s += fq.iqual()[i+window]-fq.iqual()[i]
                fq.comment = b' '.join(
                    (self.comment, b"trimmed", str(ltrim).encode('ascii'),
                     str(l-rtrim).encode(ascii), b''))
            fq.trim(ltrim, l-rtrim)
        return self

    def fsa(self):
        return Fsa(b'\n'.join((self.id, self.bseq())))

    def maskByQual(self, threshold):
        iQB = self.parent.inQualBase
        mask = [i for i in range(len(self.qual)) if self.qual[i]-iQB<threshold]
        mutseq = list(self.seq)
        for i in mask:
            mutseq[i] = ord(b'N')
        self.seq = bytes(mutseq)
        return self


def readfsa(fh):
    """Reads a file and returns an fsa object"""
    raw = list()
    seqs = list()
    for line in fh:
        if line.startswith(">") and len(raw) > 0:
            seqs.append(Fsa("".join(raw)))
            raw.clear()
        raw.append(line)
    if len(raw) > 0:
        seqs.append(Fsa("".join(raw)))
    return seqs

class FsaStream:
    def __init__(self, fh):
        self.fh = fh
        self.buf = list()

    def __iter__(self):
        return self

    def __next__(self):
        for line in self.fh:
            if type(line) == type(""):
                line = line.encode()
            if line.startswith(b">") and len(self.buf) > 0:
                fsa = Fsa(b"".join(self.buf))
                self.buf.clear()
                self.buf.append(line)
                return fsa
            self.buf.append(line)
        if len(self.buf) > 0:
            fsa = Fsa(b"".join(self.buf))
            self.buf.clear()
            return fsa
        else:
            raise StopIteration


class Fsa:
    reCodes = dict()
    ambiCodes = {ord('A'): 'A', ord('C'): 'C', ord('G'): 'G', ord('T'): 'T',
                 ord('M'): '[ACM]', ord('R'): '[AGR]', ord('W'): '[ATW]',
                 ord('S'): '[CGS]', ord('Y'): '[CTY]', ord('K'): '[GTK]',
                 ord('B'): '[CTGSYKB]', ord('D'): '[AGTRWKD]',
                 ord('H'): '[ATCMWYH]', ord('V'): '[ACGMRSV]',
                 ord('N'): '[ATGCMRWSYKBDHVN]'}

    def __init__(self, raw):
        if isinstance(raw, str):
            raw = raw.encode()
        try:
            (self.id, self.seq) = raw.replace(b'\r', b'').split(b'\n', 1)
        except ValueError:
            print(raw, sys.stderr)
        self.id = self.id.lstrip(b'>')
        self.seq = self.seq.replace(b'\n', b'')
        self.trans = tr
        self.linelength = 60
        self.PP = True

    def __len__(self):
        return len(self.seq)

    def __str__(self):
        if not self.PP:
            return b'\n'.join((b'>'+self.id, self.seq)).decode()
        linelength = self.linelength
        pos = 0
        ls = list()
        while pos+linelength < len(self.seq):
            ls.append(self.seq[pos:pos+linelength])
            pos += linelength
        ls.append(self.seq[pos:])
        return b'\n'.join((b'>'+self.id, b'\n'.join(ls))).decode()

    def sid(self, id=None):
        if id:
            self.id = id.encode(ascii)
        else:
            return self.id.decode(ascii)

    def bid(self, id=None):
        if id:
            self.id = id
        else:
            return self.id

    def sseq(self, seq=None):
        if seq:
            self.seq = seq.encode()
        else:
            return self.seq.decode()

    def bseq(self, seq=None):
        if seq:
            self.seq = seq
        else:
            return self.seq

    def tab(self):
        return "\t".join((self.sid(), self.sseq()))

    def extract(self, start, end):
        '''Extract part of the sequence in base zero'''
        return Fsa(b"\n".join((self.id,self.bseq()[start:end])))
    def reverse(self):
        return Fsa(b"\n".join((self.id,self.bseq()[::-1].translate(self.trans))))
    def re(self,flags= None):
        seq=self.bseq()
        if seq:
            s=list()
            for c in seq:
                if c in Fsa.reCodes:
                    s.append(Fsa.reCodes[c])
                else:
                    ch="[{}N]".format(chr(c))
                    Fsa.reCodes[c]=ch
                    s.append(ch)
            if flags:
                regexp=re.compile("".join(s),flags)
            else:
                regexp=re.compile("".join(s))
            return regexp

    def ambi2re(self,flags=None):
        regexp=[]
        seq=self.bseq()
        if seq:
            s = "".join([Fsa.ambiCodes[c] for c in seq])
            if flags:
                regexp=re.compile(s,flags)
            else:
                regexp=re.compile(s)
            return regexp
    def translate(dnaseq, table=11):
        proteinseq = list()
        for n in range(0, len(dnaseq), 3):
            proteinseq.append(codontable[table][dnaseq[n:n+3]])
        return "".join(proteinseq)

class Fsa2:
    def __init__(self,raw):
        lines= raw.split("\n")
        self.id=lines[0].lstrip(">")
        self.seq="".join(lines[1:])
        self.trans=tr
    def __len__(self):
        return len(self.seq)
    def __str__(self):
        linelength=60
        pos=0
        ls=list()
        while pos+linelength<len(self.seq):
            ls.append(self.seq[pos:pos+linelength])
            pos+=linelength
        ls.append(self.seq[pos:])
        return ">{}\n{}\n".format(self.id,"\n".join(ls))
    def extract(self,start,end):
        '''Extract part of the sequence in base zero'''
        return Fsa(self.id + '\n' + self.seq[start:end])
    def reverse(self):
        return Fsa(self.id + '\n' + self.seq[::-1].translate(self.trans))
    def re(self,flags=None):
        regexp=list()
        for c in self.seq:
            regexp.append("[{}N]".format(c))
        regexp=re.compile("".join(regexp),flags)
        return regexp
    def ambi2re(self,flags=None):
        ambicodes={'A':'[A]', 'C':'[C]', 'G':'[G]', 'T':'[T]', 'M':'[ACM]', 'R':'[AGR]', 'W':'[ATW]', 'S':'[CGS]', 'Y':'[CTY]', 'K':'[GTK]', 'B':'[CTGSYKB]', 'D':'[AGTRWKD]', 'H':'[ATCMWYH]', 'V':'[ACGMRSV]', 'N':'[ATGCMRWSYKBDHVN]', }
        regexp=[]
        for c in self.seq:
            regexp.append(ambicodes[c])
        regexp=re.compile("".join(regexp),flags)
        return regexp

class alignment(list):
    def __init__(self,positions=None):
        list.__init__(self)
        self.forspalte=list()
        if positions:
            self.positions=positions
        else:
            self.positions=list()
    def addpositions(self,positions):
        self.positions=positions
    def readphylip(self,string):
        lines=string.split('\n')
        (self.norgs,self.alnlength) = lines.pop(0).split()
        self.norgs=int(self.norgs)
        self.alnlength=int(self.alnlength)
        org=0
        for i in range(self.alnlength):
            self.append(list("N" * self.norgs))
        for line in lines:
            line=line.strip()
            self.forspalte.append(line[:10].strip())
            pos=0
            for char in line[10:]:
                try:
                    self[pos][org]=char
                except IndexError:
                    print("norgs={}\norg={}\nalnlength={}\ni={}".format(self.norgs,org,self.alnlength,i))
                pos+=1
            org+=1
    def readSNPmat(self,string,norgs=None):
        lines=string.split('\n')
        if(norgs):
            self.forspalte=lines.pop(0).strip().split("\t",maxsplit=norgs+1)[1:-1] 
        else:
            self.forspalte=lines.pop(0).strip().split("\t")[1:]
        for line in lines:
            if not line:
                continue
            line=line.strip()
            if norgs:
                self.append(line.split("\t",maxsplit=norgs+1)[1:-1])
            else:
                self.append(line.split("\t")[1:])
            self.positions.append(line.split("\t",maxsplit=1)[0])
        self.alnlength=len(self.positions)
        self.norgs=len(self.forspalte)
    def phylip(self):
        outstr="{} {}\n".format(self.norgs,self.alnlength)
        for org in range(self.norgs):
            outstr+="{: <9s} ".format(self.forspalte[org][:9])
            outstr+="{}\n".format("".join([self[pos][org] for pos in range(self.alnlength)]))
        return outstr.rstrip('\n')
    def nexus(self):
        outstr="#NEXUS\nBegin data;\nDimensions ntax={} nchar={};\nFormat datatype=dna missing=? gap=-;\nMatrix\n".format(self.norgs,self.alnlength)
        for org in range(self.norgs):
            outstr+="{} ".format(self.forspalte[org])
            outstr+="{}\n".format("".join([self[pos][org] for pos in range(self.alnlength)]))
        outstr += ";\nEnd;"
        return outstr
    def fasta(self):
        outlist = list()
        for org in range(self.norgs):
            outlist.append(str(Fsa(">{}\n{}".format(self.forspalte[org],"".join([self[pos][org] for pos in range(self.alnlength)])))))
        return "\n".join(outlist)
    def __str__(self,alnfrmt="phylip"):
        outputfunc=getattr(self,alnfrmt)
        return outputfunc()

def rebranch(substr):
    result=list()
    while substr:
        if substr[0]=="(":
            pos=select_parenthesis(substr)
            res=rebranch(substr[1:pos])
            result.append(res)
            
            substr=substr[pos+1:]
        elif substr[0]==",":
            substr=substr[1:]
        else:
            try:
                pos=substr.index(',')
                if len(substr[:pos])>0:
                    result.append(substr[:pos])
                substr=substr[pos+1:]
            except ValueError:
                if len(substr)>0:
                    result.append(substr)
                return result
    return result

def select_parenthesis(string):
    i=0
    for j in range(len(string)):
        if string[j]==')':
            i-=1
        elif string[j]=='(':
            i+=1
        if i==0:
            return j

def kdiff(str1,str2,k=10):
    """ Computes a kmer difference between two strings """
    cnt=0.0
    for i in range(len(str1)-k):
        cnt += 1 if str2.find(str1[i:i+k])<0 else 0
    return (cnt/(len(str1)-k))

def gzopen(name,mode='r'):
    """ Opens both gzip and plain files transparantly """
    try:
        fh=gzip.open(name)
        line=next(fh)
        fh.seek(0)
        return Fh(fh,gz=True)
    except  OSError:
        fh.close()
        fh=open(name)
    return fh

class Fh:
    """ File handle that can handle gz streams transparantly """
    def __init__(self,fh,gz=False):
        self.fh=fh
        self.gz=gz
        if self.gz:
            self.nextfunc=self.gznext
        else:
            self.nextfunc=self.nnext
    def __iter__(self):
        return self
    def __next__(self):
            return self.nextfunc()
    def gznext(self):
        return next(self.fh).decode()
    def nnext(self):
        return next(self.fh)
    def close(self):
        self.fh.close()
    def seek(self,n):
        self.fh.seek(n)


class Filename:

    def __init__(self, path):
        self.path = path
        (self.dir, self.basename) = os.path.split(path)
        self.absdir = os.path.abspath(self.dir)
        pattern = re.compile('(.*)_S(\d+)_L(\d+)_R(\d+)_(\d+)(.*)')
        mo = pattern.match(self.basename)
        if(mo):
            self.illumina = True
            self.ID = mo.group(1)
            self.S = mo.group(2)
            self.L = mo.group(3)
            self.R = mo.group(4)
            self.N = mo.group(5)
            self.ext = mo.group(6)
        else:
            self.illumina = False
            self.ID = self.basename.rstrip(".gz").rstrip(".fastq")
            self.ext = self.basename.lstrip(self.ID)

    def fullname(self):
        if self.illumina:
            return "{}_S{}_L{}_R{}_{}{}".format(
                self.ID, self.S, self.L, self.R, self.N, self.ext)
        else:
            return "{}{}".format(self.ID, self.ext)

    def __str__(self):
        return self.path


def seqs2restring(seqs):
    maxlength = max(map(len, seqs))
    regexp = []
    for i in range(maxlength):
        pos = set()
        for s in seqs:
            pos.add(s[i])
        regexp.append("[{}]".format("".join([chr(c) for c in pos])))
    return "".join(regexp)


def seqs2re(seqs, flags=None):
    return re.compile(seqs2restring(seqs), flags)


def pcr2re(start, end, minlen, maxlen, flags=None):
    startstr = seqs2restring(start)
    endstr = seqs2restring(end)
    regexp = startstr + ".{" + str(
        minlen-len(startstr)-len(endstr)
    ) + "," + str(
        maxlen-len(startstr)-len(endstr)
    ) + "}" + endstr
    print(regexp, file=sys.stderr)
    return re.compile(regexp, flags)


def get_suffix_array(str):
    """ Don't use for anything serious, as it's really slow and hogs memory """
    return sorted(range(len(str)), key=lambda i: str[i:])


def sort_bucket(str, bucket, order):
    d = collections.defaultdict(list)
    for i in bucket:
        key = str[i:i+order]
        d[key].append(i)
    result = []
    for k, v in sorted(d.iteritems()):
        if len(v) > 1:
            result += sort_bucket(str, v, order*2)
        else:
            result.append(v[0])
    return result


def suffix_array_ManberMyers(str):
    return sort_bucket(str, (i for i in range(len(str))), 1)


def translate(dnaseq, table=11):
    proteinseq = list()
    for n in range(0, len(dnaseq), 3):
        proteinseq.append(codontable[table][dnaseq[n:n+3]])
    return "".join(proteinseq)
