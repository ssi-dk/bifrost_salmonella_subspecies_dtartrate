#!/usr/bin/env python3

import argparse
import sys
from bifrost_salmonella_subspecies_dtartrate import SeqUtils
import os, glob
import re
from pathlib import Path

resource_dir = "/bifrost/components/bifrost_salmonella_subspecies_dtartrate"
def parse_args():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--st', type=int)
    group.add_argument('--file', type=argparse.FileType('r'))
    group.add_argument('--sequence', type=str)
    parser.add_argument('--references', type=argparse.FileType('r'), metavar='file',
        default=resource_dir + "/resources/salmonella_subspecies_STs.txt")
    parser.add_argument('--mlstdir',type=str,
                        default=resource_dir + "/resources/salmonella",
                        help="Directory containing MLST definitions.")
    parser.add_argument('--fetch-mlst', action='store_true')
    parser.add_argument('--species', type=str, default="Salmonella enterica")
    parser.add_argument('--print-all', action='store_true')
    parser.add_argument('--verbose', action='store_true')
    return parser.parse_args()

def updateMLST(species, mlstdir):
    curdir = os.getcwd()
    os.chdir(mlstdir)
    os.system("getmlst.py --species \"" + species + "\"")
    os.chdir(curdir)
    
def parseMLSTseqs(mlstFile):
    with open(mlstFile,'r') as fh:
        fsaStream = SeqUtils.FsaStream(fh)
        db = dict()
        for fsa in fsaStream:
            (locus,number) = re.split("[-_]",fsa.sid())
            try:
                db[locus][int(number)] = fsa.sseq()
            except KeyError:
                db[locus] = {int(number): fsa.sseq()}
    return db

def parseMLSTalleles(mlstTable):
    with open(mlstTable, 'r') as fh:
        header = next(fh)
        header = header.strip().split()
        STidx = header.index('ST')
        alleles = dict()
        for line in fh:
            fields = line.strip().split()
            for i in range(len(fields)):
                if header[i] in ('CC','ST'):
                    continue
                try:
                    alleles[int(fields[STidx])][header[i]] = int(fields[i])
                except KeyError:
                    alleles[int(fields[STidx])] = {header[i]: int(fields[i])}
        return (alleles, header)

def parseRef(reffile):
    header = next(reffile).split()
    data = dict()
    subsp = dict()
    for line in reffile:
        fields = line.split()
        data[fields[0]] = dict()
        for i in range(1,len(header)):
            if header[i]=='Subsp':
                subsp[fields[0]] = fields[i]
            else:
                data[fields[0]][header[i]] = int(fields[i])
    return (data, subsp)

def parseAlleles(fh):
    alleles = dict()
    for seq in SeqUtils.FsaStream(fh):
        locus = seq.sid().split('-')[0]
        alleles[locus] = seq.sseq()
    loci = list(alleles.keys())
    loci.sort()
    seq = list()
    for locus in loci:
        seq.append(alleles[locus])
    return "".join(seq)

def alleles2seq(mlstdb, alleles):
    seq = list()
    loci = list(alleles.keys())
    loci.sort()
    for locus in loci:
        if locus != 'ST':
            seq.append(mlstdb[locus][alleles[locus]])
    return "".join(seq)

def mlstscore(ref, target, sw):
    score = 0
    refseq = alleles2seq(mlstdb, ref)
    targetseq = alleles2seq(mlstdb, target)
    if opts.verbose:
        print(refseq, file=sys.stderr)
    score = SeqUtils.kdiff(refseq, targetseq)
    return score

def seqscore(ref, query, mlstdb, sw):
    score = 0
    refseq = alleles2seq(mlstdb, ref)
    score = SeqUtils.kdiff(refseq, query)
    return score

def read_mlst_database(mlstDir):    
    mlstdb = dict()
    assert mlstDir.exists()
    for mlstFile in mlstDir.glob("*.tfa"):
        mlstdb.update(parseMLSTseqs(mlstFile))
    for mlstTable in mlstDir.glob("*.txt"):
        (alleles, header) = parseMLSTalleles(mlstTable)
    return (mlstdb, alleles, header)

def score_subspecies(query, refalleles, subsp, mlstdb, verbose=False):
    result = list()
    refkeys = list(refalleles.keys())
    refkeys.sort()
    for ref in refkeys:
        if verbose:
            print(">{}".format(ref), file=sys.stderr)
        score = seqscore(refalleles[ref], query, mlstdb, False)
        result.append((score, subsp[ref]))
    result.sort(key=lambda x:x[0])
    return result

def query_from_ST(ST, mlstdb, alleles):
    query = alleles2seq(mlstdb,alleles[ST])
    return query

def subspecies_from_query(query, refFile, mlstdb):
    (reference_alleles, subspecies) = parseRef(refFile)
    result = score_subspecies(query, reference_alleles, subspecies, mlstdb)
    return result

def subspecies_from_st(ST, mlstDir, reffile):
    mlstdb, alleles, header = read_mlst_database(mlstDir)
    return subspecies_from_query(query_from_ST(ST, mlstdb, alleles), open(reffile,'r'), mlstdb)

if __name__ == "__main__":
    opts = parse_args()
    refFile = Path(opts.references)
    mlstDir = Path(opts.mlstdir)
    species = opts.species
    ## Update MLST database
    if opts.fetch_mlst:
        updateMLST(species, mlstDir)
    ## Read MLST database
    mlstdb, alleles, header = read_mlst_database(mlstDir)
    ## Get query Seq
    if opts.st:
        query = alleles2seq(mlstdb,alleles[opts.st])
    elif opts.sequence:
        sequence = [int(x) for x in opts.sequence.split(',')]
        input_alleles = dict()
        for i in range(len(sequence)):
            input_alleles[[x for x in header if x not in ['ST','CC']][i]]=sequence[i]
        query = alleles2seq(mlstdb,input_alleles)
    else:
        query = parseAlleles(opts.file)
    result = subspecies_from_query(query, refFile)
    if opts.print_all:
        for (score,subspecies) in result:
            print("{}\t{}".format(score,subspecies))
    else:
        print(result[0][1])


