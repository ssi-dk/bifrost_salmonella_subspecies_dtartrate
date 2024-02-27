#!/usr/bin/env python3

import argparse
import re
import io
import sys

def get_args():
    parser=argparse.ArgumentParser(description="Extracts bases at a certain position from a sam file")
    parser.add_argument('file',type=argparse.FileType())
    parser.add_argument('--pos',type=int, help="Position to extract")
    g = parser.add_mutually_exclusive_group()
    g.add_argument("--length",type=int, help="Length of sequence to extract")
    g.add_argument("--end",type=int, help="Last position to extract")
    return parser.parse_args()

def get_pos_from_sam(pos: int, end_pos: int, fh: io.TextIOBase):
    cigar_re = re.compile("(\d+)S.+")
    bases = list()
    for line in fh:
        if line.startswith("@"):
            continue
        fields = line.split()
        try:
            seq = fields[9]
        except IndexError:
            print(line,file=sys.stderr)
            continue
        cigar = fields[5]
        alignment_start = int(fields[3])
        mo = cigar_re.match(cigar)
        if alignment_start == 1 and mo is not None:
            overhang = int(mo.group(1))
            bases.append(seq[(overhang+pos-1):(overhang+end_pos-1)])
        else:
            if pos>=alignment_start:
                bases.append(seq[(pos-alignment_start):(end_pos-alignment_start)])
    return bases

def main(args):
    pos = args.pos
    if args.end:
        end_pos = args.end+1
    elif args.length:
        end_pos = pos+args.length
    else:
        end_pos = pos+1
    print(''.join(get_pos_from_sam(pos, end_pos, args.file)))

if __name__ == "__main__":
    args = get_args()
    main(args)