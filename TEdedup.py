import re 
import os 
import sys
import argparse
from collections import defaultdict

__Author__  = "BiaoFeng Zhou"
__Update__  = '20250630'
__Version__ = '0.0.2'

#-------------------------------
# Description:
#
# set level
# level 6: intact         # intact LTR
# level 5: homology       # homology TE 
# level 4: LTR/unknon     # unknown LTR 
# level 3: structural     # Structural TE 
# level 2: Unspecified    # repeats annotated by EDTA 
# level 1: Tandem_Repeats # repeats annotated by TRF
#
# Prioritizes higher-identity matches when annotations are at the same level. For matches with equal identity, longer sequences are retained.
#-------------------------------
splitpattern = r'[\r\s]+'

class LineNumberIterator:
    def __init__(self, filepath):
        self.file = open(filepath, 'r')
        self.line_num = 0
    
    def __iter__(self):
        return self
    
    def __next__(self):
        line = self.file.readline()
        if not line:
            self.file.close()
            raise StopIteration
        self.line_num += 1
        return (self.line_num, line.strip())
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.file.close()

class FileWriter:
    def __init__(self, filepath):
        self.filepath = filepath
        if os.path.exists(self.filepath):
            os.remove(self.filepath)
            print(f"File {self.filepath} was presented, and has been deleted.")

    def append_content(self, content):
        """
        :param content: wirter
        """
        with open(self.filepath, 'a') as file: 
            # Check if content is a list. If it is a list, call file.writelines(content) to write each item in the list as a separate line in the file.
            if isinstance(content, list):
                file.writelines(content)
            else:
                file.write(content)

    def wirte_content(self, content):
        with open(self.filepath, 'w') as file:
            # Check if content is a list. If it is a list, call file.writelines(content) to write each item in the list as a separate line in the file.
            if isinstance(content, list): 
                file.writelines(content)
            else:
                file.write(content)


def parse_arguments():
    usage = "python split_overlap.py INBED OUTBED --intactgff GFF3"
    parser = argparse.ArgumentParser()
    parser.add_argument(
            'inbed', help='Input BED File. 1-Based.')
    parser.add_argument(
            'outbed', help='Output BED Name. 1-Based.')
    parser.add_argument("--intact-gff", dest = "intactgff", help = 'Intact GFF3 File.')
    parser.epilog = usage
    return parser.parse_args()

def parseIntact(gff):
    intactID = []
    with open(gff, 'r') as f:
        for line in f:
            ele = re.split(r'[\t\s]+', line.strip() )
            if re.search(r'_LTR_retrotransposon', line.strip(), re.IGNORECASE):
                TE_ID = re.search(r'Name=(.*?);', ele[0])
                intactID.append(TE_ID)
    return intactID

def buildBlocks(l):
    blist = []
    binfo = defaultdict(list)
    coordall = []
    for e in l:
        coordall.extend([ int(e[1]), int(e[2]) ])
    coordall = sorted(set(coordall))

    for num, coord in enumerate(coordall):
        if num == len(coordall) - 1: continue
        else: 
            s = coordall[num]
            nexts = coordall[num+1]
            if num == len(coordall) - 2: blist.append([s, nexts])  
            else: blist.append([s, nexts - 1]) 
    for b in blist: binfo["_".join([ str(i) for i in b ] )] = []
    return blist, binfo # blist : [[start, end], ....] ; binfo :{start"_"end:[], ... }


def dealOvBlocks(l, intact, writer):
    if len(l)  == 1:
        writer.append_content("\t".join(l[0]) + "\n")
        return 

    blist, bdict = buildBlocks(l)
    
    leveldict = defaultdict()
    lendict = defaultdict()
    identitydict = defaultdict()

    for num, e in enumerate(l):
        accession = "acc" + str(num)
        if e[5] == 'structural':
            if e[3] == 'Tandem_Repeats':
                leveldict[accession] = 1
            elif e[4] == 'Unspecified':
                leveldict[accession] = 2
            else:
                leveldict[accession] = 3
        elif e[5] == 'homology':
            if e[3].startswith('TE'):
                so = e[3].split('_')
                if len(so) > 2: soname = "_".join(so[:2])
                else: soname = e[3]
            else: soname = e[3]
            if soname in intact:
                leveldict[accession] = 6
            elif re.search(r'unknown', e[4]):
                leveldict[accession] = 4
            else: 
                leveldict[accession] = 5
        else: print('unknown type.')
        
        e.append(accession)  # index accession for each
        lendict[accession] = int(e[2]) - int(e[1]) + 1
        identitydict[accession] = float(e[6]) if e[6] != "NA" else 0.0
   
    bestblocks = [] 
    for eb in blist:
        ovlist = []
        for e in l:
            if checkoverlap([ eb[0], eb[1] ], [ e[1], e[2] ]):
                ovlist.append([ e[0] ] + eb + e[3:])

        cblock = []
        clevel = -1
        clength = 0
        cidentity = 0
        for eov in ovlist:
            if leveldict[ eov[-1] ] > clevel:
                cblock = eov
                clevel = leveldict[ eov[-1] ]
                clength = lendict[ eov[-1] ]
                cidentity = identitydict[ eov[-1] ]
            elif leveldict[ eov[-1] ] == clevel:   
                if identitydict[ eov[-1] ] > cidentity:
                    cblock = eov
                    clevel = leveldict[ eov[-1] ]
                    clength = lendict[ eov[-1] ]
                    cidentity = identitydict[ eov[-1] ] 
                elif identitydict[ eov[-1] ] == cidentity:
                    if lendict[ eov[-1] ] > clength:
                        cblock = eov
                        clevel = leveldict[ eov[-1] ]
                        clength = lendict[ eov[-1] ]
                        cidentity = identitydict[ eov[-1] ]
                    else: continue
                else: continue
            else:
                 continue
        
        bestblocks.append(cblock) 

    finalblocks = []
    for num, beb in enumerate(bestblocks):
        if num == len(bestblocks) - 1:
            continue
        else:
            if num == 0: cline = beb
            nextline = bestblocks[num + 1]
            if cline[-1] == nextline[-1]: # same accession, merge
                newline = cline.copy()
                newline[2] = nextline[2]
                cline = newline.copy()
            else: # write out
                writer.append_content("\t".join([ str(i) for i in cline ]) + "\n")
                cline = nextline.copy()
    writer.append_content("\t".join([ str(i) for i in cline ]) + "\n")

def checkoverlap(interval1, interval2):
    start1, end1 = [ int(i) for i in interval1 ]
    start2, end2 = [ int(i) for i in interval2 ]
    if (end2 - start1) * (end1 - start2) >=  0: # overlap
        return True
    else: return False

def find_unique_intervals(interval1, interval2):
    start1, end1 = interval1
    start2, end2 = interval2

    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)

    unique_intervals = []
    if start2 < overlap_start:
        unique_intervals.append([start2, overlap_start - 1])
    if end2 > overlap_end:
        unique_intervals.append([overlap_end + 1, end2])
    return unique_intervals

def main(args):
    iterator1 =  LineNumberIterator(args.inbed)
    writer = FileWriter(args.outbed)
    intactID = parseIntact(args.intactgff) if args.intactgff else [] 
    try:
        overlapblocks = [] 
        cstart, cend = 0, 0 
        cchr = ''
        chrlist = []
        num = 0
        overlap_c = 0
        while True:
            num += 1
            if num % 10000 == 0: print(f"TE number:{num} has completed ! ")
            line_num, line = next(iterator1) 
            ele = re.split(splitpattern, line)
            if ele[0] not in chrlist: 
                chrlist.append(ele[0]) 
                if num != 1: dealOvBlocks(overlapblocks, intactID, writer)
                cchr = ele[0]
                overlapblocks = [ele] 
                continue
            
            for e in overlapblocks.copy():
                if ( int(ele[2]) - int(e[1]) ) * (int(ele[1]) - int(e[2])) <= 0: # Overlap
                    overlap_c+=1
                    if overlap_c == 1: overlapblocks.append(ele)

            if overlap_c == 0:
                dealOvBlocks(overlapblocks, intactID, writer)
                overlapblocks = [ele] # reset 
                overlap_c = 0         # reset
                continue
            else:
                overlap_c = 0
                continue
    except StopIteration:
        dealOvBlocks(overlapblocks, intactID, writer)
        print('File Reader Finish. Exit.')

if __name__ == "__main__":
    args = parse_arguments()
    main(args)
