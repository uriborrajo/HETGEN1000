#!/usr/bin/python

from Bio import AlignIO
import glob
import numpy as np

threshold = 5

fastas = glob.glob('*.fasta')

for i in fastas:
    handle = open(i, 'rb')
    alignment = AlignIO.read(i, 'fasta')
    align_array = np.array([list(rec) for rec in alignment], 'S1', order="F")
    OGname = i.rstrip('.fasta')
    zorro_handle = open(OGname + '.fasta.mask', 'rb')
    count = 0
    cut = align_array[:,:1]
    for j_bytes in zorro_handle:
        j_str = j_bytes.decode('utf-8') 
        if float(j_str.rstrip('\n')) >= threshold:
            cut = np.column_stack((cut, align_array[:, count]))
        count += 1
    cut = cut[:,1:]
    handle.close()
    zorro_handle.close()

    outfile = open(i + '.cut', 'w')
    for q in range(len(alignment)):
        seq = ''
        for p in cut[q,:]:
            seq += p.decode('utf-8')  
        outfile.write('>' + alignment[q].id + '\n' + seq + '\n')
    outfile.close()

