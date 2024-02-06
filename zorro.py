#!/usr/bin/python

from Bio import AlignIO
import glob
import numpy

threshold = 5 #float(raw_input('What HMM value do you which to exclude? '))

fastas = glob.glob('*.fasta')

for i in fastas:
        handle = open(i, 'rb')
        alignment = AlignIO.read(i, 'fasta')
        align_array = numpy.array([list(rec) for rec in alignment], numpy.character, order="F")
        OGname = i.rstrip('.fasta')
        zorro_handle = open(OGname + '.fasta.mask', 'rb')
        count = 0
        cut = align_array[:,:1]
        for j in zorro_handle:
                if float(j.rstrip('\n')) >= threshold:
                        cut = numpy.column_stack((cut, align_array[:,count]))
                count += 1
        cut = cut[:,1:]
        handle.close()
        zorro_handle.close()

        outfile = open(i + '.cut', 'w')
        for q in range(len(alignment)):
                seq = ''
                for p in cut[q,:]:
                        seq += p
                outfile.write('>' + alignment[q].id + '\n' + seq + '\n')
        outfile.close()
