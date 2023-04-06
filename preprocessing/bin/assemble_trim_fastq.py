#!/usr/bin/python

import sys
import dnaio

# Args
R1 = sys.argv[1]
R2 = sys.argv[1]

#Open
reader = dnaio.open(R1, file2=R2)
writer = dnaio.open('assembled.fastq.gz', mode='w')

#Iterate
for r1, r2 in reader:
    #Get cbc, umi
    cbc, umi = r1.sequence[:16], r1.sequence[16:]
    #Trim R2
    r2.sequence, r2.qualities = r2.sequence[24:], r2.qualities[24:]
    r2.name = '_'.join( r2.name.replace('/', '.').split(' ') + [cbc] + [umi] )
    writer.write(r2)

# Close
reader.close()
writer.close()

