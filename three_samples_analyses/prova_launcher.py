#!/usr/bin/python

import os
import sys

iterations = sys.argv[1]
lagging = sys.argv[2]

os.chdir('/Users/IEO5505/Desktop/MI_TO/three_samples_analyses/')
for i in range(int(iterations)):
    os.system(f'python prova.py {i} {lagging}') 