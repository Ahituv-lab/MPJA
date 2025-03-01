from collections import defaultdict
import sys
file1 = open(sys.argv[1], 'r')
lines = file1.readlines()
variants = defaultdict(lambda:defaultdict(int))
totalCount=0
import numpy as np
#import pandas as pd
for line in lines:
	field = line.split()
	fields = line.split("\t")
	count = int(field[0])
	md = eval(fields[-1])
	if md:
		md1 = np.array(md) 
		md1 = np.delete(md1,1,axis=1)
	if not md:
		md1=[]
	str1 = " "
	str1 = ','.join(map(str, md1))
	if variants[str1]:
		variants[str1] += count
	else:
		variants[str1] = count

for key in variants:
	print("%s\t%d"%(key,variants[key]))


