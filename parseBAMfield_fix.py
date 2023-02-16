#!/usr/bin/env python

import sys
def isSoftClipped(cigar):
  #Op BAM Description
  #M 0 alignment match (can be a sequence match or mismatch)
  #I 1 insertion to the reference
  #D 2 deletion from the reference
  #N 3 skipped region from the reference
  #S 4 soft clipping (clipped sequences present in SEQ)
  #H 5 hard clipping (clipped sequences NOT present in SEQ)
  #P 6 padding (silent deletion from padded reference)
  #= 7 sequence match
  #X 8 sequence mismatch
  for (op,count) in cigar:
    if op == 4: return True
  return False

def aln_length(cigarlist):
  tlength = 0
  for operation,length in cigarlist:
    if operation == 0 or operation == 2 or operation == 3 or operation >= 6: tlength += length
  return tlength

def parseMD(MD):
  MDfields = []
  value = ""
  chars = ""
  for elem in MD:
    if elem.isdigit() and chars == "":
      value+=elem
    elif not elem.isdigit() and value == "":
      chars+=elem
    elif elem.isdigit() and chars != "":
      MDfields.append(chars)
      chars = ""
      value = elem
    elif not elem.isdigit() and value != "":
      MDfields.append(int(value))
      value = ""
      chars = elem
  if value != "": MDfields.append(int(value))
  if chars != "": MDfields.append(chars)
  return MDfields

def parseMDwithCigar(MD,cigar,seq,alnStart):
  MDfields = parseMD(MD)
  res = []
  posRef = alnStart
  posRead = 0
  ind = 0 # Pointer in MDfields
  #Op BAM Description
  #M 0 alignment match (can be a sequence match or mismatch)
  #I 1 insertion to the reference
  #D 2 deletion from the reference
  #N 3 skipped region from the reference
  #S 4 soft clipping (clipped sequences present in SEQ)
  #H 5 hard clipping (clipped sequences NOT present in SEQ)
  #P 6 padding (silent deletion from padded reference)
  #= 7 sequence match
  #X 8 sequence mismatch
  for (op,count) in cigar:
    #print op,ind
    if (op in [0,7,8]) and ind % 2 == 1:
      while (count > 0):
        for base in MDfields[ind]:
          res.append((posRef,posRead,base+">"+seq[posRead]))
          posRef += 1
          posRead += 1
          count-= 1
        ind += 1 # Next posRef should be a number
        if (count > 0):
          if MDfields[ind] > count:
            posRef += count
            posRead += count
            MDfields[ind] -= count
            count = 0
          else:
            count -= MDfields[ind]
            posRef += MDfields[ind]
            posRead += MDfields[ind]
            ind += 1
      if (ind < len(MDfields)) and (MDfields[ind] == 0):
        ind += 1
    elif (op in [0,7,8]) and ind % 2 == 0:
      if MDfields[ind] > count:
        posRef += count
        posRead += count
        MDfields[ind] -= count
      else:
        count -= MDfields[ind]
        posRef += MDfields[ind]
        posRead += MDfields[ind]
        ind += 1 # Next posRef should be a base/string of bases
        while (count > 0):
          for base in MDfields[ind]:
            res.append((posRef,posRead,base+">"+seq[posRead]))
            posRef+= 1
            posRead+= 1
            count-= 1
          ind += 1 # Next posRef should be a number
          if (count > 0):
            if MDfields[ind] > count:
              posRef += count
              posRead += count
              MDfields[ind] -= count
              count = 0
            else:
              count -= MDfields[ind]
              posRef += MDfields[ind]
              posRead += MDfields[ind]
              ind += 1
      if (ind < len(MDfields)) and (MDfields[ind] == 0):
        ind += 1
    elif (op == 4): 
      posRead += count
    elif (op == 5):
      continue
    elif (op == 1):
      ins="INS:"+seq[posRead:posRead+count]
      res.append((posRef,posRead,".>"+ins))
      #res.append((posRef,posRead,("INS",seq[posRead:posRead+count])))
      posRead += count
    elif (op == 2) and (ind % 2 == 1):
      #res.append((posRef,posRead,("DEL",MDfields[ind][1:])))
      dels="DEL:"+MDfields[ind][1:]
      res.append((posRef,posRead,dels+">."))
      posRef += count
      ind += 1
    else:
      sys.stderr.write("Cigar case not implemented: (%d,%d) (%s) [%d (%s),%d]\n"%(op,count,str(cigar),ind,str(MDfields),posRef))
      #sys.exit()
      return []
  #print MD,cigar,MDfields,res,posRef
  return res

import pysam
import numpy as np

samfile = pysam.AlignmentFile(sys.argv[1], "rb",check_sq=False)
n=int(sys.argv[2])
for read in samfile:
	if not read.is_unmapped:
		cigar = read.cigar
		MD = read.get_tag("MD")		
		NM = read.get_tag("NM")
		UMI = "."
		if (len(cigar)>1) and (cigar[0][0]==4) and (cigar[0][1]<=8):
			UMI = read.seq[0:cigar[0][1]]+UMI
		if (len(cigar)>2) and (cigar[-1][0]==4) and (cigar[-1][1]<=8):
			UMI = UMI + read.seq[(-1*cigar[-1][1]):]
		a = parseMDwithCigar(MD,cigar,read.seq,read.pos)
		a=list(filter(lambda x:x[0]<n,a)) 
#		print("%d\t%d\t%s\t%s"%(read.pos,len(read.seq),UMI,parseMDwithCigar(MD,cigar,read.seq)))
#		print(parseMDwithCigar(MD,cigar,read.seq))
#		print("%d\t%d\t%s\t%d\t%s"%(read.pos,len(read.seq),UMI,NM,a))
#		print("%d\t%d\t%s\t%d\t%s"%(read.pos,RL,UMI,NM,a))
		print("%d\t%d\t%s\t%d\t%s"%(read.pos,aln_length(cigar),UMI,NM,list(map(lambda x:"%d:%s"%(x[0],x[2]),a))))
