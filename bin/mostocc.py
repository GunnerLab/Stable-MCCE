#!/usr/bin/env python

import sys


if len(sys.argv) < 2:
   print("mostocc.py pH")
   print("pH is the matching exact number in the first line of fort.38. Occupancy varies with titration point, so we have to choose a column.")
   sys.exit(0)

# get occupancy from fort.38 column 8
OccTable = {}
fort38 = open("fort.38").readlines()
pH = sys.argv[1]
header = fort38[0].split()
column = header.index(pH)
for line in fort38:
   # get confID
   fields = line.split()
   if len(fields[0]) == 14:
      confID = fields[0]
      occ = float(fields[column])
      OccTable[confID] = occ

# group into residues
OccInRes = {}
for conf in OccTable.keys():
   resID = conf[:3]+conf[5:11]
   if resID in OccInRes:
      OccInRes[resID].append((conf, OccTable[conf]))
   else:
      OccInRes[resID] = [(conf, OccTable[conf])]


# get max conf
maxConfs = []
for res in OccInRes.keys():
   confs = OccInRes[res]
   maxocc = confs[0]
   for conf in confs:
      if conf[1] > maxocc[1]:
         maxocc = conf
   maxConfs.append(maxocc[0])

# read in a file and keep only the most occupied confs
pdb = open("step3_out.pdb").readlines()
outpdb = []
for line in pdb:
   if len(line) <82: continue
   if line[26] == ' ': iCode = '_'
   else: iCode = line[26]
   confID = line[17:20]+line[80:82]+line[21:26]+iCode+line[27:30]
   if confID in maxConfs or confID[3:5] == 'BK':
      outpdb.append(line)

sys.stdout.writelines(outpdb)
