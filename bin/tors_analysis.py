#!/usr/bin/env python

lines = open("a").readlines()

confs = []
new_tors = {}
for line in lines:
    fields = line.split()
    confID = fields[1]
    tors = fields[11]
    confs.append(confID)
    new_tors[confID] = tors

print("Conformer Tors_old Tors_new")
lines = open("head3.lst").readlines()
for line in lines:
    fields = line.split()
    confID = fields[1]
    tors = fields[11]
    if confID in new_tors:
        print("%s %s %s" % (confID, tors, new_tors[confID]))



