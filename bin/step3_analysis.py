#!/usr/bin/env python

new_folder = "energies"
old_folder = "energies.bak"

def read_confnames():
    confnames = []
    lines = open("head3.lst").readlines()
    for line in lines[1:]:
        fields = line.strip()
        if len(fields) < 17:
            continue
        confnames.append(fields[1])
    return confnames

def load_old_backbone():
    backbone_ele = {}
    lines = open("head3.lst").readlines()
    for line in lines[1:]:
        fields = line.strip()
        if len(fields) < 17:
            continue
        confid = fields[1]
        epol = float(fields[12])
        backbone_ele[confid] = epol
    return backbone_ele

def load_new_backbone(confnames):
    backbone_ele = {}
    for conf in confnames:
        fpath = "%s/%s.raw" % (new_folder, conf)
        lines = open(fpath).readlines()
        for line in lines:
            fields = line.split()

    return backbone_ele

if __name__ == "__main__":
    confnames = read_confnames()
    old_backbone_ele = load_old_backbone()
    new_backbone_ele = load_new_backbone(confnames)

