#!/usr/bin/env python

import os
from datetime import datetime

record_file = "run.prm.record"

# write the entries to file run.prm
def export_runprm(runprm):
    lines = []
    for key in runprm:
        line = "%-20s    (%s)\n" % (runprm[key], key)
        lines.append(line)
    open("run.prm", "w").writelines(lines)
    return

# update the step section with new entries to run.prm.recorded
def record_runprm(runprm, section_id):
    toplines = []
    newid = "%s: %s\n" % (section_id, datetime.now().strftime("%Y/%m%d, %H:%M:%S"))
    bodylines = [newid]
    buttomlines = []
    for key in runprm:
        line = "%-20s    (%s)\n" % (runprm[key], key)
        bodylines.append(line)


    found = False
    buttom = False
    if os.path.exists(record_file):
        lines = open(record_file).readlines()
        for line in lines:
            if not found:
                if line.find(section_id) == 0:
                    found = True
                else:
                    toplines.append(line)
            else:
                if buttom:
                    buttomlines.append(line)
                elif line.find("#STEP") == 0:
                    buttom = True
                    buttomlines.append(line)

    open(record_file, "w").writelines(toplines+bodylines+buttomlines)
    #print(buttomlines)

    return

