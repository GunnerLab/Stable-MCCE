#!/usr/bin/env python

import logging

class Env:
    def __init__(self):
        # hard code values
        self.tpl = {}
        return

    def load_ftpl(self, fname):
        """Load a ftpl file"""

        # default string values, or defined as below
        float_values = ["EXTRA", "SCALING"]
        int_values = []

        logging.info("   Loading ftpl file %s" % fname)
        lines = open(fname).readlines()
        for line in lines:
            line = line.split("#")[0]
            fields = line.split(":")
            if len(fields) != 2:
                continue

            key_string = fields[0].strip()
            keys = key_string.split(",")
            keys = [x.strip().strip("\"") for x in keys]
            keys = [x for x in keys if x]
            keys = tuple(keys)

            value_string = fields[1].strip()
            if keys[0] in float_values:
                self.tpl[keys] = float(value_string)
            elif keys[0] in int_values:
                self.tpl[keys] = int(value_string)
            else:
                self.tpl[keys] = value_string

        return

class Protein:
    """Protein structure"""
    def __init__(self):
        return

