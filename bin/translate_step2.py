#!/usr/bin/env python

import sys

CONVERT = {"1H   HOH" : " H1  HOH"

}

if __name__ == "__main__":

    lines = open(sys.argv[1]).readlines()
    for line in lines:
        if line[:6] == "ATOM  " or line[:6] == "HETATM":
            if line[12:20] in CONVERT:
                print(line[:12]+CONVERT[line[12:20]]+line[20:].rstrip())

