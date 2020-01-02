#!/usr/bin/env python

import sys

if __name__ == "__main__":

    lines = open(sys.argv[1]).readlines()
    for line in lines:
        print(line, end="")