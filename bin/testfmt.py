#!/usr/bin/env python

"""
Convert unformated pdb fort.13 to formated file
typedef struct {
   int header;
/* jmao: IDL of gfortran needs a 4-byte int to store the length of a record */
   float x;
   float y;
   float z;
   float rad;
   float crg;
   int trailer;
} UNF;
"""

import sys
import struct

struct_fmt = '=ifffffi'
struct_len = struct.calcsize(struct_fmt)
struct_unpack = struct.Struct(struct_fmt).unpack_from

if __name__ == "__main__":
    if len(sys.argv) < 2:
        fname = "fort.13"  # default is fort.13
    else:
        fname = sys.argv[1]

    with open(fname, "rb") as f:
        while True:
            record = f.read(struct_len)
            if not record:
                break
            record_txt = struct_unpack(record)
            print(record_txt)


