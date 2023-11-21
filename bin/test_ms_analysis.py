#!/usr/bin/env python

import sys
import ms_analysis as msa

if __name__ == "__main__":

    print(f"{len(msa.conformers) = }")

    # redo conformer loading with non-standard path:
    h3_path = sys.argv[1]
    new_conformers = msa.read_conformers(h3_path)

    print(f"{len(new_conformers) = }")

    sys.exit(0)


