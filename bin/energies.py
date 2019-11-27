#!/usr/bin/env python
"""
This program set up and run MCCE step 3 in multiple threads.

Usage examples:

1. Run step 3 with 6 threads
    energies.py mp 6

2. Run step 3 to recreate head3.lst
    energies.py refresh

3. Run partial conformers. conformer 1 to 100
    energies.py conf 1, 100

Assuming the underlying c code will:
1. Write only opp files that are calculated (so that multi-threads won't overwrite each other). Dummies and
un-calculated conformers don't have an opp file.
2. When write opp file, no previous run is loaded.
3. A separate procedure will read all opp files and compose head3.lst at the end.
"""