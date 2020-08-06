# Calculate pKas of Lysozyme 
<small><i>Page last updated on: {{ git_revision_date }}</i></small>

## Background

Lysozyme is a small enzyme that dissolves bacterial cell walls thus kills bacteria. It was discovered as the first 
antibiotic to inhibit bacteria growth in food, before penicillin.  

![lysozyme](https://cdn.rcsb.org/images/rutgers/wu/3wum/3wum.pdb-500.jpg)

Crystal structure of hen egg white lysozyme - PDB ID 4LZT

The experimental pKas of all Lysozyme residues is provided at the bottom of this tutorial. Here we will focus on two residues in the active site that have perturbed pKas. Residue GLU 35 and ASP 52 are in the active sites. The pKas of these two residues play an important role in the enzyme's activity. In order for lysozyme to attack the glucose molecule of the substrate, GLU needs a high pKa and ASP 52 a low pKa. GLU 35 acts as a proton donor, cutting the glucose with protonation of the glycosidic oxygen and a deprotonated ASP 52 stabilizes the highly charged intermediate, making the reaction easier.

[Jens Erik Nielsen and J. Andrew McCammon, Protein Sci. 2003 Sep; 12(9): 1894–1901](https://www.ncbi.nlm.nih
.gov/pmc/articles/PMC2323987/)  

## Prepare the calculation

After the program is installed and the environment variable PATH is configured (see [link](quick.md#configure-environment) for details), make a working directory and go to the directory. MCCE will generate intermediate files and result files in 
the current directory, so it's best to prepare one directory for each calculation on one structure.
 
```
$ mkdir 4lzt
$ cd 4lzt
```
 
Then download pdb file 4LZT from rscb.org to the working directory:
```
jmao@jmao-desktop ~/4lzt $ getpdb 4lzt
Saving as 4LZT.pdb ...
Download completed. 
```

**Explicit water in pdb file.** We suggest removing explicit waters from the input pdb file. The protonation states and pKas with and without are similar, but without waters the calculation is much faster. 

```
grep -v HOH 4LZT.peb > 4LZT_noHOH.pdb 
```

## Run 4 steps of MCCE

MCCE expects to be run in a single directory with only one calculation in that directory.  There are 4 steps.  Each step expects the output from the previous steps.  You can rerun later steps without running from the beginning.  E.g. if you have run steps 1-4 you can rerun step 4.  This will overwrite the original step4. 

### Step 1. Convert PDB file into MCCE PDB
Step 1. Check for inconsistencies between PDB file and MCCE topology files.  Missing sidechain atoms will be added back.  Residues that are not in the MCCE lexicon will be flagged.  Explicit waters with greater than 5% surface exposure will be stripped off.  

Chain termini are separated to be NTR and CTR and are treated as protonatable.  Other changes in residue naming are made using name.txt.

```
$ step1.py 4LZT.pdb
```

This command generates step1_out.pdb which is required for step 2.


### Step 2. Make side chain conformers

This step makes alternative side chain locations and ionization states.

```
$ step2.py
```

This command generates step2_out.pdb which is required of step 3.

If you want to know the help information and other options of this command:
```
$ step2.py -h
```

### Step 3. Make energy table
This step calculates conformer self energy and pairwise interaction table.

```
$ step3.py
```

This command generates opp files under energies/ folder and file head3.lst which are required of step 4.

If you want to know the help information and other options of this command:
```
$ step3.py -h
```

### Step 4. Simulate a titration with Monte Carlo sampling
This step simulates a titration and writes out the conformation and ionization states of each side chain at various conditions.

```
$ step4.py --xts
```

* fort.38. is the primary output file. It gives the average occupancy of each conformer.
* sum_crg.out give the net charge on all residues.  Residues where all conformers have a zero charge are not reported in this file.
* pK.out calculates the best single pKa for each residue given the change in charge with pH in sum_crg.out.  


## Results
The pKa report is in file pK.out.

```
$ cat pK.out
```

From the result, GLU 35 (pKa = 5.13) has a higher pKa than ASP 52 (pKa = 3.33).

To analyze the ionization energy of an ionizable residue at the mid point pH=5.13 with pairwise cutoff 0.1:
```
$ mfe.py ASP-A0052_ -c 0.1
```

To analyze the ionization energy of this residue pH 7 with pairwise cutoff 0.1:
```
$ mfe.py ASP-A0052_  -p 7 -c 0.1
```

A more detailed explanation of mfe.py program can be found [here](tools.md#mfepy)

## Benchmark pKas for Lysozyme
There are 20 experimentally measured pKas in hen white lysozyme.  

* Bartik, K., C. Redfield, and C.M. Dobson, Biophys J, 1994. 66(4): p. 1180-4.
* Kuramitsu, S. and K. Hamaguchi, J Biochem (Tokyo), 1980. 87(4): p. 1215-9.
* Takahashi, T., H. Nakamura, and A. Wada, Biopolymers, 1992. 32: p. 897-909.

These pKa values have been used to benchmark MCCE and other programs that calculate pKas.  For example: 

* Sham, Y. Y., I. Muegge, and A. Warshel. 1999. Simulating proton trans- locations in proteins: probing proton transfer pathways in the Rhodobacter sphaeroides reaction center. Proteins. 36:484–500. 
* You, T. J., and D. Bashford. 1995. Conformation and hydrogen ion titration of proteins: a continuum electrostatic model with conformational flexi- bility. Biophys. J. 69:1721–1733. 
* Antosiewicz, J., J. A. McCammon, and M. K. Gilson. 1996. The determi- nants of pKa’s in proteins. Biochemistry. 35:7819–7833. 

### pKas of residues in Lysozyme

| Residue | pKa |
|---|---|
|LYS 1 | 10.8 |
|GLU 7 | 2.85 |
|LYS 13| 10.5 |
|HIS 15| 5.36 |
|ASP 18| 2.66 |
|TYR 20| 10.3 |
|TYR 23| 9.8  |
|LYS 33| 10.4 |
|GLU 35| 6.2  |
|ASP 48| 1.6  |
|ASP 52| 3.68 |
|ASP 66| 0.9  |
|ASP 87| 2.07 |
|LYS 96| 10.8 |
|LYS 97| 10.3 |
|ASP 101|4.08 |
|LYS 116|10.2 |
|ASP 119| 3.2 |
|CTR | 2.75 |
