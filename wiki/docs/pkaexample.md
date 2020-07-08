# Calculate pKas of Lysozyme 
<small><i>Page last updated on: {{ git_revision_date }}</i></small>

## Background

Lysozyme is a small enzyme that dissolves bacterial cell walls thus kills bacteria. It was discovered as the first 
antibiotic to inhibit bacteria growth in food, before penicillin.  

![lysozyme](https://cdn.rcsb.org/images/rutgers/wu/3wum/3wum.pdb-500.jpg)

Crystal structure of hen egg white lysozyme - PDB ID 3WUM

The residue GLU 35 and ASP 52 have been identified as two active sites. The pKas of these two residues play an important
 role in the enzyme's activity. In order for lysozyme to attack the glucose molecule of the substrate, 
 GLU has high pKa and ASP 52 has low pKa. This way GLU 35 acts as a proton donor, cuts the glucose with protonation of the 
 glycosidic oxygen and a deprotonated ASP 52 stabilizes the highly charged intermediate, making the reaction easier.

[Jens Erik Nielsen and J. Andrew McCammon, Protein Sci. 2003 Sep; 12(9): 1894–1901](https://www.ncbi.nlm.nih
.gov/pmc/articles/PMC2323987/)  

## Prepare the calculation

After the program is installed and the environment variable PATH is configured (see [link](quick.md#configure-environment) for details), 
make a working directory and go to the working directory. MCCE will generate intermediate files and result files in 
the current directory, so it's best to prepare one directory for calculation on one structure.
 
```
$ mkdir 4lzt
$ cd 4lzt
```
 
Then download pdb file 4LZT to the working directory:
```
jmao@jmao-desktop ~/4lzt $ getpdb 4lzt
Inquiring the remote file 4LZT.pdb ...
Saving as 4LZT.pdb ...
Download completed. 
```

## Run 4 steps of MCCE

### Step 1. Convert PDB file into MCCE PDB
This step proof reads the structure file and cut terminal residues and complex cofactors into smaller ones if necessary.
```
$ step1.py 4LZT.pdb
```

This command generates step1_out.pdb which is required of step 2.

If you want to know the help information and other options of this command:
```
$ step1.py -h
```

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

* The occupancy table is in file fort.38.
* The net charge is in file sum_crg.out
* pKas are in file pK.out

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
