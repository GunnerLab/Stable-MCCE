# Calculate Em of heme in Cytochrome
<small><i>Page last updated on: {{ git_revision_date }}</i></small>

## Background

Cytochrome c is a small protein that transports electron in mitochondria to facilitate the synthesis of ATP. It redox potential plays an imprtant role in its function. The regulation of the Cytochrome C redox potential can be explained by Continuum Electrostatic analysis.

*[Junjun Mao, Karin Hauser, and M. R. Gunner, How Cytochromes with Different Folds Control Heme Redox Potentials, Biochemistry 2003, 42, 33, 9829–9840](https://pubs.acs.org/doi/full/10.1021/bi027288k)*


Cytochrome C from E. Caballus (PDB ID 1AKK) has an experimental Em of 260 mV.

## Prepare the calculation

After the program is installed and execution path is configured (see [link](quick.md) for details),
make a working directory and go to the working directory. MCCE will generate intermediate files and result files in
the current directory, so it's best to prepare one directory for calculation on one structure.

```
$ mkdir 1akk
$ cd 1akk
```

Then download pdb file 1AKK to the working directory:
```
$ getpdb 1akk
Inquiring the remote file 1AKK.pdb ...
Saving as 1AKK.pdb ...
Download completed. 
```

## Run 4 steps of MCCE

### Step 1. Convert PDB file into MCCE PDB
This step proof reads the structure file and cut terminal residues and complex cofactors into smaller ones if necessary.

The heme in cytochrome C has two ligands HIS18 and MET80. They behave differently than HIS and MET so we have to rename them. Copy name.txt from the MCCE directoy to working directory.
```
$ cp /path/to/mcce/name.txt ./
```

Your /path/to/mcce can be found by running

```
$ var=`which mcce`; echo ${var%/*/*}
```

Now that you have a name.txt in your working directory, add these two lines to ./name.txt:
```
*****HIS**  18  *****HIL**  18
*****MET**  80  *****MEL**  80
```

Then run step 1:
```
$ step1.py -u RENAME_RULES=./name.txt 1AKK.pdb
```

The reason we specify the location is to force mcce to use our own renaming rules instead of the default one.

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
This setp simulates a titration and write out the conformation and ionization states of each side chain at various conditions.

```
$ step4.py -i 0 -d 60 -t eh
```

* The occupancy table is in file fort.38.
* The net charge is in file sum_crg.out
* Eh is in file pK.out

## Results
The pKa report is in file pK.out.

```
$ cat pK.out
```

Your will see the calculated Eh for heme is 247 mV

To analyze the ionization energy of heme at the midpoint:
```
$ mfe.py HEM+A0105_
```
