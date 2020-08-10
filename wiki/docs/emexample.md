# Calculate Em of heme in Cytochrome
<small><i>Page last updated on: {{ git_revision_date }}</i></small>

## Background

Cytochrome c is a small protein that transports electrons in mitochondria to facilitate the synthesis of ATP. Its redox potential plays an important role in its function. The regulation of the Cytochrome C redox potential can be explained by Continuum Electrostatic analysis.

*[Junjun Mao, Karin Hauser, and M. R. Gunner, How Cytochromes with Different Folds Control Heme Redox Potentials, Biochemistry 2003, 42, 33, 9829–9840](https://pubs.acs.org/doi/full/10.1021/bi027288k)*


Cytochrome C from E. Caballus (PDB ID 1AKK) has an experimental Em of 260 mV.

## Prepare the calculation

After the program is installed and execution path is configured (see [link](quick.md) for details),
make a working directory and go to the working directory. MCCE will generate intermediate files and result files in
the current directory, so it's best to prepare one directory for calculation on one structure.

```
mkdir 1akk
cd 1akk
```

Then download pdb file 1AKK to the working directory:
```
getpdb 1akk
Saving as 1AKK.pdb ...
Download completed. 
```

## Run 4 steps of MCCE

### Step 1. Convert PDB file into MCCE PDB
This step checks for inconsistencies between PDB file and MCCE topology files.  Missing sidechain atoms will be added back.  Residues that are not in the MCCE lexicon will be flagged.  Explicit waters with greater than 5% surface exposure will be stripped off.  

The heme in cytochrome C has two ligands HIS18 and MET80. They behave differently than HIS and MET so we have to rename them. step1.py can handle HIS, MET, and CYS if they are the ligands to heme. If your heme ligand is not one of these three residues, please let us know.

Run step 1:
```
step1.py 1akk.pdb
```

The reason we specify the location is to force mcce to use our own renaming rules instead of the default one.

This command generates step1_out.pdb which is required of step 2.


### Step 2. Make side chain conformers
This step makes alternative side chain locations and ionization states.

```
step2.py
```

This command generates step2_out.pdb which is required of step 3.

If you want to know the help information and other options of this command:
```
step2.py -h
```

### Step 3. Make energy table
This step calculates conformer self energy and pairwise interaction table.

```
step3.py
```

This command generates opp files under energies/ folder and file head3.lst which are required of step 4.

If you want to know the help information and other options of this command:
```
step3.py -h
```

### Step 4. Simulate a titration with Monte Carlo sampling
This setp simulates a titration and write out the conformation and ionization states of each side chain at various conditions.

```
step4.py -i 0 -d 60 -t eh
```

* The occupancy table is in file fort.38.
* The net charge is in file sum_crg.out
* Eh is in file pK.out

## Results
The pKa report is in file pK.out.

```
cat pK.out
```

Your will see the calculated Eh for heme is 247 mV

To analyze the ionization energy of heme at the midpoint:
```
mfe.py HEM+A0105_
```

A more detailed explanation of mfe.py program can be found [here](tools.md#mfepy)