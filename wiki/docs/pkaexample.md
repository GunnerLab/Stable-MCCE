# Calculate pKas of Lysozyme 
<small><i>Page last updated on: {{ git_revision_date }}</i></small>

## Background

Lysozyme is a small enzyme that dissolves bacterial cell walls thus kills bacteria. It was discovered as the first 
antibiotic to inhibit bacteria growth in food, before penicillin.  

![lysozyme](https://cdn.rcsb.org/images/rutgers/wu/3wum/3wum.pdb-500.jpg)

Crystal structure of hen egg white lysozyme - PDB ID 3WUM

The residue GLU 35 and ASP 52 have been identified as two active sites. The pKas of these two residues play an important
 role in the enzyme's activity. In order for lysozyme to attack the glucose molecule of the substrate, 
 GLU has high pKa and ASP 52 has low pKa. This way GLU 35 acts a proton donor, cuts the glucose with protonation of the 
 glycosidic oxygen and a deprotonated GLU 52 stabilizes the highly charged intermediate, making the reaction easier. 

[Jens Erik Nielsen and J. Andrew McCammon, Protein Sci. 2003 Sep; 12(9): 1894–1901](https://www.ncbi.nlm.nih
.gov/pmc/articles/PMC2323987/)  

## Prepare the calculation

After the program is installed and execution path is configured (see [link](quick.md) for details), 
make a working directory and go to the working directory. MCCE will generate intermediate files and result files in 
the current directory, so it's best to prepare one directory for calculation on one structure.
 
```
mkdir 3wum
cd 3wum
```
 
Then download pdb file 3WUM to the working directory:
```
jmao@jmao-desktop ~/3wum $ getpdb 3wum
Inquiring the remote file 3WUM.pdb ...
Saving as 3WUM.pdb ...
Download completed. 
```

## Configure run.prm for MCCE

MCCE requires a file run.prm to guide the run.

```
cp {/path/to/mcce/}run.prm.quick ./run.prm
```

```{/path/to/mcce/}``` is the path to the MCCE installation directory. The file run.prm defines how MCCE runs. If the 
line ends up with a keyword enclosed in a pair of parenthesis, this line is a definition line, 
with the key being the word in parenthesis and the value being the the first word of this line. 

### The most modified entries

```
==============================================================================
Most modified entries
------------------------------------------------------------------------------
Input and Output:
prot.pdb                                                    (INPDB)

Steps:
f        step 1: pre-run, pdb-> mcce pdb                    (DO_PREMCCE)
f        step 2: make rotatmers                             (DO_ROTAMERS)
f        step 3: do energy calculations                     (DO_ENERGY)
f        step 4: monte carlo sampling                       (DO_MONTE)
f        step 6: Hydrogen bond analysis                     (DO_ANALYSIS)
==============================================================================
```

This line defines the input structure file:
```
prot.pdb                                                    (INPDB)
```

It's a good idea to keep the downloaded pdb file 3WUM.pdb, so we make a copy named prot.pdb:
```
cp 3WUM.pdb prot.pdb
```

Since we are calculating pKas, we need to run from step 1 to 4. Mark these steps "t" for "true".
```
t        step 1: pre-run, pdb-> mcce pdb                    (DO_PREMCCE)
t        step 2: make rotatmers                             (DO_ROTAMERS)
t        step 3: do energy calculations                     (DO_ENERGY)
t        step 4: monte carlo sampling                       (DO_MONTE)
```

### The less modified entries

This section is usually modified only for the first run. It defines the location of the parameters and titration 
conditions.

```
8.0      Protein dielectric constant for DelPhi             (EPSILON_PROT)
/home/jmao/projects/Stable-MCCE/extra.tpl                   (EXTRA)
/home/jmao/projects/Stable-MCCE/name.txt MCCE renaming rule.(RENAME_RULES)

ph       "ph" for pH titration, "eh" for eh titration       (TITR_TYPE)
0.0      Initial pH                                         (TITR_PH0)
1.0      pH interval                                        (TITR_PHD)
0.0      Initial Eh                                         (TITR_EH0)
30.0     Eh interval (in mV)                                (TITR_EHD)
15       Number of titration points                         (TITR_STEPS)
```
 ..
In the above example, I set 

* the dielectric constant 8.0, which is good for pKa calculation. 
* The path to extra.tpl and name.txt is set to my mcce installation directory. 
* the titration type is set to be "ph"
* the titration range from 0 to 14 at step 1, 15 titration points in total
 
### The advanced entries
Some entries here have to be customized for your installation and system.

```
/home/jmao/projects/Stable-MCCE                                        (MCCE_HOME)
```

This should point to where your mcce is installed.


```
/home/jmao/projects/Stable-MCCE/bin/delphi DelPhi executable           (DELPHI_EXE)
```

This points to The location of delphi executable.

```
/tmp     delphi temporary file folder, "/tmp" uses node     (PBE_FOLDER)
```

For some systems, /tmp is preferred temporary directory.

## Run mcce

Make sure mcce is in your executable path. For me, I run 
```
export PATH=/home/jmao/projects/Stable-MCCE/bin:$PATH
```

To set the PATH variable.

To run mcce
```
mcce
```

Ths program takes 30 minutes to 1 hour to finish. So you can run 
```
mcce > run.log &
```

to put the program in background.

## Results
The pKa report is in file pK.out.


