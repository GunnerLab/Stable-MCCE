# MCCE tools reference
<small><i>Page last updated on: {{ git_revision_date }}</i></small>

## Improving and helping runs

---
### energies.py
*Run step3 with multiple threads.*

**Syntax:**
```
energies.py [-h] [-p processes] [-r] [-c start end] [-e /path/to/mcce]
```

**Optional arguments:**

Argument | Description
---|---
-h, --help  |      show this help message and exit
-p processes |     run mcce with number of processes, default to 1
-r          |     refresh opp files and head3.lst without running delphi
-c start end |     starting and ending conformer, default to 1 and 9999
-e /path/to/mcce | mcce executable location, default to "mcce"

**Examples:**

```
energies.py -p 6 
```
This command runs step 3 in 6 threads.

```
energies.py -r 
```
This command runs step 3 just to refresh head3.lst and vdw energy.

```
energies.py -c 1 100 -p 4 
```
This command calculates conformer from 1 to 100 in 5 threads.

**Note**

After this command, the lines:
```
1       delphi start conformer number, 0 based            (PBE_START)
99999   delphi end conformer number, self included        (PBE_END)
```
in run.prm file might be altered.

---
### getpdb
*Download a pdb file by its PDB ID.*

**Syntax:**
```
getpbd pdbid [saved name]
```

**Example:**
```
getpdb 1akk
Inquiring the remote file 1AKK.pdb ...
Saving as 1AKK.pdb ...
Download completed. 
```

Or 
```
getpdb 1akk prot.pdb
Inquiring the remote file 1AKK.pdb ...
Saving as prot.pdb ...
Download completed. 
```

to save as specified name.


---
### translate_step2.py
*Translate old atom names in step2_out to new PDB v3 names.*

**Syntax:**
```
translate_step2.py step2_out.pdb
```

Hydrogen names in PDB v3 are significantly different from PDB v2. Due to the randomness in step 2, 
it is sometimes desired to preserve step2_out.pdb. This tool translates the names so that step2_out.pdb so that it 
can be used by new mcce.

This program prints the translated content on screen. So the user is responsible to make backup of the step2_out.pdb 
file and redirect the outcome to a new step2_out.pdb.  
  
**Example:**
```
cp step2_out.pdb step2_out.pdb.bak
translate_step2.py step2_out.pdb.bak > step2_out.pdb
```


## Result Analysis
---
### state.py
*Analyze how a microstate energy is composed.*

**Syntax:**
```
state.py
```

**Required files:**

 * run.prm
 * head3.lst
 * extra.tpl for scaling factors that are not equal to 1.
 * energies/\*.opp for pairwise interactions 
 * microstate file named "states"
 
The microstate file is in a fixed format. The first line is like:
```
T=298.15, ph=7.0, eh=0.0
```

It specifies the the environment variables which affects the energy calculation. Variables are separated by ",
". The variable names are **not** case sensitive. 

 * T: Temperature. Room temperature is 298.15 K
 * ph: pH value.
 * eh: Eh value, the electric potential for potential titration.

The rest of the lines are microstates. One line per state. Each line consists the index number of conformers (0 
based). Numbers can be separated by either comma or space.
 
**Example:**

Make a states file:

```

```

Run 
```
state.py
```

---
### mfe.py
*Analyze ionization energy by mean field approach.*

**Syntax:**
```
state.py Residue_ID pH [pairwise interaction cut_off]
```

* Residue_ID: Residue ID as in pK.out
* pH: pH value at which ionization energy is calculated. Environment pH is a factor of AH <=> A- + H+ reaction. Also the calculation will use the conformation of other residues at this pH. 
* cut_off: Report only pairwise interaction bigger than this cut off value.

**Required files:**

 * run.prm
 * head3.lst
 * extra.tpl for scaling factors that are not equal to 1.
 * energies/\*.opp for pairwise interactions 

 
**Example:**

Check titration result in pKa.out:
```
cat pK.out
... 
 ASP-A0052_        2.275
...
```

ASP-A0052 is the residue_ID. Its calculated pKa is 2.275. At this point, the free energy of reaction from ASP neutral to ASP ionized should be close to 0.

Run 
```
jmao@jmao-desktop ~/3wum $ mfe.py ASP-A0052_  2.275 0.1
Residue ASP-A0052_ pKa/Em=2.275
=================================
Terms          pH     meV    Kcal
---------------------------------
vdw0        -0.01   -0.71   -0.02
vdw1         0.00    0.20    0.00
tors        -0.08   -4.47   -0.10
ebkb        -0.82  -47.40   -1.11
dsol         1.22   71.06    1.67
offset      -0.62  -36.17   -0.85
pH&pK0       2.48  143.66    3.38
Eh&Em0       0.00    0.00    0.00
residues    -2.17 -126.03   -2.96
*********************************
TOTAL        0.00    0.14    0.00  sum_crg
*********************************
ASNA0044_   -0.65  -37.99   -0.89    0.00
ARGA0045_   -0.11   -6.45   -0.15    1.00
ASNA0046_   -0.61  -35.50   -0.83    0.00
ASPA0048_    0.33   19.44    0.46   -0.58
GLNA0057_    0.15    8.84    0.21    0.00
ASNA0059_   -0.69  -40.19   -0.94    0.00
ARGA0061_   -0.25  -14.74   -0.35    1.00
ASPA0066_    0.14    8.38    0.20   -0.41
ARGA0112_   -0.18  -10.37   -0.24    1.00
=================================
```

You can do mfe calculation at pH other than mid-point.



## Parameter file preparation

## Connecting to other programs
