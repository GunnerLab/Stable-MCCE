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

## Parameter file preparation

## Connecting to other programs
