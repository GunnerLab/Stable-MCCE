# MCCE tools reference
<small><i>Page last updated on: {{ git_revision_date }}</i></small>

## Improving and helping runs

---
### energies.py
*Run step3 with multiple threads.*

**Syntax:**
energies.py [-h] [-p processes] [-r] [-c start end] [-e /path/to/mcce]

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


**Example:**


---
### translate_step2.py
*Translate old atom names in step2_out to new PDB v3 names.*

**Syntax:**


**Example:**



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

**Example:**


## Parameter file preparation

## Connecting to other programs
