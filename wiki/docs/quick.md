# MCCE quick start
<small><i>Page last updated on: {{ git_revision_date }}</i></small>

## Get the code

### To download the latest code:

git clone https://github.com/GunnerLab/Stable-MCCE.git

### Compile the code:

Compile:

```
make clean
make
```


Add the executable to your path:

```
export PATH={/path/to/mcce/}bin:$PATH
```

{/path/to/mcce/} is the path to the MCCE installation directory. 

**Troubleshooting:**

On Mac OS X, an explicit memory free for tree is not supported. You may have to change this subroutine in file lib/db.c
```C
/* release database memory */
void free_param() {
   tdestroy(param_root, free);
   return;
}
```
to
```C
/* release database memory */
void free_param() {
   return;
}
```

If you see glib.h error:

Install glib:

```sudo apt-get install libglib2.0-dev```

See the command options:

```pkg-config --cflags --libs glib-2.0```

Add the above output to makefile as compiler option (the next is an example):

```-I/usr/include/glib-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include  -lglib-2.0```


## Prepare a working directory

### Working directory:
```
mkdir test_lysozyme
cd test_lysozyme
```

### Get a pdb file
```
getpdb 1dpx
```

You now have a pdb file 1DPX in the working directory.

### run.prm
MCCE requires a file run.prm to guide the run.
```
cp {/path/to/mcce/}run.prm.quick ./run.prm
```

{/path/to/mcce/} is the path to the MCCE installation directory. Edit the following lines in run.prm:

---
```
prot.pdb                                                    (INPDB)
```
to 
```
1DPX.pdb                                                    (INPDB)
```

---
```
f        step 1: pre-run, pdb-> mcce pdb                    (DO_PREMCCE)
f        step 2: make rotatmers                             (DO_ROTAMERS)
f        step 3: do energy calculations                     (DO_ENERGY)
f        step 4: monte carlo sampling                       (DO_MONTE)
f        step 6: analysis                                   (DO_ANALYSIS)
```
to 
```
t        step 1: pre-run, pdb-> mcce pdb                    (DO_PREMCCE)
t        step 2: make rotatmers                             (DO_ROTAMERS)
t        step 3: do energy calculations                     (DO_ENERGY)
t        step 4: monte carlo sampling                       (DO_MONTE)
f        step 6: analysis                                   (DO_ANALYSIS)
```

---
```
/home/mcce/Stable-MCCE/extra.tpl                              (EXTRA)
/home/mcce/Stable-MCCE/name.txt MCCE renaming rule.           (RENAME_RULES)
```
to 
```
{/path/to/mcce/}extra.tpl                                   (EXTRA)
{/path/to/mcce/}name.txt MCCE renaming rule.                (RENAME_RULES)
```

---
```
/home/mcce/Stable-MCCE                                        (MCCE_HOME)
```

to
```
{/path/to/mcce/}                                              (MCCE_HOME)
```

---
```
/home/mcce/Stable-MCCE/bin/delphi DelPhi executable           (DELPHI_EXE)
```

to
```
{/path/to/mcce/}bin/delphi DelPhi executable                  (DELPHI_EXE)
```
---

```
/scratch     delphi temporary file folder, "/tmp" uses node     (PBE_FOLDER)
```
to
```
/tmp     delphi temporary file folder, "/tmp" uses node     (PBE_FOLDER)
```

## Run mcce
```
mcce > run.log &
```

The log will be saved in file run.log.


