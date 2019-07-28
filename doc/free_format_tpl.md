# Free Format TPL files

The format of free format tpl file is simple. It is simply a key value pair separated by ":".
Key can be up to three fields, separated by ",". Value is a string, and will be processed by the mcce program in the context of key.
If space is a part of the key or value, it must be quoted inside double quotes. ":" and "," are reserved, and can not be used as part of key and value.
The following is an example:

```CONNECT, " N  ", GLUBK: sp2, " ?  ", " CA ", " H  "```

MCCE reads in lines from ftpl files and save them in a single database. Therefore a ftpl file can define multiple residues and parameter lines of one residue can be stored in multiple files.

If the same key appears more than one time, the values are appended to be the earlier value of the same key. Look at this example:

```
ROTATE, GLU: " CA " - " CB "
ROTATE, GLU: " CB " - " CG "
ROTATE, GLU: " CG " - " CD "
```

is equal to:

```ROTATE, GLU: " CA " - " CB ", " CB " - " CG ", " CG " - " CD "```

In circumstances that one needs to overwrite previously defined entries instead of appending, one can negate an earlier record by defining a special value "!!!"
SCALING, VDW0:     !!!
SCALING, VDW0:     1.0

SCALING, VDW1:     !!!
SCALING, VDW1:     1.0

SCALING, VDW:      !!!
SCALING, VDW:      1.0

SCALING, TORS:     !!!
SCALING, TORS:     1.0

means "No matter what were defined before, delete this entry, then create an entry with new key and value pair."


## How to convert mcce tple file to free format:
Run this command to convert mcce format to free format.
>bin/tpl-mcce2free.py param04/glu.tpl > glu.ftpl

The above output file has most entries converted. Then we can manually compile rxn values in CONFORMER records.manually

First manually replace rxn= with rxn02=, rxn04=, rxn08= depending on the dielectric constant of the tpl file that was converted from.
Since mcce tpl files at different dielectric constants are different only by rxn= values, in the free format tpl files, they can be set in the same CONFORMER line.
```
CONFORMER, GLU01: Em0=   0.0, pKa0=  0.00, ne= 0, nH= 0, rxn02= -6.410, rxn04= -3.100, rxn08= -1.390
```

