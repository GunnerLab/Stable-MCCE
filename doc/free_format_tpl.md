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

In circumstances that one needs to overwrite previously defined entries instead of appending, one can negate an earlier record by defining a special value "!"
ROTATE, GLU: " CA " - " CB "
ROTATE, GLU: !
ROTATE, GLU: " CB " - " CG "
ROTATE, GLU: " CG " - " CD "



## How to convert mcce tple file to free format:
Run this command to convert mcce format to free format.
>bin/tpl-mcce2free.py tpls/glu.tpl > tpls/glu.ftpl

Run this command to complete VDW parameters of RADIUS records.
>bin/vdw-complete.py tpls/glu.ftpl > tmpfile

>mv tmpfile tpls/glu.ftpl