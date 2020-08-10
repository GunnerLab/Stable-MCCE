# How to convert an existing mcce tpl file to free format tpl file?
<small><i>Page last updated on: {{ git_revision_date }}</i></small>

There accumulated a big number of mcce tpl files (topology file of a molecule) which defines the conformers, atoms, connectivity and so on. Since mcce tpl file follows strict format, new mcce adopted a [free format](newftpl.md) to save the same information. Follow these steps to convert an existing mcce tpl file to a new free format tpl file.

This example is assuming you convert a mcce tpl file glu.tpl in param04 directory.

## Convert

Run this command to convert mcce format to free format.

```
$ bin/tpl-mcce2free.py param04/glu.tpl > glu.ftpl
```

The above output file has most entries converted. 

## Fix "rxn" record

Then we can manually compile rxn values in CONFORMER records. First manually replace ```rxn=``` with ```rxn02=```, ```rxn04=```, ```rxn08=``` depending on the dielectric constant of the tpl file that was converted from.

```
CONFORMER, GLU01: Em0=   0.0, pKa0=  0.00, ne= 0, nH= 0, rxn02= -6.410, rxn04= -3.100, rxn08= -1.390
```

## Error check

Now find the atom name difference, especially H atoms. PDB is using version 3 names and old mcce is using version 2 names. For GLU, the name comparison can be found at [here](http://ligand-expo.rcsb.org/pyapps/ldHandler.py?formid=cc-index-search&target=glu&operation=ccid).

If you have a pdb file of the residue/ligand, we can test the ftpl file against that pdb file:

```
$ bin/verify_tpl.py ftplfile pdbfile 
```

It will report any atoms that can not be loaded by the converted ftpl file definition.

To change the name in ftpl file, glu.ftp for example:

```
$ sed -i 's/1HB / HB2/g' glu.ftpl 
```

You can verify the ftpl file again after the change of names.


## Store in ```param``` directory

Finally, place the tested ftpl file under MCCE distribution's "param" directory.

The MCCE distribution directory is the ```(HOME)``` line in your run.prm.