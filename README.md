# mcce-develop
Current development version of mcce.

The task of this version of mcce is to consolidate the changes accumulated over the years since MCCE2 2009 paper.

## To install this development version

After downloading cloning this repository, go to the mcce directory, and run:

```
make clean
make
```

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

To test run the mcce code, refer [this page](https://sites.google.com/site/mccewiki/install-mcce)


## What is on the roadmap?
1. Start from Junjun's version to
  * remove dependency on gdbm,
  * use text version of energy lookup table,
  * include gfortran version of delphi
2. Merge with mcce 2.5.1 to
  * add step 5, Monte Carlo post-process
  * enable APBS solver
  * display H bond network
3. Upgrade to mcce 3.5 to
  * add an ele column to opp file for native conformers at single-conformer condition
  * separate the calculation of conformer VDW from step 3
  * atom based VDW output
4. Integrate delphi as function call   

Installation:
glib.h error

Install glib:
```sudo apt-get install libglib2.0-dev```

See the command options:
```pkg-config --cflags --libs glib-2.0```

Add the above output to makefile as compiler option (the next is an example):
```-I/usr/include/glib-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include  -lglib-2.0```
