# MCCE quick start
<small><i>Page last updated on: {{ git_revision_date }}</i></small>

## Get the code

### To download the latest code:
Under a terminal window, run:

    git clone https://github.com/GunnerLab/Stable-MCCE.git

### Compile the code:
This creates a directory Stable-MCCE, enter the directory and compile:

You need C and Fortran to compile the code.

!!! Note
    For Mac users, gcc compiler comes with Xcode from app store. Once Xcode is installed, install command line tools. This provides gcc.

    Gfortran is avalaible from https://github.com/fxcoudert/gfortran-for-macOS/releases, find your Mac OS X version and install appropriate gfortran package.


```
cd Stable-MCCE
make clean
make
```

Find the path of Stable-MCCE installation directory:
```
(base) jmao@pc:~/projects/Stable-MCCE$ pwd
/home/jmao/projects/Stable-MCCE
```

**In my case the {/path/to/mcce/} is /home/jmao/projects/Stable-MCCE**

Add the executable to your path:
```
export PATH={/path/to/mcce/}bin:$PATH
```

Also put this line at the end of .bashrc file under your home directory so that the environment is properly set every time you open a terminal window.

!!! Warning
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

## Python requirement
A significant part of the tools are written in Python3. You will need Python installed on your system. Example of installing Python and reuired modules:

### Install Python package manager Miniconda
Miniconda: https://docs.conda.io/en/latest/miniconda.html

After installing Miniconda, install these modules:
```
conda install numpy scipy pandas
```

## Run MCCE is simplest way

### Prepare a working directory:
```
mkdir test_lysozyme
cd test_lysozyme
```

### Get a pdb file
```
getpdb 1dpx
```

You now have a pdb file 1DPX.pdb in the working directory.

The simplest way to run mcce is do these four steps:

### Step 1 convert PDB file into MCCE PDB
This step proof reads the structure file and cut terminal residues and complex cofactors into smaller ones if necessary.
```
step1.py 1DPX.pdb
```

### Step 2 make side chain conformers
This step makes alternative side chain locations and ionization states.
```
step2.py
```

### Step 3 make energy table
This step calculates conformer self energy and pairwise interaction table.
```
step3.py
```

### Step 4 Simulate a titration with Monte Carlo sampling
This setp simulates a titration and write out the conformation and ionization states of each side chain at various conditions.
```
step4.py
```

## Notes

* For more detailed command usages, use "-h" switch in each command above.
* Some steps take hours to finish, so it is recommended to run at the background. For example ```step3.py > run.log &```