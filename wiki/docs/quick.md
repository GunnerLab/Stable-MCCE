# MCCE quick start
<small><i>Page last updated on: {{ git_revision_date }}</i></small>


## Prerequisites

### Compilers: C, gfortran

The C and gfortran compilers usually come with your operating system. If not, use the package manager to install.

!!! note "On Ubuntu Linux:"
    ```
    $ sudo apt-get update
    $ sudo apt-get upgrade
    $ sudo apt-get install build-essential
    $ sudo apt-get install gfortran
    ```

!!! note "On Mac OS X:"
    For Mac users, the gcc compiler comes with Xcode from the App store. Once Xcode is installed, install command line tools. This provides gcc.

    Gfortran is avalaible from https://github.com/fxcoudert/gfortran-for-macOS/releases, find your Mac OS X version and install appropriate gfortran package.


### Python and modules
We need Python3 and optionally these modules: numpy, scipy, matplotlib, pygraphviz, pandas, xlrd, and openpyxl

* Install Miniconda Python3: https://docs.conda.io/en/latest/miniconda.html
* After installing Miniconda, install these optional modules for data analysis:
```
conda install numpy scipy matplotlib pygraphviz pandas xlrd openpyxl
```


## Get the code

### Download the latest code
Under a terminal window, run:

    git clone https://github.com/GunnerLab/Stable-MCCE.git

### Compile the code:
The above command creates a directory named Stable-MCCE, enter the directory and compile:

```
$ cd Stable-MCCE
$ make clean
$ make
```

!!! warning "Compiling on Mac OS X:"
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


### Configure environment
Find the path of Stable-MCCE installation directory:
```
(base) jmao@pc:~/projects/Stable-MCCE$ pwd
/home/jmao/projects/Stable-MCCE
```

In my case the installation directory is **/home/jmao/projects/Stable-MCCE/**

Add the executable to your path:
```
export PATH=/home/jmao/projects/Stable-MCCE/bin:$PATH
```

Remember to replace the installation directory **/home/jmao/projects/Stable-MCCE/** with your own installation path.

Also put this line at the end of .bashrc file under your home directory so that the environment is properly set every time you open a terminal window.


## Run MCCE

### Prepare a working directory:
```
$ mkdir test_lysozyme
$ cd test_lysozyme
```

### Get a pdb file
```
$ getpdb 1dpx
```

You now have a pdb file 1DPX.pdb in the working directory.

The simplest way to run mcce is do these four steps:

### Step 1 convert PDB file into MCCE PDB
This step proof reads the structure file and cuts terminal residues and complex cofactors into smaller ones if necessary.
```
$ step1.py 1DPX.pdb
```

### Step 2 make side chain conformers
This step makes alternative side chain locations and ionization states.
```
$ step2.py
```

### Step 3 make energy table
This step calculates conformer self energy and pairwise interaction table.
```
$ step3.py
```

### Step 4 Simulate a titration with Monte Carlo sampling
This step simulates a titration and writes out the conformation and ionization states of each side chain at various conditions.
```
$ step4.py
```

## Notes

* For more detailed command usages, use "-h" switch in each command above.
* Some steps take hours to finish, so it is recommended to run at the background. For example:
```
step3.py > run.log &
```