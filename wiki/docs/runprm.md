# How to customize run.prm, the run configuration file?
<small><i>Page last updated on: {{ git_revision_date }}</i></small>

The mcce program reads file run.prm from the working directory and runs according to its instructions.

The two ways that run mcce are:

1. customize run.prm and call mcce directly
2. run step1.py, step2.py, step3.py, are step4.py as wrappers of mcce.

## run.prm format
The valid parameters are defined as a key-value pair, defined by a line.

At the end of the line, the key string word is enclosed in a pair of parenthesis. The key string must be one word.

If a key is detected from a line, the first word of the line is its value.

For example:
```
f        step 1: pre-run, pdb-> mcce pdb                    (DO_PREMCCE)
```
The key is "DO_MCCE" and the value is "f".

Lines other than parameter lines are considered as comment lines, therefore ignored.

## Most commonly edited entries

### Input and output
---
```
prot.pdb                                                    (INPDB)
```
Th parameter "INPDB" specifies a pdb file as the initial input structure.
-------------------------------------------------------------------------

### Major mcce steps
MCCE can be divided into four steps.

blockdiag {
   A [label = "Step 1\nProof read structure'];
   B [label = "Step 2\nMake conformers"];
   C [label = "Step 3\nMake energy table"];
   D [label = "Step 4\nSimulated titration"];

   A -> B -> C -> D;
}

The steps have keys:

 * DO_PREMCCE
 * DO_ROTAMERS
 * DO_ENERGY
 * DO_MONTE

The value 't' is for true (do), and 'f' is for false (not do).
```
Steps:
f        step 1: pre-run, pdb-> mcce pdb                    (DO_PREMCCE)
f        step 2: make rotatmers                             (DO_ROTAMERS)
f        step 3: do energy calculations                     (DO_ENERGY)
f        step 4: monte carlo sampling                       (DO_MONTE)
```

The later step depends on the out put of previous step.


## Less edited entries

Step 2 and step 4 have more flexibility and are customized more often than others.

### General
---
```
4.0      Protein dielectric constant for DelPhi             (EPSILON_PROT)
```
Dielectric constant of protein. Usually epsilon 8 is used for small soluble protein (~100 residues), and 4 is used for big trans-membrane protein (> 200 residues).
This value is used by step 2 to optimize rotamers and step 3 to compute energy table.
---

```
/home/mcce/Stable-MCCE/extra.tpl                              (EXTRA)
```
This file contains entries that overwrite those defined in ftpl and tpl files in param/ directory, and extral controls over the program.
The extra controls are in two categories:

1. offset energy of ionizable residues. These energy corrects the residue solution pKa values based on previous benchmarks.
2. The scaling factor of van der Walls energy term. The recommended values are the inverse of dielectric constant defined by EPSILON_PROT.
---

```
/home/mcce/Stable-MCCE/name.txt MCCE renaming rule.           (RENAME_RULES)
```
The renaming rule can

* divide a big group into separate ionizable groups.
* put ligands and cofactor as one group.

This is done by renaming the atom and residue.

---


### Step 2, Rotamer making

---
```
f        Use control file "head1.lst" for rotamer making    (ROT_SPECIF)
```
This file enables site specific rotamer making rules. When set as "t", it overwrites the rules defined by run.prm.
---

```
t        Do swap (stereo isotope)                           (ROT_SWAP)
```
Make stereo isotope conformers for residues like GLN and ASN.

---

```
t        Do rotate?                                         (PACK)
```
Rotate around rotatable bonds?

---

```
6        number of rotamers in a bond rotation              (ROTATIONS)
```
Number of rotamers in a bond rotation.

---

### Step 4, Monte Carlo sampling
---
```
ph       "ph" for pH titration, "eh" for eh titration       (TITR_TYPE)
```
The type of simulation, pH or Eh.
---

```
0.0      Initial pH                                         (TITR_PH0)
```
The starting pH, or the pH of Eh titration.
---
```
1.0      pH interval                                        (TITR_PHD)
```
The interval of pH titration
---
```
0.0      Initial Eh                                         (TITR_EH0)
```
The starting Eh, or the Eh of pH titration.
---
```
30.0     Eh interval (in mV)                                (TITR_EHD)
```
The interval of Eh titration.
---
```
15       Number of titration points                         (TITR_STEPS)
```
Number of titration points.
---
```
f        Output Microstate from standard monte carlo        (MS_OUT)
f        Output readable micristate file from standard monte carlo (RE_MS_OUT)
```
Write out microstates or not.
-----------------------------

## System specific entries

These entries are tailored to the system and installation.

### System and installation location

---
```
/home/mcce/Stable-MCCE                                      (MCCE_HOME)
```
The MCCE installation root. It is required by the program to find relevant parameters.
---

```
/home/mcce/Stable-MCCE/bin/delphi DelPhi executable           (DELPHI_EXE)
```
PB solver delphi location.
--------------------------
```
/tmp     delphi temporary file folder, "/tmp" uses node     (PBE_FOLDER)
```
Location of PB solver output. This should be a folder with write permission and large enough Sometimes to store GBs of files.
-----------------------------------------------------------------------------------------------------------------------------
```
new.tpl     local parameter file for unrecognized res       (NEWTPL)
```
Program mcce creates a temporary tpl file whenever it couldn't solve a new cofactor. The presumed parameters are in this file.
------------------------------------------------------------------------------------------------------------------------------

```
debug.log                                                   (DEBUG_LOG)
```
Debug log file name.
--------------------

### Detail levels of program output
---
```
t       Minimize output files                               (MINIMIZE_SIZE)
```
If this is "f", mcce will write more debug information.

---

```
t       Print trace file                                    (DO_TRACE)
```
Creates a record with MCCE revision number, time and date of run, and all user options set in run.prm and appends it to (or creates) the run.trace file.

---
```
t        clean up delphi focusing directory?                (DELPHI_CLEAN)
```
Flag of cleaning delphi temporary folder. Set to 't' only for debugging purpose.
--------------------------------------------------------------------------------


## Other entries


### Step 1
---
```
t        Label terminal residues?                           (TERMINALS)
```
If this flag is on 't', mcce program will identify the C terminus and N terminus of the polypeteides based on the geometry. Terminal residues have extra ionizable sites than residues having peptide bonds. One needs to be careful if there are missing piece of peptide in the given structure.

---

```
0.1      cut off water if % SAS exceeds this number         (H2O_SASCUTOFF)
```
The water molecules are stripped of if the surface exposure is over this cutoff. This a recursive procedure which means after some waters are removed, the new surface exposure is recaculated and the process id repeated.

---

```
2.0      distance limit of reporting clashes                (CLASH_DISTANCE)
```
The threshold to report abnormal distance found between non-bonded atoms.

---

```
t		ignore hydrogens in input structure					(IGNORE_INPUT_H)
```
Delete and rebuild input H atoms in step 1 and step 2.

---

### Step 2:
---
```
t        Rebuild sidechain based on torsion minima          (REBUILD_SC)
```
Rebuild side chain.
-------------------

```
f        Do swing?                                          (SWING)
```
Swing conformers deviate from their parent conformer by a small degree.
-------------------------

```
10.0     phi in degrees of swing                            (PHI_SWING)
```
Swing angle.
------------

```
1.00     SAS cuttoff of making fewer rotamers               (SAS_CUTOFF)
```
Surface residues get fewer rotamers. There is another mechanism to search the most exposed rotamer to compensate this.
----------------------------------------------------------------------------------------------------------------------

```
10.0     Cutoff of self vdw in kcal/mol                     (VDW_CUTOFF)
```
If a rotamer has self energy, comparing with the most favorable one, over this cutoff, it will not ne kept.
-----------------------------------------------------------------------------------------------------------

```
5000     number of repacks                                  (REPACKS)
```
The rotamer optimization starts from random initial state, then minimizes the energy on each side chain. The number of restarts of this process is rapacks.
-----------------------------------------------------------------------------------------------------------------------------------------------------------

```
0.03     occupancy cutoff of repacks                        (REPACK_CUTOFF)
```
The occupancy of side chain conformer lower than this value will not be kept.
-----------------------------------------------------------------------------

```
t        h-bond directed rotamer making                     (HDIRECTED)
```
Do search on hydrogen bond.
---------------------------

```
1.0      threshold for two conformers being different       (HDIRDIFF)
```
Only one of the similar conformers is kept. Similarity os measured by Root Mean Square Deviation.
-------------------------------------------------------------------------------------------------

```
36       Limit number of the H bond conformers              (HDIRLIMT)
```
The number limit on H bond directed conformers.
-----------------------------------------------

```
f        Do relaxation on water                             (RELAX_WAT)
```
Move water molecules to relax the distance crash.
-------------------------------------------------

```
3.2      if distance between water and heavy atom is smaller than this, move water away  (WATER_RELAX_THR)
```
The threshold of moving water.
------------------------------

```
10       number of cycles to relax rotamer clash            (HV_RELAX_NCYCLE)
```
Heavy atom rotamer relaxation by MD. Number of global restarts (cycles).
------------------------------------------------------------------------

```
10       time (in fs) for each relaxation iteration         (HV_RELAX_DT)
```
Heavy atom rotamer relaxation time.
-----------------------------------

```
10       number of iterations                               (HV_RELAX_NITER)
```
Number of local restarts.
-------------------------

```
2.       relax rotamer if vdw interaction is bigger than this (HV_RELAX_VDW_THR)
```
Lower end threshold of relaxing rotamers, Van der Waals energy value in kcal/mol.
---------------------------------------------------------------------------------

```
5.       not relax rotamer if heavy atom vdw is bigger than this (HV_RELAX_HV_VDW_THR)
```
High end threshold of relaxing rotamers, Van der Waals energy value in kcal/mol.
--------------------------------------------------------------------------------

```
-1.0     relax if electrostatic interaction is more favorable than this, and the charged groups are close  (HV_RELAX_ELEC_THR)
```
Relax heavy atoms when they show dominantly favorable electrostatic interaction.
-------------------------------------------------------------------------------

```
0.1      threshold used to define a charged atom, see hv_relax_elec_dist_thr (HV_RELAX_ELEC_CRG_THR)
```
Atoms with charge over this value are considered charged atoms.
---------------------------------------------------------------

```
3.0      only relax electrostatically favorable pairs that has charged atom distance shorter than this (HV_RELAX_ELEC_DIST_THR)
```
Relax heavy atoms when they show close proximity.
-------------------------------------------------

```
20.      scaling factor for torsion during relaxation       (HV_TORS_SCALE)
```
Scaling factor for torsion energy when searching torsion minima.
----------------------------------------------------------------

```
10000    maximun n steps of shake                           (HV_RELAX_N_SHAKE)
```
MD shake parameter.
-------------------

```
0.001    maximun allowance for bond length to deviate       (HV_RELAX_SHAKE_TOL)
```
MD bond length detaion parameter.
---------------------------------

```
1.0      constraint distance                                (HV_RELAX_CONSTRAINT)
```
MD atom constrain distance.
---------------------------

```
20.      constraint force for atom stay in original position (HV_RELAX_CONSTRAINT_FRC)
```
MD atom constraint force.
-------------------------

```
f        include neighbors when doing relaxation (expensive) (HV_RELAX_INCLUDE_NGH)
```
Include neighbor in MD relaxation.
----------------------------------

```
4.0      threshold for include as a neighbor during relaxation (HV_RELAX_NGH_THR)
```
Define a neighboring side chain.
--------------------------------

```
t        Do relaxation on hydrogens                         (RELAX_H)
```
Optimize hydrogen position.
---------------------------

```
-5.0     Energy threshold for keeping the conformer         (RELAX_E_THR)
200      Loop over N local microstates                      (RELAX_NSTATES)
6        default number of hydroxyl positions               (RELAX_N_HYD)
5.       do not relax hydrogen if vdw bwt two sidechain conformer bigger than this    (RELAX_CLASH_THR)
1.0      phi for each step of relaxation                    (RELAX_PHI)
300      Maximum number of steps of relaxation              (RELAX_NITER)
0.5      Torque threshold for keep relaxing                 (RELAX_TORQ_THR)
```
These parameters are for MD relaxation.
---------------------------------------

```
0.02     Last pruning threshold for conformers              (PRUNE_THR)
0.5      Pruning cutoff of RMSD                             (PRUNE_RMSD)
1.0      Pruning cutoff of eletrostatic pairwise            (PRUNE_ELE)
8.0      Pruning cutoff of vdw pairwise                     (PRUNE_VDW)
```
To reduce the number of final conformers, similar conformers are grouped into one. The similarity is defined by these parameters.
---------------------------------------------------------------------------------------------------------------------------------

```
0        maximum number of conformer per residue ( 0=unlimited ) (NCONF_LIMIT)
```
Set a limit to final conformers of a residue.
---------------------------------------------

### step 3:
---
```
f        Use SAS + Coulomb's law for ele interaction        (QUICK_ENERGIES)
```
Use coulomb's law instead of solving PB equation.
-------------------------------------------------

```
80.0     Solvent dielectric constant for DelPhi             (EPSILON_SOLV)
65       Grids in each DelPhi                               (GRIDS_DELPHI)
2.0      The target grids per angstrom for DelPhi           (GRIDS_PER_ANG)
1.4      Radius of the probe                                (RADIUS_PROBE)
2.0      Ion radius                                         (IONRAD)
0.15     Salt                                               (SALT)
```
PB solver parameters.
----------------------------

```
f        skip delphi in step3 (DEBUG option)                (SKIP_ELE)
```
This is a debug option.
-----------------------

```
f        Reassign charge and radii before delphi(if not true then use value in PDB file) (REASSIGN)
```
Step 3 uses atom charges from step2_out.pdb. Use this to overwrite.
-------------------------------------------------------------------

```
f        Recalculate torsions energy when write out head3.lst (RECALC_TORS)
```
Recalculate torsion energy.
---------------------------

```
1.7      defalut van der Waals radius, for SAS              (DEFAULT_RADIUS)
0.5      factor to 1-4 LJ potontial (1.0 is full)           (FACTOR_14LJ)
6.0      dielectric constant for Coulomb's law              (EPSILON_COULOMB)
```
More parameters used by the prigram in step 1, 2 and 3.
-------------------------------------------------------

```
1         delphi start conformer number, 0 based            (PBE_START)
99999     delphi end conformer number, self included        (PBE_END)
```
Start and end of PB solution on conformers.
-------------------------------------------


###  step 4:
---
```
f        Average the pairwise, "f" uses the smaller         (AVERAGE_PAIRWISE)
```
When loading energy from step 3, average the pairwise interaction.
------------------------------------------------------------------

```
20.      Warning Threshold of difference in pairwise       (WARN_PAIRWISE)
```
Issue a warning message if asymmetric pairwise interaction is detected.
-----------------------------------------------------------------------

```
5.0      Big pairwise threshold to make big list            (BIG_PAIRWISE)
```
Residues with big interactions may get synchronized flips. This is the cut off value.
-------------------------------------------------------------------------------------

```
-1       Random seed, -1 uses time as random seed           (MONTE_SEED)
```
Random seed.
------------

```
298.15   Temperature                                        (MONTE_T)
3        Number of flips                                    (MONTE_FLIPS)
100      Annealing = n_start * confs                        (MONTE_NSTART)
500      Equalibration = n_eq * confs                       (MONTE_NEQ)
0.001    Cut-off occupancy of the reduction                 (MONTE_REDUCE)
```
Monte Carlo simulation parameters.
----------------------------------

```
###  Only for Junjun's monte carlo
2        Number of independent monte carlo sampling         (MONTE_RUNS)
2000     Sampling = n_iter * confs                          (MONTE_NITER)
20000    Trace energy each MONTE_TRACE steps, 0 no trace    (MONTE_TRACE)
1000000  Maximum microstates for analytical solution        (NSTATE_MAX)
f        Specify mfe point, f=pKa/Em                        (MFE_POINT)
-1.0     MFE cutoff(Kcal), default 0.5                      (MFE_CUTOFF)
```
More Monte Carlo simulation parameters.
---------------------------------------

```
###  Only for Yifan's monte carlo
t        Using Yifan's monte carlo                          (MONTE_ADV_OPT)
f        Using format from old version                      (MONTE_OLD_INPUT)
5000     Min Sampling = n_iter * confs                      (MONTE_NITER_MIN)
-1       Max Sampling = n_iter * confs(-1 stop @ converge)  (MONTE_NITER_MAX)
10000    Number of iterations each cycle                    (MONTE_NITER_CYCLE)
1000     niter * n_conf check convergence                   (MONTE_NITER_CHK)
-1       Number of the reduce steps(-1 stop @ converge)     (MONTE_N_REDUCE)
0.01     Threshhold for convergence                         (MONTE_CONVERGE)
f        Calculate free energy                              (MONTE_DO_ENERGY)

298.15   Starting temperature for annealling                (ANNEAL_TEMP_START)
0        Number of steps of annealling                      (ANNEAL_NSTEP)
1000     Number of iterations for each annealling step      (ANNEAL_NITER_STEP)
```
More Monte Carlo simulation parameters.
---------------------------------------

### IPECE parameters
---
```
f add neutral atoms to simulate a membrane slab 	    (IPECE_ADD_ME)
33. the thichness of the membrane to be add                 (IPECE_MEM_THICKNESS)
M the chain ID of the add atoms				    (IPECE_MEM_CHAINID)
```
This is for creating a simulated membrane.
------------------------------------------