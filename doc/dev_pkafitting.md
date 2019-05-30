# Reculate sum_crg.out and pK.out with python

<<<<<<< HEAD
Python IDE:


## What do we want?
=======
## Session 1:
### What do we want?
>>>>>>> 10360f56f1b8d19aaa5255548676ce6bd4b59d3c
* Correct charge in sum_crg.out
* Multiple midpoint in pK.out

### What do we need?
* occupancy in fort.38
* charge states in head3.lst

### Program
* read fort.38
* read head3.lst
* match fort.38 and head3.lst
* get conformer list
* get charge states
* get occupancy table
* multiply occupancy table with charge states to have a charge table

### Data structure demonstrated in this script:
* variable (type)
* list (append, pop, delete, remove, slice, and index)
* class (structure in C)
* numpy array (elements have to be the same type, array addition and multiplication)

## Session 2:
### Program flow
* group conformers to residues
* pick out charged residues
* write out sum_crg.out
* loop over charged residues and fit their titration curves with scipy

### Data structure demonstrated in this script:
* tuple
* dictionary
* class (to group data)
* SciPy optimize module

## Session 3:
### Write code that could be reused
* function
* class (to embed method = function in class)
* multiple code files (global variables and cross file imports)

### Visiualize your data with matplotlib
* basic concepts and common practice
* where to start - gallary from documentation

