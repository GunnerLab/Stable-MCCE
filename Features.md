# Features to be implemented

## Must have

* [ ] Use "H" in PDB v3 nameing rules. - Junjun
   * [x] add element field for ATOM structure file mcce.h
   * [x] function to determin element name
   * change the code to determin H atom based on element field
* [ ] Keep original H position for conformer 0 if H coordinates were provided by PDB. Good for connection with MD.
   * only the completly filled conformer is labeled as native
* [ ] Step 3 pairwise interaction benchmark.
* [ ] Microstate output - Xiu Hong
* [ ] Hydrogen network report - Xiu Hong
* [ ] phi map
   * most occupied confomer pdb file
   * delphi calculation
   * apply potential on B factor field for visulization
   
## Wanted
 * [ ] Chemical potential titration
 * [ ] Parallel threads
    * Monte Carlo
    * Delphi in step 3
 * [ ] Other PB solvers
    * new delphi
    * apbs
  
