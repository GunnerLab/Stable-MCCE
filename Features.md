# Features to be implemented

## Must have

* [x] Use "H" in PDB v3 nameing rules. - Junjun
   * [x] add element field for ATOM structure file mcce.h
   * [x] function to determin element name
   * [x] change the code to determin H atom based on element field
* [ ] Keep original H position for conformer 0 if H coordinates were provided by PDB. Good for connection with MD.
   * only the completly filled conformer is labeled as native
* [ ] Step 3 pairwise interaction benchmark.
* [ ] Step 3 multiple threads support
   * [ ] write ele part only in opp at delphi
   * [ ] finish up / merge opp files
* [ ] Microstate output - Xiu Hong
  * [x] upload previous code
  * [ ] change microstate output format
* [ ] Hydrogen network report - Xiu Hong
* [ ] mfe in the analysis report
* [ ] Control of MC: standard, with complete ms out, diff ms out.
* [ ] Tools:
   * [ ] vdw clash for a specific conformer at atom level.
   * [ ] phi map
   * [ ] most occupied confomer pdb file
   * [ ] phi map file
   * [ ] Why ms hb analysis is slow? 

## Wanted
 * [ ] Chemical potential titration
 * [ ] Parallel threads
    * Monte Carlo
    * Delphi in step 3
 * [ ] Other PB solvers
    * new delphi
    * apbs
 * [ ] How does reduced MC improve acceptance?
 * [ ] Improved running energy and Average energy output
 
## TODO next week
 * [ ] Fix Taichi's TYR problem.
 * [ ] benchmark step 3
 * [ ] Native H in step 2
 * [ ] Intruction on cofactor ftpl
 * [ ] Intruction how to run to step 3
