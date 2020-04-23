# What needs to be done to have a working version (without step6)

## General
 - [ ] I’d like to have fort.38 in step 4 and sum_crg.out; pK.out etc in step 5
 - [ ] The corrected sum_crg (using head3.lst for charge)
 - [ ] the mfe analysis in pK.out (and an ability to change the mfe point from pK to a pH by retuning step 5
 - [ ] Change in chemical potential of a ligand chosen in run.prm

## Error checking
 - [ ] The program stops and returns an error message
 - [ ] more than the maximum number of conformers
 - [ ] the number of conformers in all input files (head3.lst; opp etc) don’t match
 - [ ] write out for Asp, Glu, Arg, Lys: which ones are <90% ionized and the total number (if you have a lot of neutral Arg the run is likely to be bad)

## Improvement 
 - [ ] Strip down and reorganize run.prm; have full.prm (with all possibilities) on Wiki
 - [ ] One run.prm with a toggle for run quick or default (rather than 2 different basic input files)

## Open questions (are these done or for next version?)

 - [X] Proton naming; 
 - [ ] new tpl files
 - [ ] Python scripts to show titration; 
 - [ ] Phi Map instructions
