# What needs to be done to have a working version (without step6)

## General
 - [ ] I’d like to have fort.38 in step 4 and sum_crg.out; pK.out etc in step 5 **alrady have it in step 4, script to do it again**
 - [X] The corrected sum_crg (using head3.lst for charge)
 - [ ] the mfe analysis in pK.out (and an ability to change the mfe point from pK to a pH by retuning step 5 **mfe.py**
 - [ ] Change in chemical potential of a ligand chosen in run.prm, **temp fix script, automated**

## Error checking
 - [ ] The program stops and returns an error message **From users' feedback**
 - [ ] more than the maximum number of conformers **Why the error happens and not printed out clearly**
 - [ ] the number of conformers in all input files (head3.lst; opp etc) don’t match **Check the consistency**
 - [ ] write out for Asp, Glu, Arg, Lys: which ones are <90% ionized and the total number (if you have a lot of neutral Arg the run is likely to be bad) **A tool or step 4, print out abnoramlies, chi2, n etc**

## Improvement 
 - [?] Do not pass 0 radius H to delphi to see if this fixes delphi surface error.
 - [X] Fix delphi run time error.
 - [X] Strip down and reorganize run.prm; 
 - [ ] have full.prm (with all possibilities) on Wiki. **Good idea, group options that default to a pre-defined choices**
 - [X] One run.prm with a toggle for run quick or default (rather than 2 different basic input files) **Questionaire for determining the run.prm**

## Open questions (are these done or for next version?)

 - [X] Proton naming; 
 - [X] new tpl files, **take both tpl, ftpl**
 - [ ] Python scripts to show titration; **doable** 
 - [ ] Phi Map instructions **Ask Dyvia? On a single delphi input file from microstate, or most occupied confomer structure**
