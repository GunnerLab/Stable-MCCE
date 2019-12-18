Version 1.0

This readme file outlines the basic functions of the Evolutionary Packing & Sampling algorithm.

I) Linux Operating System Requirements:
You must ensure your running kernel is at least 2.6.x.x.
A computing node with at least 4GB Ram is highly recommended.
A computing node with an Intel CPU is highly recommended.

II) MCCE Compilation:
If your cluster or computing node has an Intel CPU, we urge you to download the free Linux Intel Compiler to compile MCCE, instead of the GNU GCC Compiler.  If you checkout the source code from our CVS server, there are 2 Makefile files you will need to modify. I have already made the changes to use either the GCC or ICC Intel Compiler. One simply has to comment / uncomment the appropriate lines.  Performance increase by using the Intel compiler can vary between 20%-50%.

1st Makefile: ~./MCCE2.x/Makefile
2nd Makefile: ~./MCCE2.x/lib/Makefile
Intel Compiler Download Link: http://software.intel.com/en-us/articles/non-commercial-software-download/
(Make sure you download the 64-bit version)

III) Execution:
Our algorithm can operate in 2 different modes:

a) Regular sidechain packing and sampling: Regular approach. It generates sidechain rotamers at the specified angular increment, packs, generates a molecular ensemble distribution within the specified delta E energy range and returns the results to MCCE for electrostatics and monte-carlo sampling steps.

b) Spherical focus sidechain packing and sampling: This approach allows a user to focus the sidechain packing and sampling around a specified residue sequence number.  This is especially useful if one wishes to study mutants in a large protein like PSII for example.  All residues which have a least 1 on its backbone or sidechain atom fall within the probe radius is included.  The minimum probe radius is 20.0Angstroms.  The algorithm the strips away all other residues' sidechains and retains only their backbones. Packing and sampling proceeds as normal after this point.  This approach also allows users to generally increase the number of rotations to perform per rotatable bond for all sidechains.

For optimal results in either of these 2 approaches, we highly recommend using 12 rotations per bond for the generation of rotamers.

After generating sidechain rotamers for either case, the algorithm saves these rotamers in a binary file: "output_resume_patch0001.pdb".
This file is used so that next time the user wishes to work from the exact same pdb structure, our algorithm will simply load this rotamer file, instead of generating rotamers each time for the same protein.  Make sure to delete this file found in the folder ~./<root mcce>/pdb_patches/ if you wish to re-generate rotamers or use a modified pdb structure.

Once this file is written, the algorithm begins its packing and sampling. At the end of this process, another file "ga_output_patch0000" is written in the ~./<root mcce>/pdb_patches/ folder.  This file allows the the algorithm to avoid performing packing&sampling again by simply loading the results of a previous optimization run.  This can be useful if your step2 crashed while insert ionization states due to a power outage for example.  Make sure to delete this file if you want a brand new GA run.

IV) Algorithm Parameters
To use the Evolutionary Packing & Sampling Method, you must ensure to have all of these parameters defined in your run.prm file.

---------------------------------------------------------------------------------------------------------
param  1: 1        sidechain optimization: 0=mcce, 1=P.Comte's GA     (SIDECHAIN_OPT)

param  2: -1       Random seed value, use -1 for random value         (GA_SEED)
param  3: 5000     Population size of each GA: use up to 20k          (GA_POP_SIZE)
param  4: 0.1      Mutation rate for each GA                          (GA_MUTATION_RATE)
param  5: 0.0      Migration rate for each GA                         (GA_MIGRATION_RATE)
param  6: 0.9      Crossover rate for each GA                         (GA_CROSSOVER_RATE)
param  7: 30       Max nb. of rnd cut points for crossover            (GA_RANDOM_CUT_POINTS)

param  8: 10       Nb. of GA rounds without convergence check         (GA_PHASE)
param  9: 5 	   Nb. of GA rounds with convergence check	      (GA_SHIFT)

param 10: 15.0     Distribution Center			              (GA_DIST_CENTER)
param 11: 2.0      solutions within EPS*DIST.CENTER in bucket         (GA_DIST_CENTER_EPS)

param 12: 400000   Max nb. protein conformations in final bucket      (GA_MAX_BUCKET_POP)
param 13: 5.0      Local Residue lowest occupied energy cutoff        (GA_RESIDUE_MIN_ENERGY_CUTOFF)
param 14: 0.0      Occupied conformers in MIN_ENRG_CUTOFF are deleted (GA_OCCUPANCY_CUTOFF)
param 15: 0.7      Pop. %  for dE calculation        		      (GA_DELTA_E)

param 16: -1	   Res.Seq. number to center at			      (GA_SPHERE_FOCUS_RESID)
param 17: 20.0	   Default probe radius to scan			      (GA_SPHERE_PROBE_RADIUS)

param 18: 0        GA output format: 0 display screen, 1 file         (GA_OUTPUT)
---------------------------------------------------------------------------------------------------------

Parameter  1:  A 0 indicates the use of MCCE's random packing and sampling technique. A 1 indicates to use the Evolutionary optimization.0

Parameter  2:  The random number generator seed value. A user may wish to reproduce results of packing and sampling, or randomize it each time.

Parameter  3:  The constant number of protein conformations to keep at each generation of the evolutionary algorithm.  For larger proteins (300+ residues), we recommend increasing the population size by a factor of 1,000. A Rough rule of thumb: bigger the population, the better will the sidechain packing be.

Parameter  4:  The mutation rate of protein conformations.  A small value below 0.15 is usually preferred.

Parameter  5:  The migration rate will copy a protein conformation as is without change into the next generation's population.  Disabled by default.

Parameter  6:  The crossover rate dictates how often 2 selected protein conformations should shuffle genetic material to produce 2 new offspring conformations.

Parameter  7:  At each crossover operation between 2 protein conformations, a random number of cutting points for the shuffled exchange of genes is generated.  This parameter controls the maximum number of such cutting points.  This may help a better packing as well.

Parameter  8 & 9:  These are used to control the rate of automatic convergence.  The algorithm performs GA_PHASE number of generations (rounds) without checking whether a new best packing was found.  Subsequently, GA_SHIFT number of consecutive generations (rounds) are performed to check if a new best packing was found, gradually increasing the convergence.  Once 100% convergence is reached, the algorithm performs the evolutionary sampling.  If a new best packing is found during one GA_SHIFT round, the convergence check is reset to 0%.

Parameter  10:  The distribution center (in delta E Angstrom units) around which to generate the molecular ensemble of the protein.

Parameter  11:  Sampling range = Distribution center * epsilon.  Protein conformations which have a delta E to the best protein conformation within the sampling range are included in the molecular ensemble of the protein. Note that duplicate conformations are allowed here since we are after a statistical distribution.

Parameter  12:  The maximum number of protein conformations in the sampling bucket for the generation of the molecular ensemble. One may increase this value to 1,000,000 or 2,000,000, provided there is enough available RAM on the computing node.

Parameter  13:  After we obtain a molecular ensemble, we filter each residue's rotamers.  Rotamers of a residue that fall within 5Kcal/mol of the lowest pairwise van der Waals (VDW) interaction the residue occupied while generating the molecular ensemble, are kept.  Others are simply removed.

Parameter  14:  One may also wish to filter residues' rotamers based on the occupancy of a rotamer.  We recommend not using this filtering method, although it is available for experimentation.

Parameter  15:  The percentage of the population of protein conformations to take into account well update the average population delta E value.  Do not modify this parameter unless you have a clear and explicit complete understanding of the algorithm.


Parameter  16 & 17:  A -1 indicates that the algorithm performs a regular sidechain packing & sampling.  A different value - the residue sequence number in a protein - specifies the location at which to perform a spherical focus with specified probe radius.  The algorithm includes all residues within the sphere and then the regular packing & sampling resumes.  Other residues' sidechains outside the sphere are removed, though we keep their backbones.

Parameter  18:  Output the algorithm's output to screen or to a file.

V) Contact Information
Please report any bugs or request additional features to:

Pascal Comte
PhD Candidate
Department of Computer Science
Memorial University
St. John's, Newfoundland, Canada
e-mail: pcomte at mun.ca