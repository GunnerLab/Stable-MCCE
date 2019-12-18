The Text MFE Tool

"mfe.py" reports the ionization energy of a group at specified pH/eh point.

Syntax:
   mfe.py residue_name pH_value [threshehold]

   residue_name: the residue identifier in file "pK.out".
   pH_value:     the solution condition where the mfe analysis is performed.
   threshold:    pairwise interaction threshold of display in pH unit. This
                 argument is optional with a defualt value of 0.

Remarks:
   The command must be invoked in the working directory because the following
   files in the current directory will be used:
      run.prm
      head3.lst
      pK.out
      fort.38
      energies/

   If the titration is an eh titration, determined by the "run.prm", pH_value
   refers to an eh value, but the threshold is still in pH unit. A pH unit
   equals 58 meV at room temperature.

   When the specified titration point (pH_value) is not the simulation point, a
   linear interpolation will be performed.

Author:
   Junjun Mao (jmao@sci.ccny.cuny.edu)
