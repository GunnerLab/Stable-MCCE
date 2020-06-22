## Compute hydrogen bond network
This is a tutorial to calculate the hydrogen bond network using MCCE 
and visualizing hydrogen bond network using Cytoscape. 

Steps:
1.  First, run MCCE step 6 if you already run from step 1 to step 4 otherwise 
 you need to run from step 1 to step 6. 
step2_out.pdb, ms.dat, ms_gold (old MCCE) are required to run step 6.

Default hydrogen bond definition in run.prm is:

 ![Screenshot](run.prm.png)

You can change hydrogen bond criterion accordingly. 

2. Output file after step 6:

      hb.dat, hah.txt and hb.txt

### Result Analysis:

#### Cytoscape visualization:
We are using Cytoscape for visualizing hydrogen bond networks. 
Download and install Cytoscape (https://cytoscape.org/)

####Input file preparation for Cytoscape:
a. For single hb.txt file:
For the single hydrogen bond network file obtain from the crystal structure, 
you can visualize by directly by opening hb.txt file using Cytoscape 
and selecting different layout.

b. To merge multiple hb.txt obtained by running MCCE from Molecular Dynamics snapshots

To combine multiple hb.txt files obtain from MCCE and see water connected behavior:

Scripts need to have: 

`jhead.h, jhead-test2.cpp, jlu_new-cai.cpp, Makefile`
 and input files: `hb.txt`, `Residues_list.lst`  in the same directory.

First rename hb.txt files into 1.dat, 2.dat .... etc. Suppose if you have 5 hb.txt file then 
rename to 1.dat, 2.dat, 3.dat, 4.dat and 5.dat .
`Residues_list.lst` has information of residues interested for hydrogen bond analysis.

Residues_list.lst looks like this:

![Screenshot](res_lst.png)


a. Open `jhead.h` and here you can change `static int cutoff = 4;` for different water number bridging.
Here 4 is  up to 4 water molecules are allowed to bridge between two residues.
For example: If you want to see up to  two water bridging hydrogen bond connections between two residues, 
then you can change 4 to 2.

b. Open Makefile and you can change `a-4w.out` in makefile:
Here for up to 4 water, we write  `a-4w.out`. Suppose for  two water: change `a-4w.out` to `a-2w.out`

c. Type `make` to compile script:
This will make `a-4w.out`  file. Then run this file by using command: _./a-4w.out_ 

Here output files obtain from MCCE: `out.dat, out_opt.dat, color_map.dat`

Open `out.dat` file using the Cytoscape and play with different layout.














