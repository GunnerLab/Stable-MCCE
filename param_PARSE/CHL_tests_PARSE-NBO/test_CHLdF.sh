#!/bin/bash

rm -rf reorient_CHLdF.txt
touch reorient_CHLdF.txt
grep 'CONNECT, " MG ", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " CHA", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " CHB", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " CHC", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " CHD", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " HHB", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " HHC", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " HHD", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " N1A", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " C1A", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " C2A", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " C3A", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " C4A", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " CMA", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " CAA", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " CBA", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " CGA", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " O1A", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " O2A", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, "HMA1", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, "HMA2", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, "HMA3", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " H2A", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " H3A", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, "HAA1", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, "HAA2", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, "HBA1", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, "HBA2", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " N1B", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " C1B", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " C2B", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " C3B", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " C4B", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " CMB", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " CAB", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " CBB", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " HBB", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " OMB", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt ## HMB1 to OMB
grep 'CONNECT, " HMB", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt ## HMB2 to HMB

grep 'CONNECT, "HBB1", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, "HBB3", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " N1C", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " C1C", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " C2C", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " C3C", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " C4C", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " CMC", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " CAC", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " CBC", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, "HMC1", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, "HMC2", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, "HMC3", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, "HAC1", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, "HAC2", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, "HBC1", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, "HBC2", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, "HBC3", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " N1D", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " C1D", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " C2D", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " C3D", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " C4D", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " CMD", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " CAD", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " OBD", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " CBD", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " CGD", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " O1D", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " O2D", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " CED", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, "HMD1", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, "HMD2", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, "HMD3", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, " HBD", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, "HED1", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, "HED2", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CONNECT, "HED3", CHLdF:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt

echo  >> reorient_CHLdF.txt

grep 'CHARGE, CHLdF, " MG ":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " CHA":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " CHB":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " CHC":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " CHD":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " HHB":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " HHC":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " HHD":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " N1A":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " C1A":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " C2A":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " C3A":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " C4A":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " CMA":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " CAA":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " CBA":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " CGA":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " O1A":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " O2A":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, "HMA1":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, "HMA2":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, "HMA3":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " H2A":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " H3A":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, "HAA1":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, "HAA2":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, "HBA1":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, "HBA2":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " N1B":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " C1B":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " C2B":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " C3B":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " C4B":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " CMB":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " CAB":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " CBB":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " HBB":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " OMB":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt ## HMB1 to OMB
grep 'CHARGE, CHLdF, " HMB":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt ## HMB2 to HMB

grep 'CHARGE, CHLdF, "HBB1":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, "HBB3":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " N1C":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " C1C":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " C2C":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " C3C":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " C4C":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " CMC":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " CAC":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " CBC":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, "HMC1":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, "HMC2":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, "HMC3":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, "HAC1":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, "HAC2":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, "HBC1":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, "HBC2":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, "HBC3":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " N1D":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " C1D":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " C2D":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " C3D":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " C4D":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " CMD":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " CAD":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " OBD":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " CBD":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " CGD":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " O1D":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " O2D":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " CED":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, "HMD1":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, "HMD2":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, "HMD3":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, " HBD":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, "HED1":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, "HED2":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'CHARGE, CHLdF, "HED3":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt 

echo  >> reorient_CHLdF.txt

grep 'RADIUS,  CHLdF,  " MG ":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " CHA":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " CHB":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " CHC":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " CHD":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " N1A":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " C1A":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " C2A":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " C3A":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " C4A":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " CMA":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " CAA":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " CBA":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " CGA":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " O1A":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " O2A":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " N1B":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " C1B":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " C2B":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " C3B":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " C4B":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " CMB":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " CAB":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " CBB":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " OMB":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " N1C":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " C1C":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " C2C":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " C3C":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " C4C":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " CMC":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " CAC":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " CBC":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " N1D":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " C1D":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " C2D":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " C3D":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " C4D":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " CMD":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " CAD":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " OBD":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " CBD":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " CGD":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " O1D":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " O2D":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
grep 'RADIUS,  CHLdF,  " CED":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHLdF.txt
