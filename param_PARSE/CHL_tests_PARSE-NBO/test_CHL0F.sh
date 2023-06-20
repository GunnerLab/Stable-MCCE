#!/bin/bash

rm -rf reorient_CHL0F.txt
touch reorient_CHL0F.txt
grep 'CONNECT, " MG ", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " CHA", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " CHB", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " CHC", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " CHD", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " HHB", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " HHC", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " HHD", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " N1A", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " C1A", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " C2A", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " C3A", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " C4A", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " CMA", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " CAA", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " CBA", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " CGA", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " O1A", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " O2A", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, "HMA1", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, "HMA2", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, "HMA3", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " H2A", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " H3A", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, "HAA1", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, "HAA2", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, "HBA1", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, "HBA2", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " N1B", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " C1B", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " C2B", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " C3B", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " C4B", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " CMB", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " CAB", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " CBB", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " HBB", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " OMB", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt ## HMB1 to OMB
grep 'CONNECT, " HMB", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt ## HMB2 to HMB

grep 'CONNECT, "HBB1", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, "HBB3", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " N1C", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " C1C", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " C2C", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " C3C", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " C4C", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " CMC", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " CAC", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " CBC", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, "HMC1", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, "HMC2", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, "HMC3", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, "HAC1", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, "HAC2", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, "HBC1", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, "HBC2", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, "HBC3", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " N1D", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " C1D", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " C2D", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " C3D", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " C4D", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " CMD", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " CAD", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " OBD", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " CBD", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " CGD", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " O1D", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " O2D", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " CED", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, "HMD1", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, "HMD2", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, "HMD3", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, " HBD", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, "HED1", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, "HED2", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CONNECT, "HED3", CHL0F:' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt

echo  >> reorient_CHL0F.txt

grep 'CHARGE, CHL0F, " MG ":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " CHA":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " CHB":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " CHC":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " CHD":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " HHB":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " HHC":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " HHD":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " N1A":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " C1A":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " C2A":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " C3A":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " C4A":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " CMA":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " CAA":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " CBA":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " CGA":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " O1A":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " O2A":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, "HMA1":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, "HMA2":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, "HMA3":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " H2A":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " H3A":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, "HAA1":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, "HAA2":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, "HBA1":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, "HBA2":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " N1B":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " C1B":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " C2B":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " C3B":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " C4B":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " CMB":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " CAB":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " CBB":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " HBB":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " OMB":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt ## HMB1 to OMB
grep 'CHARGE, CHL0F, " HMB":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt ## HMB2 to HMB

grep 'CHARGE, CHL0F, "HBB1":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, "HBB3":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " N1C":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " C1C":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " C2C":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " C3C":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " C4C":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " CMC":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " CAC":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " CBC":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, "HMC1":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, "HMC2":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, "HMC3":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, "HAC1":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, "HAC2":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, "HBC1":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, "HBC2":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, "HBC3":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " N1D":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " C1D":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " C2D":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " C3D":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " C4D":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " CMD":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " CAD":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " OBD":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " CBD":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " CGD":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " O1D":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " O2D":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " CED":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, "HMD1":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, "HMD2":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, "HMD3":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, " HBD":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, "HED1":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, "HED2":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'CHARGE, CHL0F, "HED3":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt 

echo  >> reorient_CHL0F.txt

grep 'RADIUS,  CHL0F,  " MG ":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " CHA":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " CHB":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " CHC":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " CHD":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " N1A":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " C1A":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " C2A":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " C3A":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " C4A":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " CMA":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " CAA":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " CBA":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " CGA":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " O1A":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " O2A":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " N1B":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " C1B":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " C2B":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " C3B":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " C4B":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " CMB":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " CAB":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " CBB":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " OMB":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " N1C":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " C1C":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " C2C":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " C3C":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " C4C":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " CMC":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " CAC":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " CBC":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " N1D":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " C1D":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " C2D":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " C3D":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " C4D":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " CMD":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " CAD":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " OBD":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " CBD":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " CGD":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " O1D":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " O2D":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
grep 'RADIUS,  CHL0F,  " CED":' CHL_QMhybrid_head+dipole.ftpl >> reorient_CHL0F.txt
