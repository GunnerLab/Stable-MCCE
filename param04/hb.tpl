#23456789A123456789B123456789C123456789D123456789
#Edit by Cai-20160701-used in final data calculation.
HDONOR   ASP01       HD1
HDONOR   ASP02       HD2

HDONOR   ARG01      1HH1 2HH1 1HH2 2HH2
HDONOR   ARG02      1HH1 1HH2 2HH2  HE
HDONOR   ARG03      1HH1 2HH1 1HH2  HE
HDONOR   ARG+1      1HH1 2HH1 1HH2 2HH2  HE

HDONOR   GLU01       HE1
HDONOR   GLU02       HE2

HDONOR   HIS01       HE2
HDONOR   HIS02       HD1
HDONOR   HIS+1       HD1  HE2

HDONOR   SER01       HG 
#add acceptor atom for neutral ser
HACCEPT  SER01       OG

HDONOR   THR01       HG1
#add acceptor atom for neutral thr
HACCEPT  THR01       OG1

##add ASN
HDONOR   ASN01      1HD2 2HD2
HACCEPT  ASN01       OD1

##add GLN
HDONOR   GLN01      1HE2 2HE2
HACCEPT  GLN01       OE1

##add CYS
HDONOR   CYS01       HG
HACCEPT  CYS01       SG 
HACCEPT  CYS-1       SG

##add MET
HACCEPT  MET01       SD

HDONOR   TYR01       HH

##add TRP
HDONOR   TRP01       HE1

HDONOR   LYS01      1HZ  2HZ
HDONOR   LYS+1      1HZ  2HZ  3HZ
## add lys neutral donor######
HACCEPT  LYS01       NZ
############################

HDONOR   PAA01       H1A
HDONOR   PAA02       H2A

HDONOR   PDD01       H1D
HDONOR   PDD02       H2D
################################################
# hydrogen bond acceptor #######################
################################################
HACCEPT  ARG01       NE
HACCEPT  ARG02       NH1
HACCEPT  ARG03       NH2

HACCEPT  ASP01       OD1  OD2
HACCEPT  ASP02       OD1  OD2
HACCEPT  ASP-1       OD1  OD2

HACCEPT  GLU01       OE1  OE2
HACCEPT  GLU02       OE1  OE2
HACCEPT  GLU-1       OE1  OE2

HACCEPT  HIS01       ND1
HACCEPT  HIS02       NE2 

#add acceptor atom for neutral TYR 
HACCEPT  TYR01       OH
HACCEPT  TYR-1       OH

HACCEPT  PAA01       O1A  O2A
HACCEPT  PAA02       O1A  O2A
HACCEPT  PAA-1       O1A  O2A

HACCEPT  PDD01       O1D  O2D
HACCEPT  PDD02       O1D  O2D
HACCEPT  PDD-1       O1D  O2D

#for water
HDONOR   HOH01      1H   2H		# water molecule
HDONOR   HOH-1       H               	# water molecule
HDONOR   HOH+1      1H   2H   3H        # water molecule

HACCEPT  HOH01       O			# water molecule
HACCEPT  HOH-1       O			# water molecule
HACCEPT  HOH+1       O			# water molecule

# for HE3
#HDONOR   HE3+O       HD1 
HACCEPT  HE3+O       O    #OMA
HDONOR   HE30W      1H   2H
HACCEPT  HE30W       O    #OMA
HDONOR   HE3+W      1H   2H
HACCEPT  HE3+W       O    #OMA
HDONOR   HE30H       H  
HACCEPT  HE30H       O    #OMA
HDONOR   HE3+H       H
HACCEPT  HE3+H       O    #OMA


# for CUB
HDONOR   CUB1W      1H   2H   #cHD1
HACCEPT  CUB1W       O
HDONOR   CUB2W      1H   2H   #cHD1
HACCEPT  CUB2W       O
HDONOR   CUB2I      1H   2H
HACCEPT  CUB2I       O   #cND1
HDONOR   CUB1H       H   #cHD1
HACCEPT  CUB1H       O    
HDONOR   CUB2H       H   #cHD1
HACCEPT  CUB2H       O

# for TYF
HDONOR   TYF01       HH
HACCEPT  TYF01       OH
HACCEPT  TYF-1       OH
HDONOR   TYF+1       HH
HACCEPT  TYF+1       OH
HACCEPT  TYF02       OH

# for HLI
HDONOR   HLI01       HD1
HACCEPT  HLI-1       ND1

# for FA3
HDONOR   FA301       HO1
HACCEPT  FA3-1       O11
HACCEPT  FA301       O11

# for HEA
HACCEPT  HEA01       OMA
HACCEPT  HEA+1       OMA

# for FAL
HDONOR   FAL01       HO1
HACCEPT  FAL01       O11
HACCEPT  FAL-1       O11

