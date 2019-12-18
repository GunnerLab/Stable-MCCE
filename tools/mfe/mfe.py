#!/usr/bin/python
import sys

class CONFORMER:
   def __init__(self):
      self.occ = []

class E_IONIZE:
   def __init__(self):
      self.mfe = []

#Global variables and default values
ph1         = 0.0
eh1         = 0.0
scale_ele   = 1.0
scale_vdw   = 1.0
scale_vdw0  = 1.0
scale_vdw1  = 1.0
scale_tors  = 1.0
scale_dsolv = 1.0
fname_extra = ""
ph2Kcal = 1.364
mev2Kcal = 0.0235
Kcal2kT  = 1.688

conformers = []
residues = []
titration_range = []
titration_type = ""


def first_ph():
   global ph1, eh1, scale_ele, scale_vdw, scale_vdw0, scale_vdw1, scale_tors, scale_dsolv,fname_extra
   lines=open("run.prm").readlines()
   for line in lines:
      if    line.find("(TITR_PH0)")>=0:
         ph1 = float(line.split()[0])
      elif line.find("(TITR_EH0)")>=0:
         eh1 = float(line.split()[0])
      elif line.find("(EXTRA)")>=0:
         fname_extra = line.split()[0]

   
   if fname_extra:
      lines=open(fname_extra).readlines()
      for line in lines:
         fields = line.split()
         if len(fields) != 3:
	    continue  
	 elif fields[0]+fields[1] == "SCALINGVDW0":
	    scale_vdw0 = float(fields[2])
	 elif fields[0]+fields[1] == "SCALINGVDW1":
	    scale_vdw1 = float(fields[2])
	 elif fields[0]+fields[1] == "SCALINGVDW":
	    scale_vdw = float(fields[2])
	 elif fields[0]+fields[1] == "SCALINGTORS":
	    scale_tors = float(fields[2])
	 elif fields[0]+fields[1] == "SCALINGELE":
	    scale_ele = float(fields[2])
	 elif fields[0]+fields[1] == "SCALINGDSOLV":
	    scale_dsolv = float(fields[2])
	     
   for line in lines:
      if   line.find("(SCALE_ELE)")>=0:
         scale_ele = float(line.split()[0])
      elif line.find("(SCALE_VDW)")>=0:
         scale_vdw = float(line.split()[0])

   return

def read_headlst():
   global conformers
   global titration_range,titration_type,titration_unit, scale_ele, scale_vdw

   lines = open("head3.lst").readlines()
   lines.pop(0)   # remove the title line

   if len(conformers) > 0:
      print "WARNING: adding to non empty conformer list."

   for line in lines:
      fields = line.split()
      conformer = CONFORMER()
      conformer.id   = fields[1]
      conformer.crg  = float(fields[4])
      conformer.Em0  = float(fields[5])
      conformer.pKa0 = float(fields[6])
      conformer.ne   = int(fields[7])
      conformer.nH   = int(fields[8])
      conformer.vdw0 = float(fields[9]) * scale_vdw0
      conformer.vdw1 = float(fields[10])* scale_vdw1
      conformer.tors = float(fields[11])* scale_tors
      conformer.epol = float(fields[12])* scale_ele
      conformer.dsolv= float(fields[13])* scale_dsolv
      conformer.extra= float(fields[14])
      conformer.self = conformer.vdw0+conformer.vdw1+conformer.tors+conformer.epol+conformer.dsolv+conformer.extra
      conformers.append(conformer)
   return

def read_fort38():
   global conformers
   global titration_range,titration_type,titration_unit

   lines = open("fort.38").readlines()

   temp = lines[0].split()
   ttype = temp[0].strip().upper()
   if   ttype == 'EH':
      titration_type = 'Eh'
      titration_unit = 'mV'
   elif ttype == 'PH':
      titration_type = 'pH'
      titration_unit = 'pH'
   else:
      titration_type = ttype
      titration_unit = '?'
   titration_range = [float(x) for x in temp[1:]]
   lines.pop(0)

   for i in range(len(conformers)):
      line = lines[i].split()
      if conformers[i].id != line[0]:
         print "ERROR, %s in head3.lst doesn't match %s in fort.38" %(conformers[i].id, line[0])
         return
      conformers[i].occ = [float(x) for x in line[1:]]

   return

def group_residues():
   global conformers
   global titration_range,titration_type,titration_unit, residues

   if len(conformers) < 1: return residues    # no conformers
   old_resid = conformers[0].id[:3] + conformers[0].id[5:11]
   residue = [old_resid, [], []]
   if conformers[0].id[3] == '0':
      residue[1].append(conformers[0])
   else:
      residue[2].append(conformers[0])

   for conformer in conformers[1:]:
      current_resid = conformer.id[:3] + conformer.id[5:11]
      if current_resid == old_resid:
         if conformer.id[3] == '0':
            residue[1].append(conformer)
         else:
            residue[2].append(conformer)
      else: # a new residue
         residues.append(residue)
         old_resid = current_resid
         residue = [old_resid, [], []]
         if conformer.id[3] == '0':
            residue[1].append(conformer)
         else:
            residue[2].append(conformer)

   residues.append(residue)

   return

def read_pK():
   lines = open("pK.out").readlines()
   pK = [[line[:10], line[10:]] for line in lines]
   return pK

def E_ionize(res_pKa_name):
   global residues
   import math

   res_name = res_pKa_name[:3] + res_pKa_name[4:]
   found = 0
   for residue in residues:
      #print res_name, residue[0]   
      if res_name == residue[0]:
         found = 1
         break

   if found==0:
      print "Residue %s not found in fort.38" % res_name
      sys.exit(0)

   dG_ionize = E_IONIZE()

   # mfe for each conformer
   for conformer in residue[1]+residue[2]:
      conformer.pHeffect = [0.0 for i in range(len(titration_range))]
      conformer.Eheffect = [0.0 for i in range(len(titration_range))]
      conformer.res_mfe =  [[0.0 for i in range(len(titration_range))] for x in residues]
      conformer.mfe_total= [0.0 for i in range(len(titration_range))]
      conformer.E_total =  [0.0 for i in range(len(titration_range))]

      # pairwise energy table
      lines = open("energies/"+conformer.id+".opp").readlines()
      pairwise = {}
      for line in lines:
         line = line.split()
         if len(line) == 4:
            pairwise[line[1]] = float(line[2])*scale_ele + float(line[3])*scale_vdw

      # make conformer mfe
      for i in range(len(titration_range)):
         point = titration_range[i]
         conf_mfe = [0.0 for x in titration_range]

         # pH effect in Kcal/mol
         if titration_type == 'pH':
            conformer.pHeffect[i] = (point-conformer.pKa0)*conformer.nH*ph2Kcal
         else:
            conformer.pHeffect[i] = (ph1-conformer.pKa0)*conformer.nH*ph2Kcal

         # Eh effect in Kcal/mol
         if titration_type == 'Eh':
            conformer.Eheffect[i] = (point-conformer.Em0)*conformer.ne*mev2Kcal
         else:
            conformer.Eheffect[i] = (eh1-conformer.Em0)*conformer.ne*mev2Kcal

         for j in range(len(residues)):
            res = residues[j]
            if res[0] == residue[0]:
               mfe = 0.0
            else:
               mfe = 0.0
               for conf in res[1]+res[2]:
                  if not pairwise.has_key(conf.id):pairwise[conf.id] = 0.0
                  mfe+=pairwise[conf.id]*conf.occ[i]
            # This mfe is at 1 titration point, from one residue
            conformer.res_mfe[j][i] = mfe
            conformer.mfe_total[i] += mfe

         # update conformer E_total
         conformer.E_total[i] = conformer.mfe_total[i]\
                                +conformer.pHeffect[i]\
                                +conformer.Eheffect[i]\
                                +conformer.self


   Eref=residue[1][0].E_total[0]    # reference E (lowest of this res)
   for conf in residue[1]+residue[2]:
      for i in range(len(titration_range)):
         if Eref > conf.E_total[i]: Eref=conf.E_total[i]

   # Calculate mfe occupancy of each conformer
   SigmaE = [0.0 for i in range(len(titration_range))]
   for i in range(len(titration_range)):
      Ei = [math.exp(-(conformer.E_total[i]-Eref)*Kcal2kT) for conformer in residue[1]+residue[2]]
      for x in Ei: SigmaE[i] += x

   for conformer in residue[1]+residue[2]:
      conformer.rocc = [0.0 for x in titration_range]
      for i in range(len(titration_range)):
         conformer.rocc[i] = math.exp(-(conformer.E_total[i]-Eref)*Kcal2kT)/SigmaE[i]  # recovered occ

   SigmaOcc = [0.0 for i in range(len(titration_range))]
   for conformer in residue[1]:
      for i in range(len(titration_range)):SigmaOcc[i] += conformer.rocc[i]
   for conformer in residue[1]:
      conformer.nocc = [conformer.rocc[i]/SigmaOcc[i] for i in range(len(titration_range))]

   SigmaOcc = [0.0 for i in range(len(titration_range))]
   for conformer in residue[2]:
      for i in range(len(titration_range)):SigmaOcc[i] += conformer.rocc[i]
   for conformer in residue[2]:
      conformer.nocc = [conformer.rocc[i]/SigmaOcc[i] for i in range(len(titration_range))]


   # energy terms of ground state
   ground_state = E_IONIZE()
   ground_state.vdw0  = [0.0 for x in titration_range]
   ground_state.vdw1  = [0.0 for x in titration_range]
   ground_state.tors  = [0.0 for x in titration_range]
   ground_state.epol  = [0.0 for x in titration_range]
   ground_state.dsolv = [0.0 for x in titration_range]
   ground_state.extra = [0.0 for x in titration_range]
   ground_state.pHeffect= [0.0 for x in titration_range]
   ground_state.Eheffect= [0.0 for x in titration_range]
   ground_state.mfe_total=[0.0 for x in titration_range]
   ground_state.res_mfe = [[0.0 for x in titration_range] for x in residues]
   ground_state.E_total  =[0.0 for x in titration_range]
   ground_state.TS       =[0.0 for x in titration_range]
   for conformer in residue[1]:
      for i in range(len(titration_range)):
         ground_state.vdw0[i]     += conformer.nocc[i]*conformer.vdw0
         ground_state.vdw1[i]     += conformer.nocc[i]*conformer.vdw1
         ground_state.tors[i]     += conformer.nocc[i]*conformer.tors
         ground_state.epol[i]     += conformer.nocc[i]*conformer.epol
         ground_state.dsolv[i]    += conformer.nocc[i]*conformer.dsolv
         ground_state.extra[i]   += conformer.nocc[i]*conformer.extra
         ground_state.pHeffect[i] += conformer.nocc[i]*conformer.pHeffect[i]
         ground_state.Eheffect[i] += conformer.nocc[i]*conformer.Eheffect[i]
         ground_state.mfe_total[i]+= conformer.nocc[i]*conformer.mfe_total[i]
         ground_state.E_total[i]  += conformer.nocc[i]*conformer.E_total[i]
         if conformer.nocc[i]>0.000001:
            ground_state.TS[i]       +=-conformer.nocc[i]*math.log(conformer.nocc[i])/Kcal2kT
      for j in range(len(residues)):
         for i in range(len(titration_range)):
            ground_state.res_mfe[j][i] += conformer.nocc[i]*conformer.res_mfe[j][i]

   ground_state.G = [ground_state.E_total[i] - ground_state.TS[i] for i in range(len(titration_range))]

   # energy terms of charged state
   charged_state = E_IONIZE()
   charged_state.vdw0  = [0.0 for x in titration_range]
   charged_state.vdw1  = [0.0 for x in titration_range]
   charged_state.tors  = [0.0 for x in titration_range]
   charged_state.epol  = [0.0 for x in titration_range]
   charged_state.dsolv = [0.0 for x in titration_range]
   charged_state.extra = [0.0 for x in titration_range]
   charged_state.pHeffect= [0.0 for x in titration_range]
   charged_state.Eheffect= [0.0 for x in titration_range]
   charged_state.res_mfe = [[0.0 for x in titration_range] for x in residues]
   charged_state.mfe_total=[0.0 for x in titration_range]
   charged_state.E_total  =[0.0 for x in titration_range]
   charged_state.TS       =[0.0 for x in titration_range]
   for conformer in residue[2]:
      for i in range(len(titration_range)):
         charged_state.vdw0[i]     += conformer.nocc[i]*conformer.vdw0
         charged_state.vdw1[i]     += conformer.nocc[i]*conformer.vdw1
         charged_state.tors[i]     += conformer.nocc[i]*conformer.tors
         charged_state.epol[i]     += conformer.nocc[i]*conformer.epol
         charged_state.dsolv[i]    += conformer.nocc[i]*conformer.dsolv
         charged_state.extra[i]   += conformer.nocc[i]*conformer.extra
         charged_state.pHeffect[i] += conformer.nocc[i]*conformer.pHeffect[i]
         charged_state.Eheffect[i] += conformer.nocc[i]*conformer.Eheffect[i]
         charged_state.mfe_total[i]+= conformer.nocc[i]*conformer.mfe_total[i]
         charged_state.E_total[i]  += conformer.nocc[i]*conformer.E_total[i]
         if conformer.nocc[i]>0.000001:
            charged_state.TS[i]       +=-conformer.nocc[i]*math.log(conformer.nocc[i])/Kcal2kT
      for j in range(len(residues)):
         for i in range(len(titration_range)):
            charged_state.res_mfe[j][i] += conformer.nocc[i]*conformer.res_mfe[j][i]

   charged_state.G = [charged_state.E_total[i] - charged_state.TS[i] for i in range(len(titration_range))]

   dG_ionize.ground_state = ground_state
   dG_ionize.charged_state = charged_state
   dG_ionize.ground_confs = residue[1]
   dG_ionize.charged_confs = residue[2]
   dG_ionize.resID = residue[0]

   return dG_ionize


if __name__ == '__main__':

   if (len(sys.argv) < 3):
      print "mfe.py res_id titration_point [pH_cutoff]"
      print "   res_id:          The residue ID in pK.out"
      print "   titration_point: Titration point that mfe is calculated."
      print "                    If this number is between two calculated titration"
      print "                    values, a linear average will be performed."
      print "   pH_cutoff:       display pairwise interaction bigger than this value"
      sys.exit(0)
   else:
      t_point = float(sys.argv[2])
      
   pH_cutoff = -0.001
   if len(sys.argv) > 3:
      pH_cutoff = float(sys.argv[3])

   # read run.prm
   first_ph()

   #read head list
   read_headlst()

   # read pK.out
   pK = read_pK()

   # read fort.38
   read_fort38()

   group_residues()

   # decide which two columns will be used to get residue mfe
   t_low_found = t_high_found = 0
   for i in range(len(titration_range)):
      if t_point > titration_range[i]-0.001:
         i_low = i
         t_low_found = 1

   for i in range(len(titration_range)):
      if t_point < titration_range[i]+0.001:
         i_high = i
         t_high_found = 1
         break

   if not t_low_found or not t_high_found:
      print "titration_point out off range"
      sys.exit(0)

   dG = E_ionize(sys.argv[1]);
   dG_point = E_IONIZE()

   if abs(titration_range[i_low] -titration_range[i_high]) < 0.01: # one point
      dG_point.vdw0 = dG.charged_state.vdw0[i_low] - dG.ground_state.vdw0[i_low]
      dG_point.vdw1 = dG.charged_state.vdw1[i_low] - dG.ground_state.vdw1[i_low]
      dG_point.tors = dG.charged_state.tors[i_low] - dG.ground_state.tors[i_low]
      dG_point.epol = dG.charged_state.epol[i_low] - dG.ground_state.epol[i_low]
      dG_point.dsolv = dG.charged_state.dsolv[i_low] - dG.ground_state.dsolv[i_low]
      dG_point.extra = dG.charged_state.extra[i_low] - dG.ground_state.extra[i_low]
      dG_point.pHeffect = dG.charged_state.pHeffect[i_low] - dG.ground_state.pHeffect[i_low]
      dG_point.Eheffect = dG.charged_state.Eheffect[i_low] - dG.ground_state.Eheffect[i_low]
      dG_point.mfe_total = dG.charged_state.mfe_total[i_low] - dG.ground_state.mfe_total[i_low]
      dG_point.E_total = dG.charged_state.E_total[i_low] - dG.ground_state.E_total[i_low]
      dG_point.TS = dG.charged_state.TS[i_low] - dG.ground_state.TS[i_low]
      for i in range(len(dG.charged_state.res_mfe)):
         dG_point.mfe.append(dG.charged_state.res_mfe[i][i_low] - dG.ground_state.res_mfe[i][i_low])
      dG_point.G = dG.charged_state.G[i_low] - dG.ground_state.G[i_low]

   else: # scale average of two points
      k = (t_point - titration_range[i_low])/(titration_range[i_high] - titration_range[i_low])
      dG_point.vdw0 = (1-k)*(dG.charged_state.vdw0[i_low] - dG.ground_state.vdw0[i_low])\
                      +   k*(dG.charged_state.vdw0[i_high] - dG.ground_state.vdw0[i_high])
      dG_point.vdw1 = (1-k)*(dG.charged_state.vdw1[i_low] - dG.ground_state.vdw1[i_low])\
                      +   k*(dG.charged_state.vdw1[i_high] - dG.ground_state.vdw1[i_high])
      dG_point.tors = (1-k)*(dG.charged_state.tors[i_low] - dG.ground_state.tors[i_low])\
                      +   k*(dG.charged_state.tors[i_high] - dG.ground_state.tors[i_high])
      dG_point.epol = (1-k)*(dG.charged_state.epol[i_low] - dG.ground_state.epol[i_low])\
                      +   k*(dG.charged_state.epol[i_high] - dG.ground_state.epol[i_high])
      dG_point.dsolv =(1-k)*(dG.charged_state.dsolv[i_low] - dG.ground_state.dsolv[i_low])\
                      +   k*(dG.charged_state.dsolv[i_high] - dG.ground_state.dsolv[i_high])
      dG_point.extra =(1-k)*(dG.charged_state.extra[i_low] - dG.ground_state.extra[i_low])\
                      +   k*(dG.charged_state.extra[i_high] - dG.ground_state.extra[i_high])
      dG_point.pHeffect =(1-k)*(dG.charged_state.pHeffect[i_low] - dG.ground_state.pHeffect[i_low])\
                      +      k*(dG.charged_state.pHeffect[i_high] - dG.ground_state.pHeffect[i_high])
      dG_point.Eheffect =(1-k)*(dG.charged_state.Eheffect[i_low] - dG.ground_state.Eheffect[i_low])\
                      +      k*(dG.charged_state.Eheffect[i_high] - dG.ground_state.Eheffect[i_high])
      dG_point.mfe_total =(1-k)*(dG.charged_state.mfe_total[i_low] - dG.ground_state.mfe_total[i_low])\
                      +      k*(dG.charged_state.mfe_total[i_high] - dG.ground_state.mfe_total[i_high])
      dG_point.E_total =(1-k)*(dG.charged_state.E_total[i_low] - dG.ground_state.E_total[i_low])\
                      +      k*(dG.charged_state.E_total[i_high] - dG.ground_state.E_total[i_high])
      dG_point.TS = (1-k)*(dG.charged_state.TS[i_low] - dG.ground_state.TS[i_low])\
                      + k*(dG.charged_state.TS[i_high] - dG.ground_state.TS[i_high])
      for i in range(len(dG.charged_state.res_mfe)):
         dG_point.mfe.append((1-k)*(dG.charged_state.res_mfe[i][i_low] - dG.ground_state.res_mfe[i][i_low])\
                            +k*(dG.charged_state.res_mfe[i][i_high] - dG.ground_state.res_mfe[i][i_high]))

      dG_point.G = (1-k)*(dG.charged_state.G[i_low] - dG.ground_state.G[i_low])\
                      +k*(dG.charged_state.G[i_high] - dG.ground_state.G[i_high])



   # print
   for x in pK:
      if (sys.argv[1][:3] == x[0][:3] and sys.argv[1][4:] == x[0][4:]):
         res_pKa = x[1].split()[0]
   print "Residue %s pKa/Em=%s" % (sys.argv[1], res_pKa)
   print "================================="
   print "Terms          pH     meV    Kcal"
   print "---------------------------------"
   print "vdw0     %8.2f%8.2f%8.2f" % (dG_point.vdw0/ph2Kcal, dG_point.vdw0/mev2Kcal, dG_point.vdw0)
   print "vdw1     %8.2f%8.2f%8.2f" % (dG_point.vdw1/ph2Kcal, dG_point.vdw1/mev2Kcal, dG_point.vdw1)
   print "tors     %8.2f%8.2f%8.2f" % (dG_point.tors/ph2Kcal, dG_point.tors/mev2Kcal, dG_point.tors)
   print "ebkb     %8.2f%8.2f%8.2f" % (dG_point.epol/ph2Kcal, dG_point.epol/mev2Kcal, dG_point.epol)
   print "dsol     %8.2f%8.2f%8.2f" % (dG_point.dsolv/ph2Kcal, dG_point.dsolv/mev2Kcal, dG_point.dsolv)
   print "offset   %8.2f%8.2f%8.2f" % (dG_point.extra/ph2Kcal, dG_point.extra/mev2Kcal, dG_point.extra)
   print "pH&pK0   %8.2f%8.2f%8.2f" % (dG_point.pHeffect/ph2Kcal, dG_point.pHeffect/mev2Kcal, dG_point.pHeffect)
   print "Eh&Em0   %8.2f%8.2f%8.2f" % (dG_point.Eheffect/ph2Kcal, dG_point.Eheffect/mev2Kcal, dG_point.Eheffect)
   print "-TS      %8.2f%8.2f%8.2f" % (-dG_point.TS/ph2Kcal, -dG_point.TS/mev2Kcal, -dG_point.TS)
   print "residues %8.2f%8.2f%8.2f" % (dG_point.mfe_total/ph2Kcal, dG_point.mfe_total/mev2Kcal, dG_point.mfe_total)
   print "*********************************"
   print "TOTAL    %8.2f%8.2f%8.2f" % (dG_point.G/ph2Kcal, dG_point.G/mev2Kcal, dG_point.G)
   print "*********************************"
   for i in range(len(dG_point.mfe)):
      if abs(dG_point.mfe[i]/ph2Kcal) > pH_cutoff:
         print "%-9s%8.2f%8.2f%8.2f" % (residues[i][0], dG_point.mfe[i]/ph2Kcal, dG_point.mfe[i]/mev2Kcal, dG_point.mfe[i])
   print "================================="
