#!/usr/bin/env python

from data import *

PW_PRINT_CUT = 0.001

def analyze_state_energy(state, ph=7.0, eh=0.0, T=ROOMT, cutoff = PW_PRINT_CUT):
    print("Environment: pH = %.2f  eh = %.f  Temperature = %.2f K" % (ph, eh, T))
    print("Microstate: %s" % (",".join(["%d" % x for x in state])))
    print("Self energy in kCal/mol:")
    print("iConf  CONFORMER    FL  occ eh0-Em0 pKa-pKa0    vdw0    vdw1    epol    tors   dsolv   extra entropy "
          "E_self*occ")

    moving_occ = [x.occ for x in head3lst]
    # Self energy
    E_self = []
    for ic in state:
        conf = head3lst[ic]
        E_ph = T / ROOMT * conf.nh * (ph - conf.pk0) * PH2KCAL
        E_eh = T / ROOMT * conf.ne * (eh - conf.em0) * PH2KCAL / 58.0
        if conf.flag.upper() == "T":
            occ = moving_occ[ic] = conf.occ
        else:
            occ = moving_occ[ic] = 1.0

        E_self_ic = occ * (E_eh + E_ph + conf.vdw0 + conf.vdw1 + conf.epol + conf.tors + conf.dsolv +
                                   conf.extra +
                           conf.entropy)
        E_self.append(E_self_ic)
        print("%05d %s %c %4.2f %7.3f  %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f   %8.3f" % (ic + 1,
                                                                                                   conf.confname,
                                                                                                   conf.flag, occ, E_eh,
                                                                                                   E_ph,
                                                                                                   conf.vdw0, conf.vdw1,
                                                                                                   conf.epol, conf.tors,
                                                                                                   conf.dsolv,
                                                                                                   conf.extra,
                                                                                                   conf.entropy,
                                                                                                   E_self_ic))
    E_self_total = sum(E_self)
    print("Total E_self   %96.3f" % (E_self_total))
    print("\n")

    # Pairwise energy
    E_pw = 0.0
    for ic in state:
        for jc in state:
            E_pw_icjc = pairwise[ic][jc] * moving_occ[ic] * moving_occ[jc]
            E_pw += E_pw_icjc * 0.5    # This is because the interaction will be counted twice A <- B, B <- A
            if abs(E_pw_icjc) > cutoff:
                print("%s <- %s: %5.2f" % (head3lst[ic].confname, head3lst[jc].confname, E_pw_icjc))

    print("Total pairwise interaction: %6.2f" % E_pw)

    print("State energy = E_self + E_pairwise = %.2f + %.2f = %.2f" %(E_self_total, E_pw, E_self_total+E_pw))

    print("\n")
    return
