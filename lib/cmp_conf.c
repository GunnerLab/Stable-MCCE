#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mcce.h"

int cmp_conf(CONF conf1, CONF conf2, float IDEN_THR) {
    int iatom, jatom,n_on1,n_on2;
    float IDEN_THR2 = IDEN_THR*IDEN_THR;
    
    //if (strcmp(conf1.confName, conf2.confName)) return -1;
    if (conf1.n_atom != conf2.n_atom) return -1;
    
    n_on1 = 0;
    for (iatom = 0; iatom<conf1.n_atom; iatom++) {
        if (conf1.atom[iatom].on) n_on1++;
    }
    n_on2 = 0;
    for (jatom = 0; jatom<conf2.n_atom; jatom++) {
        if (conf2.atom[jatom].on) n_on2++;
    }
    if (n_on1 != n_on2) return -2;
    
    for (iatom = 0; iatom<conf1.n_atom; iatom++) {
        if (!conf1.atom[iatom].on) continue;
        for (jatom = 0; jatom<conf2.n_atom; jatom++) {
            if (!conf2.atom[jatom].on) continue;
            
            if (strncmp(conf1.atom[iatom].name+1, conf2.atom[jatom].name+1, 2) ) continue;
            if (fabs(conf1.atom[iatom].crg - conf2.atom[jatom].crg) > 1e-3) continue;
            if (ddvv(conf1.atom[iatom].xyz, conf2.atom[jatom].xyz) < IDEN_THR2) break;
        }
        if (jatom == conf2.n_atom) return -1;
    }
    return 0;
}

int cmp_conf_hv(CONF conf1, CONF conf2, float IDEN_THR) {
    int iatom, jatom;
    float IDEN_THR2 = IDEN_THR*IDEN_THR;
    int n_hv1=0, n_hv2=0;
    //if (strcmp(conf1.confName, conf2.confName)) return -1;
    //if (conf1.n_atom != conf2.n_atom) return -1;
    /* check if they same number of heavy atoms */
    for (iatom = 0; iatom<conf1.n_atom; iatom++) {
        if (!conf1.atom[iatom].on) continue;
        if (conf1.atom[iatom].name[1] == 'H') continue;
        n_hv1++;
    }
    for (jatom = 0; jatom<conf2.n_atom; jatom++) {
        if (!conf2.atom[jatom].on) continue;
        if (conf2.atom[jatom].name[1] == 'H') continue;
        n_hv2++;
    }
    if (n_hv1!=n_hv2) return -1;
    
    for (iatom = 0; iatom<conf1.n_atom; iatom++) {
        if (!conf1.atom[iatom].on) continue;
        if (conf1.atom[iatom].name[1] == 'H') continue;
        
        for (jatom = 0; jatom<conf2.n_atom; jatom++) {
            if (!conf2.atom[jatom].on) continue;
            if (conf2.atom[jatom].name[1] == 'H') continue;
            
            if (strncmp(conf1.atom[iatom].name+1, conf2.atom[jatom].name+1, 2) ) continue;
            if (ddvv(conf1.atom[iatom].xyz, conf2.atom[jatom].xyz) < IDEN_THR2) break;
        }
        if (jatom == conf2.n_atom) return -1;
    }
    return 0;
}

float dist_conf_hv(CONF conf1, CONF conf2) {
    int iatom, jatom;
    float max_dist = 999.;
    float next_test_dist, test_dist2;
    int n_hv1=0, n_hv2=0;

    /* check if they have the same number of heavy atoms */
    for (iatom = 0; iatom<conf1.n_atom; iatom++) {
        if (!conf1.atom[iatom].on) continue;
        if (conf1.atom[iatom].name[1] == 'H') continue;
        n_hv1++;
    }
    for (jatom = 0; jatom<conf2.n_atom; jatom++) {
        if (!conf2.atom[jatom].on) continue;
        if (conf2.atom[jatom].name[1] == 'H') continue;
        n_hv2++;
    }
    if (n_hv1!=n_hv2) return 999.;
    
    next_test_dist = max_dist;
    while (1) {
        test_dist2 = next_test_dist*next_test_dist - 1e-4;
        next_test_dist = 0.;
        for (iatom = 0; iatom<conf1.n_atom; iatom++) {
            if (!conf1.atom[iatom].on) continue;
            if (conf1.atom[iatom].name[1] == 'H') continue;
            
            for (jatom = 0; jatom<conf2.n_atom; jatom++) {
                if (!conf2.atom[jatom].on) continue;
                if (conf2.atom[jatom].name[1] == 'H') continue;
                
                if (strncmp(conf1.atom[iatom].name+1, conf2.atom[jatom].name+1, 2) ) continue;
                if (ddvv(conf1.atom[iatom].xyz, conf2.atom[jatom].xyz) < test_dist2) break;
            }
            
            if (jatom >= conf2.n_atom) break; /* can't find a jatom within test_dist to match iatom */
            else {
                /* next test threshold is the maximum of all distances btw matched iatom and jatom */
                if (dvv(conf1.atom[iatom].xyz, conf2.atom[jatom].xyz) > next_test_dist) {
                    next_test_dist = dvv(conf1.atom[iatom].xyz, conf2.atom[jatom].xyz);
                }
            }
        }
        if (iatom < conf1.n_atom) { /* test_dist in the current loop is shorter than the distance btw i and j */
            break;
        }
    }
    return sqrt(test_dist2+1e-4);
}

float rmsd_conf_hv(CONF conf1, CONF conf2) {
    int iatom, jatom;
    int n_hv1=0, n_hv2=0;
    float sum_distsq;
    
    /* check if they have the same number of heavy atoms */
    for (iatom = 0; iatom<conf1.n_atom; iatom++) {
        if (!conf1.atom[iatom].on) continue;
        if (conf1.atom[iatom].name[1] == 'H') continue;
        n_hv1++;
    }
    for (jatom = 0; jatom<conf2.n_atom; jatom++) {
        if (!conf2.atom[jatom].on) continue;
        if (conf2.atom[jatom].name[1] == 'H') continue;
        n_hv2++;
    }
    if (n_hv1!=n_hv2) return 999.;
    if (n_hv1 == 0) return 0.;
    
    sum_distsq = 0.;
    for (iatom = 0; iatom<conf1.n_atom; iatom++) {
        if (!conf1.atom[iatom].on) continue;
        if (conf1.atom[iatom].name[1] == 'H') continue;
        
        for (jatom = 0; jatom<conf2.n_atom; jatom++) {
            if (!conf2.atom[jatom].on) continue;
            if (conf2.atom[jatom].name[1] == 'H') continue;
            
            if (!strcmp(conf1.atom[iatom].name, conf2.atom[jatom].name)) break;
        }
        if (jatom >= conf2.n_atom) return 999.; /* can't find a jatom to match iatom */
        
        sum_distsq += ddvv(conf1.atom[iatom].xyz, conf2.atom[jatom].xyz);
    }
    
    return sqrt(sum_distsq/(float) n_hv1);
}

