#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mcce.h"

float torsion_conf(CONF *conf_p) {
    float e = 0, phi;
    int i_atom, i_term;
    ATOM *atom0_p, *atom1_p, *atom2_p, *atom3_p;
    TORS tors;
    
    for (i_atom=0; i_atom<conf_p->n_atom; i_atom++) {
        float k_tor_e = 0.;
        if (!conf_p->atom[i_atom].on) continue;
        if (torsion_atoms(conf_p, i_atom, &atom0_p, &atom1_p, &atom2_p, &atom3_p, &tors, 1) == -1) 
            continue;
        phi = torsion_angle(atom0_p->xyz,atom1_p->xyz,atom2_p->xyz,atom3_p->xyz);
        for (i_term=0;i_term<tors.n_term; i_term++) {
            k_tor_e += torsion(phi, tors.V2[i_term], tors.n_fold[i_term], tors.gamma[i_term]);
        }
        e += k_tor_e;

        /* print large single bond torsion 
        if (k_tor_e > 0.5)
            printf("Large torsion: %s-%s-%s-%s, %8.3f\n",atom0_p->name,atom1_p->name,atom2_p->name,atom3_p->name,k_tor_e);
        */
    }
    return e;
}

float torsion_conf_print(CONF *conf_p)
{
    float e = 0, phi;
    int i_atom, i_term;
    ATOM *atom0_p, *atom1_p, *atom2_p, *atom3_p;
    TORS tors;
    
    for (i_atom=0; i_atom<conf_p->n_atom; i_atom++) {
        float k_tor_e = 0.;
        if (!conf_p->atom[i_atom].on) continue;
        if (torsion_atoms(conf_p, i_atom, &atom0_p, &atom1_p, &atom2_p, &atom3_p, &tors, 1) == -1) 
            continue;
        phi = torsion_angle(atom0_p->xyz,atom1_p->xyz,atom2_p->xyz,atom3_p->xyz);
        for (i_term=0;i_term<tors.n_term; i_term++) {
            float e_term = torsion(phi, tors.V2[i_term], tors.n_fold[i_term], tors.gamma[i_term]);
            k_tor_e += e_term;
            
            /* print single bond torsion */
            printf("Torsion: %s %s-%s-%s-%s, angle=%7.1f, torsion term# %2d, e=%8.3f\n",conf_p->uniqID,atom0_p->name,atom1_p->name,atom2_p->name,atom3_p->name,phi/env.d2r,i_term,e_term);
        }
        e += k_tor_e;
    }
    return e;
}

float torsion(float phi, float V2, float n_fold, float gamma) {
    float e;
    e = V2 * (1. + cos(n_fold*phi - gamma));
    return e;
}

VECTOR torsion_torq(float phi, float V2, float n_fold, float gamma, VECTOR k) {
    VECTOR torq;
    torq = vector_rescale(k, n_fold * V2 * sin(n_fold*phi - gamma));
    return torq;
}

int torsion_atoms(CONF *conf_p, int i_atom, ATOM **atom0_p, ATOM **atom1_p, ATOM **atom2_p, ATOM **atom3_p, TORS *tors, int handle) {
    int i_connect;
    char resName[4];
    
    *atom0_p = &(conf_p->atom[i_atom]);
    if (!(*atom0_p)->on) return -1;
    *atom1_p = NULL;
    *atom2_p = NULL;
    *atom3_p = NULL;
    
    /* handle = 1 get atoms according to TORSION parameter, handle=0 just follow the connectivity */
    if (handle == 1) {
        if ( param_get("TORSION",conf_p->confName, (*atom0_p)->name, tors) ) {
            strncpy(resName, conf_p->confName, 3); resName[3] = '\0';
            if ( param_get("TORSION",resName, (*atom0_p)->name, tors) ) {
                return -1;
            }
        }
        for (i_connect=0; i_connect<MAX_CONNECTED; i_connect++) {
            if (!(*atom0_p)->connect12[i_connect]) break;
            if (!(*atom0_p)->connect12[i_connect]->on) continue;
            if (!strcmp((*atom0_p)->connect12[i_connect]->name, tors->atom1)) {
                (*atom1_p) = (*atom0_p)->connect12[i_connect];
                break;
            }
        }
        if (!(*atom1_p)) return -1;
        
        for (i_connect=0; i_connect<MAX_CONNECTED; i_connect++) {
            if (!(*atom1_p)->connect12[i_connect]) break;
            if (!(*atom1_p)->connect12[i_connect]->on) continue;
            if (!strcmp((*atom1_p)->connect12[i_connect]->name, tors->atom2)) {
                (*atom2_p) = (*atom1_p)->connect12[i_connect];
                break;
            }
        }
        if (!(*atom2_p)) return -1;
        
        for (i_connect=0; i_connect<MAX_CONNECTED; i_connect++) {
            if (!(*atom2_p)->connect12[i_connect]) break;
            if (!(*atom2_p)->connect12[i_connect]->on) continue;
            if (!strcmp((*atom2_p)->connect12[i_connect]->name, tors->atom3)) {
                (*atom3_p) = (*atom2_p)->connect12[i_connect];
                break;
            }
        }
        if (!(*atom3_p)) return -1;
    }
    else if (handle == 0) {
        for (i_connect=0; i_connect<MAX_CONNECTED; i_connect++) {
            if (!(*atom0_p)->connect12[i_connect]) break;
            if (!(*atom0_p)->connect12[i_connect]->on) continue;
            if ( (*atom0_p)->connect12[i_connect]->name[1] == 'H') continue;
            (*atom1_p) = (*atom0_p)->connect12[i_connect];
            break;
        }
        if (!(*atom1_p)) return -1;
        
        for (i_connect=0; i_connect<MAX_CONNECTED; i_connect++) {
            if (!(*atom1_p)->connect12[i_connect]) break;
            if (!(*atom1_p)->connect12[i_connect]->on) continue;
            if ((*atom1_p)->connect12[i_connect] == (*atom0_p)) continue;
            if ( (*atom1_p)->connect12[i_connect]->name[1] == 'H') continue;
            (*atom2_p) = (*atom1_p)->connect12[i_connect];
            break;
        }
        if (!(*atom2_p)) return -1;
        
        for (i_connect=0; i_connect<MAX_CONNECTED; i_connect++) {
            if (!(*atom2_p)->connect12[i_connect]) break;
            if (!(*atom2_p)->connect12[i_connect]->on) continue;
            if ( (*atom2_p)->connect12[i_connect]->name[1] == 'H') continue;
            if ((*atom2_p)->connect12[i_connect] == (*atom0_p)) continue;
            if ((*atom2_p)->connect12[i_connect] == (*atom1_p)) continue;
            (*atom3_p) = (*atom2_p)->connect12[i_connect];
            break;
        }
        if (!(*atom3_p)) return -1;
    }
    return 0;
}

float torsion_angle(VECTOR v0, VECTOR v1, VECTOR v2, VECTOR v3) {
    float angle,cos_theta;
    VECTOR i,j,k;
    VECTOR r21,r10,r10_p,r23;
    
    r21 = vector_vminusv(v1,v2);
    r10 = vector_vminusv(v0,v1);
    r23 = vector_vminusv(v3,v2);
    k   = vector_normalize(r21);
    i   = vector_normalize(vector_vminusv(r23, vector_rescale(k, vdotv(k,r23))));
    j   = vector_vxv(k,i);
    r10_p = vector_normalize(vector_vminusv(r10, vector_rescale(k, vdotv(k,r10))));
    
    cos_theta = vdotv(r10_p,i);
    //printf("%8.3f,%8.3f,%8.3f,%8.3f,%8.3f,%8.3f,%8.3f\n",cos_theta,i.x,i.y,i.z,k.x,k.y,k.z);
    if (cos_theta > 1.) cos_theta = 1.;
    else if (cos_theta < -1.) cos_theta = -1.;
    angle = acos(cos_theta);
    
    if (vdotv(r10_p,j) < 0.)
        angle = 2.*env.PI - angle;
    
    return angle;
}


