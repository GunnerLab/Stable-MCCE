#include <stdio.h>
#include <string.h>
#include "mcce.h"
#include <sys/timeb.h>

//int     n_elem;
//float   C6_matrix[N_ELEM_MAX][N_ELEM_MAX];
//float   C12_matrix[N_ELEM_MAX][N_ELEM_MAX];
FILE *vdwf;

void get_vdw0(PROT prot)
{
    int i_res, i_conf, i_atom;
    
    get_connect12(prot);
    setup_vdw_fast(prot);
    
    /* set vdw rad, temporary solution, need to modified surf() directly later -Yifan */
    for (i_res=0; i_res<prot.n_res; i_res++) {
        for (i_conf=1; i_conf<prot.res[i_res].n_conf; i_conf++) {
            for (i_atom=0; i_atom<prot.res[i_res].conf[i_conf].n_atom; i_atom++) {
                if (!prot.res[i_res].conf[i_conf].atom[i_atom].on) continue;
                prot.res[i_res].conf[i_conf].atom[i_atom].rad_backup = prot.res[i_res].conf[i_conf].atom[i_atom].rad;
                prot.res[i_res].conf[i_conf].atom[i_atom].rad = prot.res[i_res].conf[i_conf].atom[i_atom].vdw_rad;
            }
        }
    }
    
    if ((vdwf = fopen("vdw0.lst", "r")) != NULL) {
        fclose(vdwf);
        remove("vdw0.lst");
    }
    if ((vdwf = fopen("vdw0.lst", "a")) == NULL) {
        printf("Can't open vdw0.lst to write\n");
    }
    else {
        fprintf(vdwf, "   confID        vdw0    sas     sum\n");
    }

    for (i_res=0; i_res<prot.n_res; i_res++) {
        get_vdw0_res(i_res,prot);
        //setup_connect_res(prot, i);
        //prot.res[i].conf[0].E_vdw0 = 0.0; /* conformer 0 is defined to have 0 torsion */
        
        //for (j=1; j<prot.res[i].n_conf; j++) {
        //    prot.res[i].conf[j].E_vdw0 = vdw_conf_fast(i, j, i, j, prot, 0);
        //}
        //free_connect_res(prot, i);
    }

    /* resset vdw rad */
    for (i_res=0; i_res<prot.n_res; i_res++) {
        for (i_conf=1; i_conf<prot.res[i_res].n_conf; i_conf++) {
            for (i_atom=0; i_atom<prot.res[i_res].conf[i_conf].n_atom; i_atom++) {
                if (!prot.res[i_res].conf[i_conf].atom[i_atom].on) continue;
                prot.res[i_res].conf[i_conf].atom[i_atom].rad = prot.res[i_res].conf[i_conf].atom[i_atom].rad_backup;
            }
        }
    }
    return;
}

void get_vdw0_no_sas(PROT prot)
{
    int i, j;
    
    get_connect12(prot);
    setup_vdw_fast(prot);
    
    for (i=0; i<prot.n_res; i++) {
        setup_connect_res(prot, i);
        prot.res[i].conf[0].E_vdw0 = 0.0; /* conformer 0 is defined to have 0 torsion */
        
	#pragma omp parallel for shared(i) private(j)
        for (j=1; j<prot.res[i].n_conf; j++) {
            prot.res[i].conf[j].E_vdw0 = vdw_conf_fast(i, j, i, j, prot, 0);
        }
        free_connect_res(prot, i);
    }
    return;
}

void get_vdw0_res(int i_res, PROT prot)
{
    int i_conf, i_atom;
    int k_res, k_conf;
    float raw_vdw=0.0, sas_vdw=0.0, sum_vdw=0.0;

    setup_vdw_fast_res(i_res, prot);
    setup_connect_res(prot, i_res);
    prot.res[i_res].conf[0].E_vdw0 = 0.0; /* conformer 0 is backbone */
    
    for (i_conf=1; i_conf<prot.res[i_res].n_conf; i_conf++) {
        prot.res[i_res].conf[i_conf].E_vdw0 = vdw_conf_fast(i_res, i_conf, i_res, i_conf, prot, 0);
    }
    free_connect_res(prot, i_res);
    
    /* calculate surface term: */
    /* set vdw rad for asp and glu */
    if (!strcmp(prot.res[i_res].resName, "ASP") || !strcmp(prot.res[i_res].resName, "GLU") || !strcmp(prot.res[i_res].resName, "ASN") || !strcmp(prot.res[i_res].resName, "GLN")) {
        for (i_conf=1; i_conf<prot.res[i_res].n_conf; i_conf++) {
            for (i_atom=0; i_atom<prot.res[i_res].conf[i_conf].n_atom; i_atom++) {
                if (!prot.res[i_res].conf[i_conf].atom[i_atom].on) continue;
                if (prot.res[i_res].conf[i_conf].atom[i_atom].name[1] == 'O') {
                    prot.res[i_res].conf[i_conf].atom[i_atom].vdw_rad *= 2.0;
                }
                if (prot.res[i_res].conf[i_conf].atom[i_atom].name[1] == 'N') {
                    prot.res[i_res].conf[i_conf].atom[i_atom].vdw_rad *= 2.0;
                }
            }
        }
    }
    
    /* set conf.on flag */
    for (k_res=0; k_res<prot.n_res; k_res++) {
        prot.res[k_res].conf[0].on = 1;
        if (prot.res[k_res].n_conf > 1) prot.res[k_res].conf[1].on = 1;
        for (k_conf=2; k_conf<prot.res[k_res].n_conf; k_conf++) prot.res[k_res].conf[k_conf].on = 0;
    }
    
    surfw_res(prot, i_res, env.radius_probe);
    for (i_conf=1; i_conf<prot.res[i_res].n_conf; i_conf++) {
        raw_vdw = prot.res[i_res].conf[i_conf].E_vdw0;
        if (!strcmp(prot.res[i_res].resName, "ASP") || !strcmp(prot.res[i_res].resName, "GLU") || !strcmp(prot.res[i_res].resName, "ASN") || !strcmp(prot.res[i_res].resName, "GLN")) {
            prot.res[i_res].conf[i_conf].E_vdw0 += env.sas2vdw * prot.res[i_res].conf[i_conf].sas * 0.25;
        }
        else {
            prot.res[i_res].conf[i_conf].E_vdw0 += env.sas2vdw * prot.res[i_res].conf[i_conf].sas;
        }
        sum_vdw = prot.res[i_res].conf[i_conf].E_vdw0;
        sas_vdw = sum_vdw - raw_vdw;
        fprintf(vdwf, "%s%8.3f%8.3f%8.3f\n", prot.res[i_res].conf[i_conf].uniqID, raw_vdw, sas_vdw, sum_vdw);        
    }
    
    /* reset conf.on flag */
    for (k_res=0; k_res<prot.n_res; k_res++) {
        for (k_conf=0; k_conf<prot.res[k_res].n_conf; k_conf++)
            prot.res[k_res].conf[k_conf].on = 1;
    }
    
    /* reset vdw rad for asp and glu */
    if (!strcmp(prot.res[i_res].resName, "ASP") || !strcmp(prot.res[i_res].resName, "GLU") || !strcmp(prot.res[i_res].resName, "ASN") || !strcmp(prot.res[i_res].resName, "GLN")) {
        for (i_conf=1; i_conf<prot.res[i_res].n_conf; i_conf++) {
            for (i_atom=0; i_atom<prot.res[i_res].conf[i_conf].n_atom; i_atom++) {
                if (!prot.res[i_res].conf[i_conf].atom[i_atom].on) continue;
                if (prot.res[i_res].conf[i_conf].atom[i_atom].name[1] == 'O') {
                    prot.res[i_res].conf[i_conf].atom[i_atom].vdw_rad /= 2.0;
                }
                if (prot.res[i_res].conf[i_conf].atom[i_atom].name[1] == 'N') {
                    prot.res[i_res].conf[i_conf].atom[i_atom].vdw_rad /= 2.0;
                }
            }
        }
    }
}

void get_vdw0_res_no_sas(int i_res, PROT prot)
{
    int i_conf;
    setup_vdw_fast_res(i_res, prot);
    setup_connect_res(prot, i_res);
    prot.res[i_res].conf[0].E_vdw0 = 0.0; /* conformer 0 is backbone */
    
    for (i_conf=1; i_conf<prot.res[i_res].n_conf; i_conf++) {
        prot.res[i_res].conf[i_conf].E_vdw0 = vdw_conf_fast(i_res, i_conf, i_res, i_conf, prot, 0);
    }
    free_connect_res(prot, i_res);
}

