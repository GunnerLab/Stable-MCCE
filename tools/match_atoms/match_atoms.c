#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mcce.h"
#define  MATCH_THR 1.
#define  BOND_THR 2

extern FILE *premcce_confname(FILE *fp);
extern int create_param(FILE *pdb_fp, int k_line);
extern int place_missing(PROT prot, int handle_addconf);
GEOM geom_fit(ATOM *atoms1, ATOM *atoms2, int na);

typedef struct {
    char name[3];
    int  na;
    int  *ia;
    int  iswitch;
    int  jswitch;
} ELEM;

ENV env;

int main(int argc, char *argv[]) {
    
    FILE *fp1,*fp2, *tmp_fp, *tmp_fp2;
    int  na1=0, na2=0;
    int  ia1, ia2, ja1;
    int  nelem1=0, nelem2=0;
    int  ielem1, ielem2;
    int  jelem1, jelem2;
    int  kelem1, kelem2;
    int  na_ielem, na_jelem;
    int  ia_elem1;
    int  nelem_match, natom_match;
    int  katom_match;
    ELEM *elems1=NULL, *elems2=NULL, elem_buffer;
    PROT prot;
    ATOM *atoms1=NULL, *atoms2=NULL;
    ATOM *atoms_to_match1=NULL, *atoms_to_match2=NULL;
    char  line[MAXCHAR_LINE];   /* line buffer */
    char  elem_name[3];
    
    if (argc<3) {
        printf("   Usage: match_atoms file1.pdb file2.pdb [param_dir]\n");
        exit(-1);
    }
    
    fp1 = fopen(argv[1],"r");
    while ( memset(line,0,MAXCHAR_LINE*sizeof(char)),
    fgets(line, MAXCHAR_LINE, fp1)) {
        if (strncmp("ATOM  ", line, 6) && strncmp("HETATM", line, 6)) continue;
        if (line[strlen(line) - 1] == '\n') line[strlen(line) - 1] = '\0'; /* this is to make sure conformer history not end with a new line */
        
        na1++;
        atoms1 = realloc(atoms1,na1*sizeof(ATOM));
        atoms1[na1-1] = pdbline2atom(line);
    }
    
    if (argc>=4) {
        db_open();
        if (load_all_param(argv[3])) {printf("   FATAL: \"failed.\"\n"); db_close(); return USERERR;}
        
        fp2 = fopen(argv[2],"r");
        tmp_fp = tmpfile();
        while ( memset(line,0,MAXCHAR_LINE*sizeof(char)),
        fgets(line, MAXCHAR_LINE, fp2)) {
            int i;
            if (strncmp("ATOM  ", line, 6) && strncmp("HETATM", line, 6)) continue;
            for (i=17;i<19;i++) if (line[i] == ' ') line[i]='_';
            fprintf(tmp_fp,"%s\n",line);
        }
        fclose(fp2);
        
        strcpy(env.new_tpl,"new.tpl");
        strcpy(env.debug_log,"debug.log");
        while (!(tmp_fp2 = premcce_confname(tmp_fp))) {
            if (env.err_msg >= 0) {
                printf("   Creating temporary parameter file for unrecognized residue...\n");
                if (create_param(tmp_fp, env.err_msg)) printf("   parameter file %s not made.\n",env.new_tpl),exit(-1);
                load_param(env.new_tpl);
                rewind(tmp_fp);
                printf("   Trying labeling again...\n");
            }
            else {
                printf("   STOP: fatal errors reported by premcce_confname()\n");
                return USERERR;
            }
        }
        fclose(tmp_fp);
        prot = load_pdb(tmp_fp2);
        place_missing(prot,0);
        tmp_fp = tmpfile();
        write_pdb(tmp_fp,prot);
        rewind(tmp_fp);
        db_close();
        
        while ( memset(line,0,MAXCHAR_LINE*sizeof(char)),
        fgets(line, MAXCHAR_LINE, tmp_fp)) {
            if (strncmp("ATOM  ", line, 6) && strncmp("HETATM", line, 6)) continue;
            if (line[strlen(line) - 1] == '\n') line[strlen(line) - 1] = '\0'; /* this is to make sure conformer history not end with a new line */
            
            na2++;
            atoms2 = realloc(atoms2,na2*sizeof(ATOM));
            atoms2[na2-1] = pdbline2atom(line);
        }
    }
    else {
        fp2 = fopen(argv[2],"r");
        while ( memset(line,0,MAXCHAR_LINE*sizeof(char)),
        fgets(line, MAXCHAR_LINE, fp2)) {
            if (strncmp("ATOM  ", line, 6) && strncmp("HETATM", line, 6)) continue;
            if (line[strlen(line) - 1] == '\n') line[strlen(line) - 1] = '\0'; /* this is to make sure conformer history not end with a new line */
            
            na2++;
            atoms2 = realloc(atoms2,na2*sizeof(ATOM));
            atoms2[na2-1] = pdbline2atom(line);
        }
    }
    
    for (ia1=0;ia1<na1;ia1++) {
        strncpy(elem_name,atoms1[ia1].name,2);elem_name[2]='\0';
        if (strchr("1234567890",elem_name[0])) elem_name[0] = ' ';
        
        for (ielem1=0;ielem1<nelem1;ielem1++) {
            if (!strcmp(elem_name,elems1[ielem1].name)) break;
        }
        if (ielem1<nelem1) {
            elems1[ielem1].na++;
            elems1[ielem1].ia = realloc(elems1[ielem1].ia,elems1[ielem1].na*sizeof(int));
            elems1[ielem1].ia[elems1[ielem1].na-1] = ia1;
        }
        else {
            nelem1++;
            elems1 = realloc(elems1,nelem1*sizeof(ELEM));
            strcpy(elems1[ielem1].name,elem_name);
            elems1[ielem1].na=1;
            elems1[ielem1].ia = malloc(sizeof(int));
            elems1[ielem1].ia[0] = ia1;

            elems1[ielem1].iswitch = 0;
            elems1[ielem1].jswitch = 0;
        }
    }

    for (ia2=0;ia2<na2;ia2++) {
        strncpy(elem_name,atoms2[ia2].name,2);elem_name[2]='\0';
        if (strchr("1234567890",elem_name[0])) elem_name[0] = ' ';
        
        for (ielem2=0;ielem2<nelem2;ielem2++) {
            if (!strcmp(elem_name,elems2[ielem2].name)) break;
        }
        if (ielem2<nelem2) {
            elems2[ielem2].na++;
            elems2[ielem2].ia = realloc(elems2[ielem2].ia,elems2[ielem2].na*sizeof(int));
            elems2[ielem2].ia[elems2[ielem2].na-1] = ia2;
        }
        else {
            nelem2++;
            elems2 = realloc(elems2,nelem2*sizeof(ELEM));
            strcpy(elems2[ielem2].name,elem_name);
            elems2[ielem2].na=1;
            elems2[ielem2].ia = malloc(sizeof(int));
            elems2[ielem2].ia[0] = ia2;
            
            elems2[ielem2].iswitch = 0;
            elems2[ielem2].jswitch = 0;
        }
    }
    
    nelem_match = nelem1;
    for (ielem1=0;ielem1<nelem_match;ielem1++) {
        for (ielem2=ielem1;ielem2<nelem2;ielem2++) {
            if (!strcmp(elems1[ielem1].name,elems2[ielem2].name)) {
                break;
            }
        }
        if (ielem2 == nelem2) {
            elem_buffer = elems1[ielem1];
            for (jelem1=ielem1;jelem1<nelem1-1;jelem1++) {
                elems1[jelem1] = elems1[jelem1+1];
            }
            elems1[nelem1-1] = elem_buffer;
            ielem1--;
            nelem_match--;
            continue;
        }
        
        na_ielem = elems1[ielem1].na;
        if (elems2[ielem2].na > na_ielem) na_ielem = elems2[ielem2].na;
        for (jelem1=0;jelem1<ielem1;jelem1++) {
            na_jelem = elems1[jelem1].na;
            if (elems2[jelem1].na > na_jelem) na_jelem = elems2[jelem1].na;
            if (na_ielem < na_jelem) break;
        }
        
        elem_buffer = elems1[ielem1];
        for (kelem1=ielem1;kelem1>jelem1;kelem1--) {
            elems1[kelem1] = elems1[kelem1-1];
        }
        elems1[jelem1]=elem_buffer;
        
        elem_buffer = elems2[ielem2];
        for (kelem2=ielem2;kelem2>jelem1;kelem2--) {
            elems2[kelem2] = elems2[kelem2-1];
        }
        elems2[jelem1]=elem_buffer;
    }
    
    int total_matchable = 0;
    for (ielem1=0;ielem1<nelem_match;ielem1++) {
        na_ielem = elems1[ielem1].na;
        if (elems2[ielem2].na < na_ielem) na_ielem = elems2[ielem2].na;
        total_matchable += na_ielem;
    }
        
    sprintf(line," File \"%s\":",argv[1]);
    while(strlen(line) < 40) strcat(line," ");
    printf("\n%s File \"%s\":\n",line,argv[2]);
    
    printf(" %2d kind(s) of element:                | %2d kind(s) of element:\n\n",nelem1,nelem2);
    for (ielem1=0;ielem1<nelem_match;ielem1++) {
        printf("%2d %s atom(s)                          | %2d %s atom(s)\n",elems1[ielem1].na,elems1[ielem1].name,elems2[ielem1].na,elems2[ielem1].name);
    }
    for (ielem1=nelem_match;ielem1<nelem1;ielem1++) {
        printf("%2d %s atom(s)                          |\n",elems1[ielem1].na,elems1[ielem1].name);
    }
    for (ielem2=nelem_match;ielem2<nelem2;ielem2++) {
        printf("                                       | %2d %s atom(s)\n",elems2[ielem2].na,elems2[ielem2].name);
    }
    
    printf("\nGuessed connectivity in file %s\n",argv[1]);
    for (ia1=0;ia1<na1;ia1++) {
        printf("(%4d)%5d %4s%c%3s %c%4d%c|",
        ia1,atoms1[ia1].serial,atoms1[ia1].name,atoms1[ia1].altLoc,atoms1[ia1].resName,atoms1[ia1].chainID,atoms1[ia1].resSeq,atoms1[ia1].iCode);
        for (ja1=0;ja1<na1;ja1++) {
            if (ia1==ja1) continue;
            if (atoms1[ia1].name[1]=='H' && atoms1[ja1].name[1]=='H') continue;
            if (dvv(atoms1[ia1].xyz, atoms1[ja1].xyz)<BOND_THR) {
                printf("(%4d)%5d %4s|",
                ja1,atoms1[ja1].serial,atoms1[ja1].name);
            }
        }
        printf("\n");
    }

    int counter=0;
    while(1){
        counter++;
        
        natom_match = 0;
        int change_order=1;
        for (ielem1=0;ielem1<nelem_match;ielem1++) {
            if (natom_match > 3) {
                break;
            }
            katom_match = natom_match;
            
            /* switching orders of ia list for this element */
            if (elems1[ielem1].na>=elems2[ielem1].na) {
                int ia_buffer = elems1[ielem1].ia[elems1[ielem1].iswitch];
                elems1[ielem1].ia[elems1[ielem1].iswitch] = elems1[ielem1].ia[elems1[ielem1].jswitch];
                elems1[ielem1].ia[elems1[ielem1].jswitch] = ia_buffer;
                if (change_order) {
                    change_order = 0;
                    elems1[ielem1].jswitch++;
                    if (elems1[ielem1].jswitch == elems1[ielem1].na) {
                        elems1[ielem1].jswitch = 0;
                        
                        elems1[ielem1].iswitch++;
                        if (elems1[ielem1].iswitch == elems1[ielem1].na) {
                            elems1[ielem1].iswitch = 0;
                            change_order = 1;
                        }
                    }
                    /*
                    for (ia1=0; ia1<elems1[ielem1].na; ia1++) {
                        printf("%3d,",elems1[ielem1].ia[ia1]);
                    }
                    printf("|%d,%d,%d\n",ielem1,elems1[ielem1].iswitch,elems1[ielem1].jswitch);
                    */
                }
            }
            else {
                int ia_buffer = elems2[ielem1].ia[elems2[ielem1].iswitch];
                elems2[ielem1].ia[elems2[ielem1].iswitch] = elems2[ielem1].ia[elems2[ielem1].jswitch];
                elems2[ielem1].ia[elems2[ielem1].jswitch] = ia_buffer;
                if (change_order) {
                    change_order = 0;
                    elems2[ielem1].jswitch++;
                    if (elems2[ielem1].jswitch == elems2[ielem1].na) {
                        elems2[ielem1].jswitch = 0;
                        
                        elems2[ielem1].iswitch++;
                        if (elems2[ielem1].iswitch == elems2[ielem1].na) {
                            elems2[ielem1].iswitch = 0;
                            change_order = 1;
                        }
                    }
                    /*
                    for (ia1=0; ia1<elems2[ielem1].na; ia1++) {
                        printf("%3d,",elems2[ielem1].ia[ia1]);
                    }
                    printf("|%d,%d,%d\n",ielem1,elems2[ielem1].iswitch,elems2[ielem1].jswitch);
                    */
                }
            }
            
            /*
            if (change_order) {
                change_order = 0;
                if (elems1[ielem1].na>=elems2[ielem1].na) {
                    elems1[ielem1].jswitch++;
                    if (elems1[ielem1].jswitch == elems1[ielem1].na) {
                        elems1[ielem1].jswitch = 0;
                        
                        elems1[ielem1].iswitch++;
                        if (elems1[ielem1].iswitch == elems1[ielem1].na) {
                            elems1[ielem1].iswitch = 0;
                            change_order = 1;
                        }
                    }
                    for (ia1=0; ia1<elems1[ielem1].na; ia1++) {
                        printf("%3d,",elems1[ielem1].ia[ia1]);
                    }
                    printf("|%d,%d,%d\n",ielem1,elems1[ielem1].iswitch,elems1[ielem1].jswitch);
                }
                else {
                    elems2[ielem1].jswitch++;
                    if (elems2[ielem1].jswitch == elems2[ielem1].na) {
                        elems2[ielem1].jswitch = 0;
                        
                        elems2[ielem1].iswitch++;
                        if (elems2[ielem1].iswitch == elems2[ielem1].na) {
                            elems2[ielem1].iswitch = 0;
                            change_order = 1;
                        }
                    }
                    for (ia1=0; ia1<elems2[ielem1].na; ia1++) {
                        printf("%3d,",elems2[ielem1].ia[ia1]);
                    }
                    printf("|%d,%d,%d\n",ielem1,elems2[ielem1].iswitch,elems2[ielem1].jswitch);
                }
            }
            
            if (change_order) {
                change_order = 0;
                if (elems1[ielem1].na<elems1[ielem1].na) {
                    elems1[ielem1].jswitch++;
                    if (elems1[ielem1].jswitch == elems1[ielem1].na) {
                        elems1[ielem1].jswitch = 0;
                        
                        elems1[ielem1].iswitch++;
                        if (elems1[ielem1].iswitch == elems1[ielem1].na) {
                            elems1[ielem1].iswitch = 0;
                            change_order = 1;
                        }
                    }
                }
                else {
                    elems2[ielem1].jswitch++;
                    if (elems2[ielem1].jswitch == elems2[ielem1].na) {
                        elems2[ielem1].jswitch = 0;
                        
                        elems2[ielem1].iswitch++;
                        if (elems2[ielem1].iswitch == elems2[ielem1].na) {
                            elems2[ielem1].iswitch = 0;
                            change_order = 1;
                        }
                    }
                    printf("%d,%d,%d\n",ielem1,elems2[ielem1].iswitch,elems2[ielem1].jswitch);
                }
            }
            */

            na_ielem = elems1[ielem1].na;
            if (elems2[ielem1].na < na_ielem) na_ielem = elems2[ielem1].na;
            natom_match += na_ielem;
            
            atoms_to_match1 = realloc(atoms_to_match1,natom_match*sizeof(ATOM));
            atoms_to_match2 = realloc(atoms_to_match2,natom_match*sizeof(ATOM));
            
            for (ia1=0; ia1<elems1[ielem1].na&&ia1<elems2[ielem1].na; ia1++,katom_match++) {
                atoms_to_match1[katom_match] = atoms1[elems1[ielem1].ia[ia1]];
                atoms_to_match2[katom_match] = atoms2[elems2[ielem1].ia[ia1]];
            }
        }
        
        //if (change_order) exit(-1);
        //printf("Trying matching first %2d atoms...\n",natom_match);
        
        int ia;
        for (ia=2;ia<natom_match;ia++) {
            GEOM op = geom_3v_onto_3v(atoms_to_match2[0].xyz,atoms_to_match2[1].xyz,atoms_to_match2[ia].xyz,atoms_to_match1[0].xyz,atoms_to_match1[1].xyz,atoms_to_match1[ia].xyz);
            int natom_matched=0;
            for (ia2=0;ia2<natom_match;ia2++) {
                VECTOR xyz2 = atoms_to_match2[ia2].xyz;
                geom_apply(op, &xyz2);

                if (dvv(atoms_to_match1[ia2].xyz, xyz2) < MATCH_THR) {
                    natom_matched++;
                    //printf("%d: %8.3f,%8.3f,%8.3f  %8.3f,%8.3f,%8.3f, %8.3f\n",ia2,atoms_to_match1[ia2].xyz.x,atoms_to_match1[ia2].xyz.y,atoms_to_match1[ia2].xyz.z,atoms_to_match2[ia2].xyz.x,atoms_to_match2[ia2].xyz.y,atoms_to_match2[ia2].xyz.z,dvv(atoms_to_match1[ia2].xyz,atoms_to_match2[ia2].xyz));
                }
            }
            
            if (natom_matched>=3) {
                natom_matched = 0;
                for (ia2=0;ia2<na2;ia2++) atoms2[ia2].on = 0;
                for (ia1=0;ia1<na1;ia1++) {
                    for (ia2=0;ia2<na2;ia2++) {
                        if (atoms1[ia1].name[1] != atoms2[ia2].name[1]) continue;
                        VECTOR xyz2 = atoms2[ia2].xyz;
                        geom_apply(op, &xyz2);
                        if (dvv(atoms1[ia1].xyz,xyz2) <= MATCH_THR) {natom_matched++;break;}
                    }
                }
                if (natom_matched<total_matchable/2.) continue;
                
                /* Last fitting */
                natom_matched = 0;
                for (ia2=0;ia2<na2;ia2++) atoms2[ia2].on = 0;
                for (ia1=0;ia1<na1;ia1++) {
                    int matched_ia2 = -1;
                    float matched_dist = MATCH_THR;
                    for (ia2=0;ia2<na2;ia2++) {
                        if (atoms1[ia1].name[1] != atoms2[ia2].name[1]) continue;
                        if (atoms2[ia2].on) continue;
                        VECTOR xyz2 = atoms2[ia2].xyz;
                        geom_apply(op, &xyz2);
                        if (dvv(atoms1[ia1].xyz,xyz2) >= matched_dist) continue;
                        matched_ia2 = ia2;
                        matched_dist = dvv(atoms1[ia1].xyz, xyz2);
                    }
                    if (matched_ia2 != -1) {
                        atoms2[matched_ia2].on = 1;
                        natom_matched++;
                        atoms_to_match1 = realloc(atoms_to_match1,natom_matched*sizeof(ATOM));
                        atoms_to_match2 = realloc(atoms_to_match2,natom_matched*sizeof(ATOM));
                        atoms_to_match1[natom_matched-1] = atoms1[ia1];
                        atoms_to_match2[natom_matched-1] = atoms2[matched_ia2];
                        geom_apply(op, &atoms_to_match2[natom_matched-1].xyz);
                    }
                }
                GEOM op_correct = geom_fit(atoms_to_match1,atoms_to_match2,natom_matched);
                
                printf("\n  Matched atoms:\n");
                sprintf(line," File \"%s\":",argv[1]);
                while(strlen(line) < 40) strcat(line," ");
                printf("\n%s File \"%s\":\n",line,argv[2]);
                natom_matched = 0;
                for (ia2=0;ia2<na2;ia2++) atoms2[ia2].on = 0;
                for (ia1=0;ia1<na1;ia1++) {
                    int matched_ia2 = -1;
                    float matched_dist = MATCH_THR;
                    for (ia2=0;ia2<na2;ia2++) {
                        if (atoms1[ia1].name[1] != atoms2[ia2].name[1]) continue;
                        if (atoms2[ia2].on) continue;
                        //printf("natom_matched %d %d\n",natom_matched,na2);
                        VECTOR xyz2 = atoms2[ia2].xyz;
                        geom_apply(op, &xyz2);
                        geom_apply(op_correct, &xyz2);
                        if (dvv(atoms1[ia1].xyz,xyz2) >= matched_dist) continue;
                        matched_ia2 = ia2;
                        matched_dist = dvv(atoms1[ia1].xyz, xyz2);
                    }
                    if (matched_ia2 != -1) {
                        atoms2[matched_ia2].on = 1;
                        printf("(%4d)%5d %4s%c%3s %c%4d%c              |%5d %4s%c%3s %c%4d%c\n",
                        ia1,atoms1[ia1].serial,atoms1[ia1].name,atoms1[ia1].altLoc,atoms1[ia1].resName,atoms1[ia1].chainID,atoms1[ia1].resSeq,atoms1[ia1].iCode,
                        atoms2[matched_ia2].serial,atoms2[matched_ia2].name,atoms2[matched_ia2].altLoc,atoms2[matched_ia2].resName,atoms2[matched_ia2].chainID,atoms2[matched_ia2].resSeq,atoms2[matched_ia2].iCode);
                        
                        natom_matched++;
                    }
                    else {
                        printf("(%4d)%5d %4s%c%3s %c%4d%c\n",
                        ia1,atoms1[ia1].serial,atoms1[ia1].name,atoms1[ia1].altLoc,atoms1[ia1].resName,atoms1[ia1].chainID,atoms1[ia1].resSeq,atoms1[ia1].iCode);
                    }
                }
                printf("Number of atoms matched: %d out of %d\n",natom_matched,na1);
                
                printf("Check matched.pdb file for rotated PDB.\n");
                FILE *out_fp = fopen("matched.pdb","w");
                for (ia2=0;ia2<na2;ia2++) {
                    geom_apply(op, &atoms2[ia2].xyz);
                    geom_apply(op_correct, &atoms2[ia2].xyz);
                    fprintf(out_fp, "ATOM  %5d %4s%c%3s %c%04d%c   %8.3f%8.3f%8.3f %7.3f      %6.3f      %-10s\n",
                                    atoms2[ia2].serial,
                                    atoms2[ia2].name,
                                    atoms2[ia2].altLoc,
                                    atoms2[ia2].resName,
                                    atoms2[ia2].chainID,
                                    atoms2[ia2].resSeq,
                                    atoms2[ia2].iCode,
                                    atoms2[ia2].xyz.x,
                                    atoms2[ia2].xyz.y,
                                    atoms2[ia2].xyz.z,
                                    atoms2[ia2].rad,
                                    atoms2[ia2].crg,
                                    atoms2[ia2].history);
                }
                exit(0);
            }
            //printf("Counter %d\n",counter);
        }
    }
    //if (ielem1 == 0) printf("Too much freedom for matching!\n");
    
    /*
    for (ielem2=0;ielem2<nelem2;ielem2++) {
        for (jelem1=0;jelem1<nelem1;jelem1++) {
    }}
    */
    return 0;
}

GEOM geom_fit(ATOM *atoms1, ATOM *atoms2, int na) {
    GEOM op;
    int ia;
    VECTOR sum1,sum2;
    VECTOR *tmp_v = malloc(na*sizeof(VECTOR));
    sum1.x=0;sum1.y=0;sum1.z=0;
    sum2.x=0;sum2.y=0;sum2.z=0;
    for (ia=0;ia<na;ia++) {
        sum1 = vector_vplusv(sum1,atoms1[ia].xyz);
        sum2 = vector_vplusv(sum2,atoms2[ia].xyz);
    }
    VECTOR move = vector_vminusv(sum1,sum2);
    move = vector_rescale(move,1./na);
    for (ia=0;ia<na;ia++) {
        tmp_v[ia] = vector_vplusv(atoms2[ia].xyz,move);
    }
    geom_reset(&op);
    geom_move(&op, move);
    
    VECTOR origin = vector_rescale(sum1, 1./na);
    float rmsd=0.,avg_r=0.;
    for (ia=0;ia<na;ia++) {
        rmsd += ddvv(atoms1[ia].xyz,tmp_v[ia]);
        avg_r += ddvv(atoms1[ia].xyz,origin);
    }
    rmsd = sqrt(rmsd/na);
    avg_r = sqrt(avg_r / na);
    
    float step = rmsd/avg_r;
    
    while (step*avg_r>0.001) {
        VECTOR  torq;
        torq.x=0;torq.y=0;torq.z=0;
        for (ia=0;ia<na;ia++) {
            VECTOR t  = vector_vminusv(atoms1[ia].xyz,tmp_v[ia]);
            VECTOR r2 = vector_vminusv(tmp_v[ia],origin);
            torq = vector_vplusv(torq,vector_vxv(r2,t));
        }
        LINE axis = line_2v(origin,vector_vplusv(origin,torq));
        geom_roll(&op, step, axis);
        
        GEOM tmp_op;
        geom_reset(&tmp_op);
        geom_roll(&tmp_op, step, axis);
        for (ia=0;ia<na;ia++) {
            geom_apply(tmp_op, &tmp_v[ia]);
        }
        
        float new_rmsd=0.;
        for (ia=0;ia<na;ia++) {
            new_rmsd += ddvv(atoms1[ia].xyz,tmp_v[ia]);
        }
        new_rmsd = sqrt(new_rmsd/na);
        step = (rmsd-new_rmsd)/avg_r;
        //printf("step %8.3f\n",step);
        if (step<0) step = -step;
        
        rmsd = new_rmsd;
    }
    
    free(tmp_v);
    return op;
}
