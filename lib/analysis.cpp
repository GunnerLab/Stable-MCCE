//#include <glib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <zlib.h>
//#include "mcce.h"
extern "C" {
#include "mcce.h"
}

#include <cstdlib>
#include <cstring>
#include <cmath>
#include <set>
#include <vector>
#include <string>

#include <iostream>
#include <fstream>
#include <map>
#include <iomanip>

using namespace std;

/* public variable  */
const string PDBF = "step2_out.pdb";
const string HBOUT = "hb.dat";
float DFAR;       // the distance used to determine a hydrogen bond
float DNEAR;       // default: 1.2 < d < 3.2
float ANGCUT;

static RES conflist;
FILE *h_fp;
//const float PI = 3.1415926;


/* function */
int hbond_matrix();
int hbond_network();
int is_hb(CONF *conf1, CONF *conf2);
int load_head3lst();


int analysis()
{   
    
    /* Load step2_out.pdb and get residue list for the microstate output  */
    /*
    if (env.get_ms_gold){
        printf("   Load step2_out.pdb and get ms_gold ...\n"); fflush(stdout);
        if (ms_gold()) {db_close(); return USERERR;}
        else printf("File ms_gold is obtained.\n\n");
    }
    */


    /*Load step3_out.pdb and get hydrogen bond matrix */
    if (env.get_hbond_matrix){
        printf("   Load step2_out.pdb and get hydrogen bond matrix ...\n"); fflush(stdout);
        if (hbond_matrix()) {
	printf("   Fatal error detected in hbond_matrix().\n"); 
	return USERERR;
	}
        else printf("   Files hb.dat and hah.txt are obtained.\n\n");
    }
    
    /*Load hb.dat and ms.dat and combine them to get hydrogen bond network */
    if (env.get_hbond_network){
        printf("   Load hb.dat and ms.dat and combine them to get hydrogen bond network ...\n"); fflush(stdout);
        if (hbond_network()) {
	printf("   Fatal error detected in hbond_network().\n");
	return USERERR;
	}
        else printf("   Files hb.txt are obtained.\n\n");
    }

    return 0;
}





int hbond_matrix()
{
/*	db_open();
	if (init()) {
		db_close();
		printf("Help message: double check file \"run.prm\" in current directory.\n");
		return USERERR;
	}
*/

	DNEAR=env.hbond_lower_limit;
        DFAR=env.hbond_upper_limit;
        ANGCUT=env.hbond_ang_cutoff;
	printf("   Hydrogen Bond are created when distance is between %.3f and %.3f, angle is equal or larger than %f.\n", DNEAR, DFAR, ANGCUT);
	FILE *fp;
	PROT prot;
	if (!(fp=fopen(STEP2_OUT, "r"))) {
		printf("   \"No step 2 output \"%s\".\n", STEP2_OUT);
		return USERERR;
	}
	prot = load_pdb(fp);
	if (prot.n_res == 0) {
		printf("   There are errors in pdb file, quiting ...\n");
		return USERERR;
	}

	id_conf(prot);
	get_connect12(prot);

	load_head3lst();
	int n_conf = conflist.n_conf;

	int i, j, k;
	//    for (k=0; k<n_conf; k++) printf("%s\t%d\n", conflist.conf[k].uniqID, conflist.conf[k].iConf);

	for (i=0; i<prot.n_res; i++) {
		for (j=1; j<prot.res[i].n_conf; j++) {
			strncpy(prot.res[i].conf[j].confName, prot.res[i].conf[j].uniqID, 5);
			prot.res[i].conf[j].confName[5] = '\0';
			// adjust the conf id to avoid id diff between step2 and head3.lst
			for (k=0; k<n_conf; k++) {
				if (!strcmp(conflist.conf[k].uniqID, prot.res[i].conf[j].uniqID)) {
					prot.res[i].conf[j].iConf = conflist.conf[k].iConf;
					break;
				}
			}
		}
	}



	h_fp = fopen("hah.txt", "w");

	// hb_fp is a binary file to store the hbpw matrix, the first 4 bytes give the conformer number n_conf
	FILE *hb_fp = fopen(HBOUT.c_str(), "wb");
	fwrite(&n_conf, sizeof(int), 1, hb_fp);

	int i_res, j_res, i_conf, j_conf;

	// hbpw is a matrix to store the hb connection between each two conf, the elem is 0 or 1
	FILE *fp_res = fopen("reshbond.txt", "w");
	char *resHb = (char *) calloc(prot.n_res * prot.n_res, sizeof(char));

	char *hbpw = (char *) calloc(n_conf * n_conf, sizeof(char));
	for (i_res=0; i_res<prot.n_res; i_res++) {
		for (i_conf=1; i_conf<prot.res[i_res].n_conf; i_conf++) {
			for (j_res=0; j_res<prot.n_res; j_res++) {
				if (j_res == i_res) continue;
				for (j_conf=1; j_conf<prot.res[j_res].n_conf; j_conf++) {

					/* there is hbond from donor i_conf to accepter j_conf */
					if (is_hb(&prot.res[i_res].conf[i_conf], &prot.res[j_res].conf[j_conf])) {
						hbpw[prot.res[i_res].conf[i_conf].iConf * n_conf + prot.res[j_res].conf[j_conf].iConf] = 1;
						resHb[i_res * prot.n_res + j_res] = 1;
					}
				}
			}
		}
	}
	/*
    for (i=0; i<n_conf; i++) {
        for (j=0; j<i; j++) {
            hbpw[i * n_conf + j] = hbpw[j * n_conf + i];
        }
    }
	 */
	fwrite(hbpw, sizeof(char), n_conf * n_conf, hb_fp);
	fclose(hb_fp);

	set <int> resInHbNet;
	for (i_res=0; i_res<prot.n_res; i_res++) {
		for (j_res=0; j_res<prot.n_res; j_res++) {
			if (j_res == i_res) continue;
			if (resHb[i_res * prot.n_res + j_res] == 1) {
				fprintf(fp_res, "%s%c%04d\t%s%c%04d\n", prot.res[i_res].resName, prot.res[i_res].chainID, prot.res[i_res].resSeq,
						prot.res[j_res].resName, prot.res[j_res].chainID, prot.res[j_res].resSeq);
				resInHbNet.insert(i_res);
				resInHbNet.insert(j_res);
			}
		}
	}
	fclose(fp_res);

	FILE *fp_resInHbNet = fopen("resInHbNet.txt", "w");
	for (i_res=0; i_res<prot.n_res; i_res++) {
		if (resInHbNet.find(i_res) != resInHbNet.end()) {
			fprintf(fp_resInHbNet, "%s%c%04d\n", prot.res[i_res].resName, prot.res[i_res].chainID, prot.res[i_res].resSeq);
		}
	}
	fclose(fp_resInHbNet);

	printf("   n_res: %d\n", prot.n_res);
//	db_close();
	return 0;
}

int is_hb( CONF *conf1, CONF *conf2)
{
	STRINGS Datoms,Aatoms;
	int iD, iA, Dseq, Aseq;
	float d=0.0;
	float ang = 0.0;
	int isHb = 0;
	/** To establish hydrogen bond between two conformers, the following criteria must be met:
	 * 1. one conformer is h donor, the other is acceptor.
	 * 2. the distance between one donor atom and one acceptor atom is between DNEAR and DFAR.
	 * 3. the angle of h bond should be no less than 90.
	 */
	if (!param_get((char *) "HDONOR", conf1->confName,(char *) "", &Datoms)) {
//		printf("conf1->confName: %s, conf2->confName: %s, Datoms.n: %d\n", conf1->confName,
//			conf2->confName, Datoms.n);                    //test--Cai
		if (!param_get((char *)"HACCEPT", conf2->confName,(char *) "", &Aatoms)) {
			for (iD=0; iD<Datoms.n; iD++) {
				param_get((char *)"IATOM", conf1->confName, Datoms.strings[iD], &Dseq);
//				if (!strcmp(conf1->confName, "HOH-1")) {
//					printf("Donor: %s, Datoms.n: %d, Dseq: %d, atom name %s\n", conf1->confName,
//						Datoms.n, Dseq, Datoms.strings[iD]);
//				}
				for (iA=0; iA<Aatoms.n; iA++) {
					param_get((char *) "IATOM", conf2->confName, Aatoms.strings[iA], &Aseq);
					d = ddvv(conf1->atom[Dseq].xyz, conf2->atom[Aseq].xyz);

//					d = ddvv(conf1->atom[Dseq].connect12[0]->xyz, conf2->atom[Aseq].xyz);  //*****
					if (d > DNEAR * DNEAR && d < DFAR * DFAR) {
						ang = avv(vector_vminusv((conf1->atom[Dseq].connect12[0])->xyz, conf1->atom[Dseq].xyz),
								vector_vminusv(conf2->atom[Aseq].xyz, conf1->atom[Dseq].xyz)) * 180.0 / env.PI;
						if (abs(ang) < ANGCUT) continue;
						fprintf(h_fp, "%s\t%s\t%s~%s--%s\t%.2f\t%.0f\n", conf1->uniqID, conf2->uniqID,
								(conf1->atom[Dseq].connect12[0])->name, conf1->atom[Dseq].name, conf2->atom[Aseq].name, sqrt(d), ang);
						isHb = 1;
					}
				}
			}
		}
	}
	
	return isHb;
}

int load_head3lst()
{
	FILE *fp;
	char sbuff[MAXCHAR_LINE];
	char stemp[MAXCHAR_LINE];
	CONF conf_temp;
	char notfound;
	int iconf, ires;
	int kr;
	int counter;

	conflist.n_conf = 0;
	conflist.conf   = NULL;

	if (!(fp=fopen(FN_CONFLIST3, "r"))) {
		printf("   FATAL: Can't open file %s\n", FN_CONFLIST3);
		return USERERR;
	}
	fgets(sbuff, sizeof(sbuff), fp); /* skip the first line */
	counter = 0;
	while(fgets(sbuff, sizeof(sbuff), fp)) {
		/* load this line to a conf template */
		if (strlen(sbuff) < 20) continue;
		sscanf(sbuff, "%d %s %c %f %f %f %f %d %d %f %f %f %f %f %f %s", &conf_temp.iConf,
				conf_temp.uniqID,
				&conf_temp.on,
				&conf_temp.occ,
				&conf_temp.netcrg,
				&conf_temp.Em,
				&conf_temp.pKa,
				&conf_temp.e,
				&conf_temp.H,
				&conf_temp.E_vdw0,
				&conf_temp.E_vdw1,
				&conf_temp.E_tors,
				&conf_temp.E_epol,
				&conf_temp.E_dsolv,
				&conf_temp.E_extra,
				conf_temp.history);

		conf_temp.E_TS = 0.0; /* initialize entropy effect at the time of loading conflist */

		/* rescale */
		conf_temp.E_vdw0 *= env.scale_vdw0;
		conf_temp.E_vdw1 *= env.scale_vdw1;
		conf_temp.E_epol *= env.scale_ele;
		conf_temp.E_tors *= env.scale_tor;
		conf_temp.E_dsolv*= env.scale_dsolv;

		strncpy(conf_temp.resName, conf_temp.uniqID, 3); conf_temp.resName[3] = '\0';
		strncpy(conf_temp.confName, conf_temp.uniqID, 5); conf_temp.confName[5] = '\0';
		conf_temp.chainID = conf_temp.uniqID[5];
		strncpy(stemp, conf_temp.uniqID+6, 4); stemp[4] = '\0';
		conf_temp.resSeq = atoi(stemp);
		conf_temp.iCode = conf_temp.uniqID[10];
		conf_temp.n_atom = 0;
		if (conf_temp.on == 't' || conf_temp.on == 'T') conf_temp.on = 't';
		else conf_temp.on = 'f';
		conf_temp.iConf = counter;
		/* creating conflist */
		iconf = ins_conf(&conflist, conflist.n_conf, 0);
		cpy_conf(&conflist.conf[iconf], &conf_temp);
		counter++;
	}

	return 0;
}




int hbond_network()
{
	time_t time_start = time(NULL);

	ifstream iHbfile("hb.dat", ios::in | ios::binary);
        if (!iHbfile.is_open()) {
            cerr << "   Can't open file hb.dat." << endl;
            return -1;
        }
	int n_conf;
	vector<vector<char> > hinter;

	iHbfile.read((char *) &n_conf, sizeof(int));
	hinter.resize(n_conf);
	for (int i=0; i<n_conf; i++) hinter[i].resize(n_conf);

	for (int i=0; i<n_conf; i++) {
		for (int j=0; j<n_conf; j++) {
			iHbfile.read((char *) &hinter[i][j], 1);
		}
	}
	iHbfile.close();

	int n_spe, count;
	unsigned short conf;
	char buffer[9];
        char method[9];
	//double H_state, Hsq_state;
        double H_state, Hav_state;
	vector <string> resName;
	vector <unsigned short> iconf;

	ifstream iMsfile("ms.dat", ios::in | ios::binary);

	// n_spe gives the number of special residues
	iMsfile.read((char *) &n_spe, 4);
	cout << "   There are " << n_spe << " key residues." << endl;
	// resName are the names of the special residues
	for (int i_spe=0; i_spe<n_spe; i_spe++) {
		iMsfile.read(buffer, 8);
		buffer[8] = '\0';
		resName.push_back(string(buffer));
	}

        iMsfile.read(method, 9);
        cout << "   The microstates are obtained from: " << method << "." << endl;
	vector<vector<int> > resinter;
	resinter.resize(n_spe);
	for (size_t i=0; i<n_spe; i++) resinter[i].resize(n_spe);

	int totalState = 0;
	int totalRecords = 0;

	const int PRINT_INTERVAL = 10000;
	while (iMsfile.read((char *) &conf, 2)) {
		iconf.push_back(conf);
		for (int i_spe=1; i_spe<n_spe; i_spe++) {
			iMsfile.read((char *) &conf, 2);
			iconf.push_back(conf);
		}
		iMsfile.read((char *) &H_state, 8);
		//iMsfile.read((char *) &Hsq_state, 8);
                iMsfile.read((char *) &Hav_state, 8);
		iMsfile.read((char *) &count, 4);

		totalRecords++;
		totalState += count;
//		if (totalRecords % PRINT_INTERVAL == 0) cout << totalRecords << " records have been loaded." << endl;

		for (int i_res=0; i_res<n_spe; i_res++) {
			for (int j_res=0; j_res<n_spe; j_res++) {
				if (j_res == i_res) continue;
				resinter[i_res][j_res] += hinter[iconf[i_res]][iconf[j_res]] * count;
			}
		}
		iconf.clear();
	}
	iMsfile.close();
	cout << "   "  << totalRecords << " records and " << totalState << " states has been loaded." << endl;
	cout << "   Total time to load ms.dat: " << time(NULL) - time_start << "." << endl;
	const float THRESHOLD_TO_WRITE = 0.001;

	ofstream ofile("hb.txt", ios::out);
	for (int i=0; i<n_spe; i++) {
		for (int j=0; j<n_spe; j++) {
			if ((float) resinter[i][j]/totalState >= THRESHOLD_TO_WRITE) {
				ofile << resName[i] << '\t' << resName[j] << '\t' << fixed << setprecision(3) << ((float) resinter[i][j]/totalState) << endl;;
			}
		}
	}
	ofile.close();



	cout << "   Total time: " << time(NULL) - time_start << "." << endl;

    return 0;
}

