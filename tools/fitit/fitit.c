#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mcce.h"

struct STAT {
    float a;
    float b;
    float chi2;
};

struct STAT fit(float a, float b);
float score(float v[]);
void dhill(float **p, float *y, int ndim, float ftol, float (*funk)(float []), int *nfunk);

int Nx;
float *xp, *yp;   /* titration points */

int main()
{   FILE *fp, *fp2;
    int i, j;
    int N_res;
    struct STAT stat;   /* a structure of statistcs */
    char **shead;
    float **ysp;
    float a, b, mid;
    float n;
    char sbuff[128], sbuff2[512];

    /*<<< Get titration points >>>*/

    /* read in sum_crg.out file */
    if (!(fp = fopen("sum_crg.out", "r"))) {
        printf("File sum_crg.out not found\n");
        return 0; 
    }
    else {
       fgets(sbuff, sizeof(sbuff), fp);

       /* titration types */
       if (strstr(sbuff, "pH") || strstr(sbuff, "ph")) env.titr_type = 'p';
       else env.titr_type = 'e';

       /*--- Determine the titration points, x values ---*/
       Nx = (strlen(sbuff)-11)/6;
       xp = (float *) malloc(Nx * sizeof(float));     /* store x points */
       yp = (float *) malloc(Nx * sizeof(float));     /* store y points */

       /*--- Convert x points to float ---*/
       for (i=0; i<Nx; i++) {
          strncpy(sbuff2, sbuff+i*6+11, 6); sbuff2[6] = '\0';
          xp[i] = atof(sbuff2);
       }

       N_res = 0;
       fp2 = tmpfile();
       while (fgets(sbuff2, sizeof(sbuff2), fp)) {
         fputs(sbuff2, fp2);
         N_res++;
       }
       fclose(fp);
       N_res -= 4; /* extra lines */
       
       /*--- Load y values ---*/
       ysp = (float **) malloc(N_res * sizeof(float *));
       for (i=0; i<N_res; i++) ysp[i] = (float *) malloc(Nx * sizeof(float));
       shead = (char **) malloc(N_res * sizeof(char *));
       for (i=0; i<N_res; i++) shead[i] = (char *) malloc(20 * sizeof(char));

       rewind(fp2);
       for (i=0; i<N_res; i++) {
         fgets(sbuff2, sizeof(sbuff2), fp2);
         strncpy(shead[i], sbuff2, 10); shead[i][10] = '\0';
         for (j=0; j<Nx; j++) {
            strncpy(sbuff, sbuff2+j*6+11, 6); sbuff[6] = '\0';
            ysp[i][j] = fabs(atof(sbuff));
            /* printf("%8.3f", ysp[i][j]); */
         }
         /* printf("\n"); */
       }
       fclose(fp2);
       
    }

    /*<<< Loop over y values >>>*/
    if (env.titr_type == 'p') {   /* pH titration */
        printf("  pH      ");
    }
    else {      /* Eh titration assumed */
        printf("  Eh      ");
    }
    printf("       pKa/Em  n(slope) 1000*chi2\n");
    for (i=0; i <N_res; i++) {
        /*--- Convert y points to float ---*/

        /* a reasonable guess makes optimization easier */
        a = 0.0;
        mid = 0.6;
        for (j=0; j<Nx; j++) {
            yp[j] = ysp[i][j];
            if (fabs(yp[j] - 0.5) < mid) {
                mid = fabs(yp[j] - 0.5);
                b = xp[j];
            }
            /* printf("%.3f %.3f\n", yp[j], mid); */
        }
        if (mid >= 0.485) {
            if (fabs(yp[0]-yp[Nx-1])>0.5) {  /* jumps form <0.015 to >0.985 */
               printf("%s          pKa titration curve too sharp\n", shead[i]);
            }
            else {                           /* all < 0.015 or all >0.985 */
               printf("%s          pKa or Em out of range   \n", shead[i]);
            }
            continue;
        }

        /* printf("a=%.3f; b=%.3f\n",a,b); */
        stat = fit(a, b);
        if (env.titr_type == 'p') n = fabs(stat.a * 8.617342E-2 * 273.15 / 58.0);
        else if (env.titr_type == 'e') n = fabs(stat.a * 8.617342E-2 * 273.15);
        else n = stat.a;

        if (stat.b < xp[0] || stat.b > xp[Nx-1]) printf("%s        pKa or Em out of range   \n", shead[i]);
        else printf("%s    %9.3f %9.3f %9.3f\n", shead[i], stat.b, n, 1000*stat.chi2);
    }


    free(xp);
    free(yp);
    for (i=0; i<N_res; i++) free(ysp[i]);
    free(ysp);
    for (i=0; i<N_res; i++) free(shead[i]);
    free(shead);
    
    return 0;
}


struct STAT fit(float a, float b)
/* initialize the simplex and set up optimization */
{  float **p;
    float y[3];
    int neval;
    struct STAT result;
    
    p = (float **) malloc(3 * sizeof(float *));
    p[0] = (float *) malloc(2 * sizeof(float));
    p[1] = (float *) malloc(2 * sizeof(float));
    p[2] = (float *) malloc(2 * sizeof(float));
    
    p[0][0] = a;      p[0][1] = b;      y[0] = score(p[0]);
    p[1][0] = a+0.01; p[1][1] = b;      y[1] = score(p[1]);
    p[2][0] = a;      p[2][1] = b+1.0;  y[2] = score(p[2]);
    
    dhill(p, y, 2, 0.0001, &score, &neval);
    
    /* DEBUG
    printf("(%8.3f%8.3f)=%8.3f exit at %d\n", p[0][0], p[0][1], y[0], neval);
    */
    
    /* fit this eq. y = exp(a(x-b))/(1+exp(a(x-b))) */
    
    result.a = p[0][0];
    result.b = p[0][1];
    result.chi2 = y[0];
    
    free(p[0]);
    free(p[1]);
    free(p[2]);
    free(p);
    
    return result;
}


float score(float v[])
{  float S = 0.0;
    float yt;
    float T;
    int i;
    
    for (i=0; i<Nx; i++) {
        T = exp(v[0]*(xp[i]-v[1]));
        yt = T/(1.0+T);
        S += (yt-yp[i])*(yt-yp[i]);
    }
    
    return S;
}


#define TINY 1.0E-10
#define NMAX 5000
#define GET_PSUM for (j=0; j<ndim; j++) {\
                        for (sum=0.0, i=0; i<mpts; i++) sum += p[i][j];\
                       psum[j] = sum;}
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}
                 
void dhill(float **p, float *y, int ndim, float ftol, float (*funk)(float []), int *nfunk)
{  float dhtry(float **p, float *y, float *psum, int ndim, float (*funk)(float []), int ihi, float fac);
    int i, ihi, ilo, inhi,j, mpts = ndim+1;
    float rtol,sum,swap,ysave,ytry,*psum;
    
    
    psum = (float *)malloc(ndim * sizeof(float));
    *nfunk = 0;
    GET_PSUM
    
    for (;;) {
        ilo = 0;
        ihi = y[0] > y[1] ? (inhi = 1,0):(inhi = 0,1);
        for(i=0; i<mpts;i++) {
            if(y[i] <= y[ilo]) ilo=i;
            if(y[i] > y[ihi]) {
                inhi = ihi;
                ihi  = i;
            }
            else if (y[i] > y[inhi] && i != ihi) inhi = i;
        }
        
        /* DEBUG
        for (i=0; i<mpts; i++) {
            printf("%8.3f at (%8.3f, %8.3f)\n", y[i], p[i][0], p[i][1]);
        }
        printf("\n");
        */
        
        rtol = 2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo])+TINY);
        
        if(rtol < ftol || *nfunk >= NMAX) {
            SWAP(y[0],y[ilo])
            for (i=0; i<ndim; i++) SWAP(p[0][i],p[ilo][i])
                break;
        }
        
        *nfunk += 2;
        
        ytry = dhtry(p,y,psum,ndim,funk,ihi,-1.0);
        
        if (ytry <= y[ilo]) ytry=dhtry(p,y,psum,ndim,funk,ihi,2.0);
        else if (ytry >= y[inhi]) {
            ysave = y[ihi];
            ytry=dhtry(p,y,psum,ndim,funk,ihi,0.5);
            if (ytry >= ysave) {
                for (i=0; i<mpts; i++) {
                    if (i !=ilo) {
                        for (j=0; j<ndim; j++)
                            p[i][j] = psum[j] = 0.5*(p[i][j]+p[ilo][j]);
                        y[i] = (*funk)(psum);
                    }
                }
                *nfunk += ndim;
                GET_PSUM
            }
        }
        else --(*nfunk);
    }
    free(psum);
}

float dhtry(float **p, float *y, float *psum, int ndim, float (*funk)(float []), int ihi, float fac)
{  int j;
    float fac1, fac2, ytry, *ptry;
    
    ptry = (float *)malloc(ndim*sizeof(float));
    fac1 = (1.0-fac)/ndim;
    fac2 = fac1 - fac;
    for (j=0; j<ndim; j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
    ytry = (*funk)(ptry);   /* evaluate the function at the trial point */
    if (ytry < y[ihi]) {    /* if it's better than the highest, then replace the highest */
        y[ihi] = ytry;
        for (j=0; j<ndim; j++) {
            psum[j] += ptry[j] - p[ihi][j];
            p[ihi][j] = ptry[j];
        }
    }
    free(ptry);
    return ytry;
}

