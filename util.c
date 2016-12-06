#include <stdlib.h>
#include <stdio.h>

#include "def.h"
#include "methods.h"

/************************************************************************
 *  Memory Handling
 ************************************************************************/
int allocate_arrays() {
    ptc = malloc(np*sizeof(vector));
    vtx = malloc(nv*sizeof(vector));
    K = malloc(nv*sizeof(double));
    
    return 0;
}
int free_arrays() {
    free(ptc);
    free(vtx);
    free(K);
    
    return 0;
}
/************************************************************************
 *  File IO
 ************************************************************************/
int initFile(char *fname, int npts) {
    FILE *fp;
    
    fp = fopen(fname,"w");
    fprintf(fp,"%d\n%d\n%d\n",nv,np,npts);
    fclose(fp);
    
    return 0;
}
int frqwriteio(char *fname) {
    int i;
    FILE *fp;
    
    fp = fopen(fname,"a");
    
    for(i=0; i<nv; i++) {
        fprintf(fp,"%.17g\n%.17g\n",vtx[i].x,vtx[i].y);
    }
    
    fclose(fp);
    
    return 0;
}

/************************************************************************
 *  Log Messages
 ************************************************************************/
int printHeader() {
    
    int i;
    
    calcInvars();
    printf("Initial values of integrals of involution\n");
    printf("H: %.17g\tQ: %.17g\tP: %.17g\tL2: %.17g\n\n",H, Q, P, L2);
    
    for(i=0; i<np; i++) {
        printf("Initial particle position: \t%8.6lf, %8.6lf\n",ptc[i].x,ptc[i].y);
    }
    printf("\n");
    
    for(i=0; i<nv; i++) {
        printf("Initial vortex positions: \t%8.6lf, %8.6lf\n",vtx[i].x,vtx[i].y);
    }
    printf("\n");
    
    return 0;
}
int printFooter() {
    calcInvars();
    printf("Final values of integrals of involution\n");
    printf("H: %.17g\tQ: %.17g\tP: %.17g\tL2: %.17g\n\n",H, Q, P, L2);
    
    return 0;
}
