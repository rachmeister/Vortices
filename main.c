#include <stdio.h>
#include <stdlib.h>

#include "def.h"
#include "methods.h"

int main(int argc, char *argv[]) {
    
    initSRand();
    
    printf("VORTEX MOTION AND ADVECTED PARTICLE MODEL\n");
    printf("--------------------------------------------------------------------\n\n");
    
    /*printf("2 IDENTICAL VORTICES\n");
    printf("--------------------------------------------------------------------\n");
    init2Same();
    run("2same_sia.txt");
    
    printf("2 VORTICES WITH OPPOSITE VORTICITY\n");
    printf("--------------------------------------------------------------------\n");
    init2Opp();
    run("2opp_sia.txt");
    
    printf("2 VORTICES WITH DIFFERENT VORTICITY\n");
    printf("--------------------------------------------------------------------\n");
    init2Opp();
    run("2skew_sia.txt");
    
    printf("3 VORTICES ON ISOSCELES RIGHT TRIANGLE\n");
    printf("--------------------------------------------------------------------\n");
    init3right();
    run("3right_sia.txt");*/
    
    /*printf("4 VORTICES ON SQUARE\n");
    printf("--------------------------------------------------------------------\n");
    init4sq();
    run("4sq_rk4_unperturbed.txt");*/
    
    initCircle(8);
    run("unperturbed.txt");
    
    initCircle(8);
    perturb(7,1e-3);
    run("perturbed.txt");
    
    calcLyapunov("unperturbed.txt","perturbed.txt",5000);
	
    /*initRandom(4);
    run("path.txt");*/
    //init2Opp();
    //run("sia2.txt");
	
    return 0;
}

/************************************************************************
 * Run case
 ************************************************************************/
int run(char *fname) {
    int i;
    int frqw = 1000;
    double t, tottime;
    FILE *fp;
    
    vector *vtx_temp;
    
    /* Initilization */
    tottime = 1000.;
    dt = 1e-4;
    t = 0.;
    
    printHeader();
    
    vtx_temp = malloc(nv*sizeof(vector));
    /************************
     * Begin Time Loop
     ************************/
    //initFile(fname, (int)(tottime/dt/frqw));
    
    
    fp = fopen(fname,"w");
    fprintf(fp,"%d\n%d\n",nv,(int)(tottime/dt/frqw));
    while( t<=tottime ){
        /* Track Advected Particle */

        /* Advect Vortices */
        for(i=0; i<nv; i++) vtx_temp[i] = rk4_vstep(vtx[i],i);
        
        /* Output */
        if( (int)(t/dt)%frqw == 0 ) {
            //frqwriteio(fname);
            for(i=0; i<nv; i++) {
                fprintf(fp,"%.17g\n%.17g\n",vtx[i].x,vtx[i].y);
            }
        }
        
        //sia_1();
        
        /* Update */
        t += dt;
        //pos0 = pos;
        for(i=0; i<nv; i++) vtx[i] = vtx_temp[i];
        
    }
    fclose(fp);
    
    printFooter();
    free_arrays();
    
    return 0;
}

