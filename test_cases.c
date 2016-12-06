#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "def.h"
#include "methods.h"

int initSRand() {
    time_t stime;
    
    srand((unsigned)time(&stime));
    return 0;
}
/************************************************************************
 * 2 Vortices with same vorticity
 ************************************************************************/
int init2Same() {
    nv = 2;
    np = 0;
    
    allocate_arrays();
    
    K[0]=W;
    vtx[0].x=1.;
    vtx[0].y=1+SEP;
    K[1]=W;
    vtx[1].x=1.;
    vtx[1].y=1.+2*SEP;
    
    return 0;
}
/************************************************************************
 * 2 Vortices with opposite vorticity
 ************************************************************************/
int init2Opp() {
    nv = 2;
    np = 0;
    
    allocate_arrays();
    
    K[0]=W;
    vtx[0].x=0.;
    vtx[0].y=0.;
    K[1]=-W;
    vtx[1].x=0.;
    vtx[1].y=1.;
    
    return 0;
}
/************************************************************************
 * 2 Vortices 
 ************************************************************************/
int init2skew() {
    nv = 2;
    np = 0;
    
    allocate_arrays();
    
    K[0]=3;
    vtx[0].x=0.;
    vtx[0].y=0.;
    K[1]=2;
    vtx[1].x=1.;
    vtx[1].y=0.;
    
    return 0;
}
/************************************************************************
 * 3 Vortices on isoceles right triangle
 ************************************************************************/
int init3right() {
    nv = 3;
    np = 0;
    
    allocate_arrays();
    
    K[0]=3;
    vtx[0].x=1.;
    vtx[0].y=0.;
    K[1]=2;
    vtx[1].x=0.;
    vtx[1].y=1.;
    K[2]=1;
    vtx[2].x=-1.;
    vtx[2].y=0.;
    
    return 0;
}
/************************************************************************
 * 4 Vortices on square
 ************************************************************************/
int init4sq() {
    nv = 4;
    np = 0;
    
    allocate_arrays();
    
    K[0]=W;
    vtx[0].x=0.;
    vtx[0].y=0.;
    K[1]=W;
    vtx[1].x=1.;
    vtx[1].y=0.;
    K[2]=W;
    vtx[2].x=1.;
    vtx[2].y=1.;
    K[3]=W;
    vtx[3].x=0.;
    vtx[3].y=1.;
    
    return 0;
}
/************************************************************************
 * n Vortices equally spaced on circle
 ************************************************************************/
int initCircle(int n) {
    int i;
    double angle;
    
    nv = n;
    np = 0;
    allocate_arrays();
    
    angle = 2.*PI/(double)n;
    
    for(i=0; i<nv; i++) {
        K[i] = W;
        vtx[i].x = sin(i*angle)*4.;
        vtx[i].y = cos(i*angle)*4.;
    }
    return 0;
    
}
/************************************************************************
 * n Vortices randomly determined
 ************************************************************************/
int initRandom(int n) {
    
    int i;
    
    nv = n;
    np = 0;
    allocate_arrays();
    
    //pos0 = (vector){rand()%5,rand()%5};
    
    for(i=0; i<nv; i++) {
        K[i] = W*pow(-1,i);
        vtx[i].x = (rand()%20)/4.;
        vtx[i].y = (rand()%20)/4.;
        //printf("Vortex at:\t\t\t%lf, %lf\t Vorticity: %lf\n",vtx[i].x,vtx[i].y,K[i]);
    }
    
    return 0;
}
int perturb(int idx, double delta) {
    
    double angle;
    angle = (rand()%100)/100.*2*PI;
    
    vtx[idx].x += sin(angle)*delta;
    vtx[idx].y += cos(angle)*delta;
    
    return 0;
}
