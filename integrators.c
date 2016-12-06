#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "def.h"

double absolute(double x) {
    if( x < 0. ) return -x;
    return x;
}

/****************************************************************
 * Functions for standard RK4 Integrator
 ***************************************************************/
vector advFunc(vector pos) {
    int i;
    double Rij;
    vector u={0,0};
    
    for(i=0; i<nv; i++) {
        Rij = (pow(pos.x-vtx[i].x,2)+pow(pos.y-vtx[i].y,2));
        
        u.x -= K[i]*(pos.y-vtx[i].y)/Rij;
        u.y += K[i]*(pos.x-vtx[i].x)/Rij;
    }
    u.x /= (2.*PI);
    u.y /= (2.*PI);
    
    return u;
    
}
vector advVortex(vector pos, int idx) {
    int i;
    vector u={0,0};
    double Rij;
    
    for(i=0; i<nv; i++) {
        if( i == idx ) continue;
        Rij = (pow(pos.x-vtx[i].x,2)+pow(pos.y-vtx[i].y,2));
        
        u.x -= K[i]*(pos.y-vtx[i].y)/Rij;
        u.y += K[i]*(pos.x-vtx[i].x)/Rij;
    }
    u.x /= (4.*PI);
    u.y /= (4.*PI);
    
    return u;
}

vector rk4_pstep(vector pos0) {
    vector k1, k2, k3, k4;
    vector step;
    vector new_pos;
    
    k1 = advFunc(pos0);
    step = (vector){pos0.x+0.5*k1.x*dt,pos0.y+0.5*k1.y*dt};
    k2 = advFunc(step);
    step = (vector){pos0.x+0.5*k2.x*dt,pos0.y+0.5*k2.y*dt};
    k3 = advFunc(step);
    step = (vector){pos0.x+k3.x*dt,pos0.y+k3.y*dt};
    k4 = advFunc(step);
    
    new_pos = (vector){pos0.x+dt/6.*(k1.x+2*k2.x+2*k3.x+k4.x),pos0.y+dt/6.*(k1.y+2*k2.y+2*k3.y+k4.y)};
    
    return new_pos;
    
}
vector rk4_vstep(vector pos0, int idx) {
    vector k1, k2, k3, k4;
    vector step;
    vector new_pos;
    
    k1 = advVortex(pos0, idx);
    step = (vector){pos0.x+0.5*k1.x*dt,pos0.y+0.5*k1.y*dt};
    k2 = advVortex(step, idx);
    step = (vector){pos0.x+0.5*k2.x*dt,pos0.y+0.5*k2.y*dt};
    k3 = advVortex(step, idx);
    step = (vector){pos0.x+k3.x*dt,pos0.y+k3.y*dt};
    k4 = advVortex(step, idx);
    
    new_pos = (vector){pos0.x+dt/6.*(k1.x+2*k2.x+2*k3.x+k4.x),pos0.y+dt/6.*(k1.y+2*k2.y+2*k3.y+k4.y)};
    
    return new_pos;
    
}

/****************************************************************
 * Functions for symplectic integrator
 ***************************************************************/
int sia_1() {
	int i, j, it;
	double k0, k1;
	double *p, *q, *k0q, *k1q;
	double r2ij, ierr;

    p = malloc(nv*sizeof(double));
	q = malloc(nv*sizeof(double));
	k0q = malloc(nv*sizeof(double));
	k1q = malloc(nv*sizeof(double));

	/* Solve for p=K*y */
	for(i=0; i<nv; i++) {
		k0 = K[i]*vtx[i].y;
		for(j=0,k1=0; j<nv; j++) {
			if(j == i) continue;
			r2ij = pow(vtx[i].x-vtx[j].x,2)+pow(vtx[i].y-vtx[j].y,2);
			k1 += K[j]*(vtx[i].x-vtx[j].x)/r2ij;
		}
		k1 = k1*dt*K[i]/(4.*PI);
        
        

 		p[i] = k0 + k1;
	}

	/* iteratively solve for q=x */
	it = 0;
	for(i=0; i<nv; i++) q[i] = vtx[i].x;
    do {
        ierr = 0.;
        for(i=0; i<nv; i++) {
            k0q[i] = q[i];
            for(j=0,k1=0; j<nv; j++) {
                if(j == i) continue;
                r2ij = pow(k0q[i]-vtx[j].x,2)+pow(vtx[i].y-vtx[j].y,2);
                k1q[i] += K[j]*(vtx[i].y-vtx[j].y)/r2ij;
            }
            k1q[i] = -k1q[i]*dt/(4.*PI);
            q[i] = vtx[i].x+k1q[i];
            ierr += fabs(q[i] - k0q[i]);
        }
        it++;
    }while(it < 1e7 && ierr > 1e-16);

	for(i=0; i<nv; i++) {
		vtx[i].x = q[i];
        vtx[i].y = p[i]/K[i];
	}

	free(q);
    free(p);
    free(k0q);
    free(k1q);

	return 0;
}
int sia_2() {
    int i, j, it;
    double k0, k1, k2;
    double k2a, k2b, k2c, k2d;
    double *p, *q, *k0q, *k1q, *k2q;
    double r2ij, ierr;
    
    p = malloc(nv*sizeof(double));
    q = malloc(nv*sizeof(double));
    k0q = malloc(nv*sizeof(double));
    k1q = malloc(nv*sizeof(double));
    k1q = malloc(nv*sizeof(double));
    k2q = malloc(nv*sizeof(double));
    
    /* Solve for p=K*y */
    for(i=0; i<nv; i++) {
        k0 = K[i]*vtx[i].y;
        for(j=0,k1=0; j<nv; j++) {
            if(j == i) continue;
            r2ij = pow(vtx[i].x-vtx[j].x,2)+pow(vtx[i].y-vtx[j].y,2);
            k1 += K[j]*(vtx[i].x-vtx[j].x)/r2ij;
        }
        k1 = k1*dt*K[i]/(4.*PI);
        
        k2a=0.; k2b=0.; k2c=0.; k2d=0.;
        for(j=0; j<nv; j++) {
            if(j == i) continue;
            r2ij = (pow(vtx[i].x-vtx[j].x,2)+pow(vtx[i].y-vtx[j].y,2));
            k2a += K[i]*K[j]*(pow(vtx[i].y-vtx[j].y,2)-pow(vtx[i].x-vtx[j].x,2))/(r2ij*r2ij);
            k2b += K[j]*(vtx[i].y-vtx[j].y)/r2ij;
            k2c += K[j]*(vtx[i].x-vtx[j].x)/r2ij;
            k2d += K[j]*(r2ij-2*pow(vtx[i].y-vtx[j].y,2))/(r2ij*r2ij);
        }
        k2 = dt*dt/(32.*PI*PI)*(k2a*k2b+k2c*k2d);
        
        p[i] = k0 + k1 + k2;
    }
    
    /* iteratively solve for q=x */
    it = 0;
    for(i=0; i<nv; i++) q[i] = vtx[i].x;
    do {
        ierr = 0.;
        for(i=0; i<nv; i++) {
            k0q[i] = q[i];
            for(j=0,k1=0; j<nv; j++) {
                if(j == i) continue;
                r2ij = pow(k0q[i]-vtx[j].x,2)+pow(vtx[i].y-vtx[j].y,2);
                k1q[i] += K[j]*(vtx[i].y-vtx[j].y)/r2ij;
            }
            k1q[i] = -k1q[i]*dt/(4.*PI);
            k2a=0.; k2b=0.; k2c=0.; k2d=0.;
            for(j=0; j<nv; j++) {
                if(j == i) continue;
                r2ij = (pow(k0q[i]-vtx[j].x,2)+pow(vtx[i].y-vtx[j].y,2));
                k2a += K[j]*(-2)*(k0q[i]-vtx[j].x)*(vtx[i].y-vtx[j].y)/(r2ij*r2ij);
                k2b += K[j]*(vtx[i].y-vtx[j].y)/r2ij;
                k2c += K[j]*(k0q[i]-vtx[j].x)/r2ij;
                k2d += K[j]*(r2ij-2*pow(vtx[i].y-vtx[j].y,2))/(r2ij*r2ij);
            }
            k2q[i] = dt*dt/(32.*PI*PI)*(k2a*k2b+k2c*k2d);
            
            q[i] = vtx[i].x+k1q[i]+k2q[i];
            ierr += fabs(q[i] - k0q[i]);
            //printf("%lf\t", q[i]);
        }
        //printf("\n");
        it++;
    }while(it < 1e7 && ierr > 1e-16);
    
    for(i=0; i<nv; i++) {
        vtx[i].x = q[i];
        vtx[i].y = p[i]/K[i];
    }
    
    free(q);
    free(p);
    free(k0q);
    free(k1q);
        free(k2q);
    
    return 0;
}
vector sia_2old(vector pos, int idx) {

    int j, it;
    double k0, k1, k2;
    double k2a, k2b, k2c, k2d;
    double r2ij;
    double p, q, qt;
    
    k0 = pos.y;
    
    for(j=0, k1=0; j<nv; j++) {
        if(j == idx) continue;
        r2ij = (pow(pos.x-vtx[j].x,2)+pow(pos.y-vtx[j].y,2));
        k1 += K[j]*(pos.x-vtx[j].x)/r2ij;
    }
    k1 = k1*dt*K[idx]/(4.*PI);
    
    k2a=0.; k2b=0.; k2c=0.; k2d=0.;
    for(j=0; j<nv; j++) {
        if(j == idx) continue;
        r2ij = (pow(pos.x-vtx[j].x,2)+pow(pos.y-vtx[j].y,2));
        k2a += K[idx]*K[j]*(pow(pos.y-vtx[j].y,2)-pow(pos.x-vtx[j].x,2))/(r2ij*r2ij);
        k2b += K[j]*(pos.y-vtx[j].y)/r2ij;
        k2c += K[j]*(pos.x-vtx[j].x)/r2ij;
        k2d += K[j]*(r2ij-2*pow(pos.y-vtx[j].y,2))/(r2ij*r2ij);
    }
    k2 = dt*dt/(32.*PI*PI)*(k2a*k2b+k2c*k2d);
    p = k0 + k1 + k2;
    
    it = 0;
    q = pos.x;
    do {
        qt = q;
        
        for(j=0, k1=0; j<nv; j++) {
            if(j == idx) continue;
            r2ij = (pow(qt-vtx[j].x,2)+pow(pos.y-vtx[j].y,2));
            k1 += K[j]*(qt-vtx[j].y)/r2ij;
        }
        k1 = -k1*dt/(4.*PI);
        
        k2a=0.; k2b=0.; k2c=0.; k2d=0.;
        for(j=0; j<nv; j++) {
            if(j == idx) continue;
            r2ij = (pow(qt-vtx[j].x,2)+pow(pos.y-vtx[j].y,2));
            k2a += K[j]*(-2)*(qt-vtx[j].x)*(pos.y-vtx[j].y)/(r2ij*r2ij);
            k2b += K[j]*(pos.y-vtx[j].y)/r2ij;
            k2c += K[j]*(qt-vtx[j].x)/r2ij;
            k2d += K[j]*(r2ij-2*pow(pos.y-vtx[j].y,2))/(r2ij*r2ij);
        }
        k2 = dt*dt/(32.*PI*PI)*(k2a*k2b+k2c*k2d);
        q = k1+k2;
        
        it++;

    } while( fabs(q-qt) < 1e-10 );
    q = -q+pos.x;
    
    return (vector){q,p};

}
void calcLyapunov(char *fname1, char *fname2, int nt) {
    int i, it;
    double dtrack;
    FILE *fp1, *fp2, *fp3;
    double x1, y1, x2, y2;
    char line[50];
    
    fp1 = fopen(fname1,"r");
    fp2 = fopen(fname2,"r");
    fp3 = fopen("Lyapunov.txt","w");
	 fgets(line,sizeof(line), fp1);
	 fgets(line,sizeof(line), fp1);
	 fgets(line,sizeof(line), fp2);
	 fgets(line,sizeof(line), fp2);
    
    for(it=0; it<nt; it++) {
        for(i=0; i<nv; i++) {
            fgets(line, sizeof(line), fp1);
            sscanf(line,"%lf",&x1);
            fgets(line, sizeof(line), fp1);
            sscanf(line,"%lf",&y1);
            fgets(line, sizeof(line), fp2);
            sscanf(line,"%lf",&x2);
            fgets(line, sizeof(line), fp2);
            sscanf(line,"%lf",&y2);
            
            dtrack = sqrt(pow(x1-x2,2)+pow(y1-y2,2));
            if(i==7) fprintf(fp3,"%.17g\n",dtrack);
        }
    }
}
void calcInvars() {
    int i,j;
    double rij;
    
    H=0;
    for(i=0; i<nv; i++) {
        for(j=i+1; j<nv; j++) {
            if(j==i) continue;
            rij = sqrt(pow(vtx[j].x-vtx[i].x,2)+pow(vtx[j].y-vtx[i].y,2));
            H += K[i]*K[j]*log(rij);
        }
    }
    H /= -(4*PI);
    
    Q=0; P=0; L2=0;
    for(i=0; i<nv; i++) {
        Q  += K[i]*vtx[i].x;
        P  += K[i]*vtx[i].y;
        L2 += K[i]*(pow(vtx[i].x,2)+pow(vtx[i].y,2));
    }
    
}





