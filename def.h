#define PI (acos(-1.))
#define EPS (1e-8)
#define W 10.
#define SEP 2.

int nv, np;                         // number of vortices [#] and number of particles [#]

typedef struct vector {
    double x;
    double y;
}vector;

vector *vtx;                        // array of vortex (x,y) locations
vector *ptc;
double *K;                          // array of vortex vorticities
double dt;                          // global time step
double H, Q, P, L2;                 // integrals of involution
