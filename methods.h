vector advFunc(vector pos);
vector advVortex(vector pos, int idx);
vector rk4_pstep(vector pos0);
vector rk4_vstep(vector pos0, int idx);
int sia_1();
vector sia_2();
void calcLyapunov(char *fname1, char *fname2, int nt);
void calcInvars();

int run(char *fname);
int initSRand();
int init2Same();
int init2Opp();
int init2skew();
int init3right();
int init4sq();
int initCircle(int n);
int initRandom(int n);
int perturb(int idx, double delta);

int allocate_arrays();
int free_arrays();

int initFile(char *fname, int npts);
int frqwriteio(char *fname);
int printHeader();
int printFooter();
