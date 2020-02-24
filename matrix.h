
void PrintMatrix(int n, double* A);

int InvMatrix(int n, double* A, double* B);

void InitV(double *v, double *u);

void InitM(double *V, double *N);

void MxV ( double* V, double* v);

void MxM ( double* V, double* N );

void intxvec(double c, double* v, double* res);

void plus(double *v, double *u, double* res);

void minus(double *v, double *u, double* res);

void Init(double *v, double x, double y, double z);

void MatrixId ( double *N );

double inner(double *v, double* u, int n);

void cross(double* a, double* b, double *c);

void normalize(double *a, int n);

void reflect(double* a, double* b, double *c);
