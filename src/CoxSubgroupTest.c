#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>

#define QEPS 1e-6
#define GAMMA_M 5
#define MPI 3.1415926
#define MPI1 0.1591549   //1.0/(2*MPI)
#define MPI2 0.3989423   //1.0/sqrt(2*MPI)
#define SQRT2 0.7071068  //1.0/sqrt(2.0)
#define LOG2PI 1.837877  //log(2pi)
#define MEPS 1e-10
#define ITMAX 100
#define EPS 1.0e-8
#define FPMIN 1.0e-30
#define MIN(a, b) (((a)<(b))?(a):(b))
#define MAX(a, b) (((a)>(b))?(a):(b))
#define IDEX(a, b) (((a)<(b))?1:0)

void printArrayDouble(const double *arr, int n, int ncol){
    for (int i=0; i<n; i++){
        Rprintf("%f   ",arr[i]);
        if(!((i+1)%ncol)) Rprintf("\n");
    }
    printf("\n");
}

void printArrayDouble2(double **arr, int n, int ncol){
    for (int i=0; i<n; i++){
        for (int j = 0; j < ncol; j++ ) {
            Rprintf("%f   ",arr[i][j]);
        }
        Rprintf("\n");
    }
    printf("\n");
}

void printArrayDoubleInt(const int *arr, int n, int ncol){
    for (int i=0; i<n; i++){
        printf("%d   ",arr[i]);
        if(!((i+1)%ncol)) printf("\n");
    }
    printf("\n");
}

void printArrayDouble1(const double *arr, int n, int nrow, int ncol){
    if(nrow>n) fprintf(stderr, "nrow must not be greater than n!");
    for (int i=0; i<nrow; i++){
        for(int j=0; j<ncol; j++){
            printf("%f   ",arr[j*n+i]);
        }
        printf("\n");
    }
}

void AbyB(double *outVector, double *A, double *v, int n, int p, int q){
    int i,j,k;
    double tmp;
    for (i=0;i<n;i++){
        for(k=0;k<q;k++){
            tmp = 0;
            for(j=0;j<p;j++)
                tmp += A[j*n + i]*v[k*p + j];
            outVector[k*n+i] = tmp;
        }
    }
}

double incbeta(double x, double a, double b) {
    if (x < 0.0 || x > 1.0) return 1.0/0.0;

    /*The continued fraction converges nicely for x < (a+1)/(a+b+2)*/
    if (x > (a+1.0)/(a+b+2.0)) {
        return (1.0-incbeta(1.0-x,b,a)); /*Use the fact that beta is symmetrical.*/
    }

    /*Find the first part before the continued fraction.*/
    const double lbeta_ab = lgamma(a)+lgamma(b)-lgamma(a+b);
    const double front = exp(log(x)*a+log(1.0-x)*b-lbeta_ab) / a;

    /*Use Lentz's algorithm to evaluate the continued fraction.*/
    double f = 1.0, c = 1.0, d = 0.0;

    int i, m;
    for (i = 0; i <= ITMAX; ++i) {
        m = i/2;

        double numerator;
        if (i == 0) {
            numerator = 1.0; /*First numerator is 1.0.*/
        } else if (i % 2 == 0) {
            numerator = (m*(b-m)*x)/((a+2.0*m-1.0)*(a+2.0*m)); /*Even term.*/
        } else {
            numerator = -((a+m)*(a+b+m)*x)/((a+2.0*m)*(a+2.0*m+1)); /*Odd term.*/
        }

        /*Do an iteration of Lentz's algorithm.*/
        d = 1.0 + numerator * d;
        if (fabs(d) < FPMIN) d = FPMIN;
        d = 1.0 / d;

        c = 1.0 + numerator / c;
        if (fabs(c) < FPMIN) c = FPMIN;

        const double cd = c*d;
        f *= cd;

        /*Check for stop.*/
        if (fabs(1.0-cd) < EPS) {
            return front * (f-1.0);
        }
    }

    return 1.0/0.0; /*Needed more loops, did not converge.*/
}

void EstCoxphR(double *beta0, double *x, int *status, double *p_num,
                double *p_den, int n, int p, int maxstep, double eps, double tol){
    int i, j, k, step;
    double *d, *xb, *x_i;

    d       = (double*)malloc(sizeof(double)*p);
    xb      = (double*)malloc(sizeof(double)*n);
    double theta;
    double temp;
    tol /= eps;

    for (j = 0; j < p; j++) beta0[j] = 0.0;
    for (j = 0; j < n; j++) xb[j] = 0.0;

    for (step = 0; step < maxstep; ++step)
    {
        // calculate derivative
        for (i = 0; i < p; ++i) d[i] = 0.0;
        theta = 0.0;
        for (i = n-1; i >= 0; --i) {
            theta += exp(xb[i]);
            if (status[i]) {
                for (j = i; j < n; ++j) {
                    x_i = x;
                    temp = exp(xb[j])/theta;
                    for (k = 0; k < p; ++k) {
                        d[k] -= temp * (x_i[i]-x_i[j]);
                        x_i += n;
                    }

                }
            }
        }

        // update beta0
        for (j = 0; j < p; j++) beta0[j] -= d[j]*eps;
        // update xb = x %*% beta
        for (i = 0; i < n; ++i) {
            xb[i] = 0.0;
            x_i = x;
            for (j = 0; j < p; j++) {
                xb[i] += beta0[j] * x_i[i];
                x_i += n;
            }
        }

        // whether covergence
        temp = 0.0;
        for (i = 0; i < p; ++i) {
            theta = fabs(d[i]);
            if(theta > temp) temp = theta;
        }
        if (temp <= tol) break;

        if (step == maxstep-1) {
            Rprintf("Solution path unfinished, more iterations are needed.\n");
            break;
        }
    }

    p_num[n-1] = exp(xb[n-1]);
    p_den[n-1] = p_num[n-1];
    for (i = n-2; i >= 0; --i) {
        p_num[i] = exp(xb[i]);
        p_den[i] = p_den[i+1] + p_num[i];
    }

    free(d);
    free(xb);
}

SEXP _EST_COXPH(SEXP X_, SEXP STATUS, SEXP PARA_INT_, SEXP EPS0, SEXP TOL){
    // dimensions
    int *para   = INTEGER(PARA_INT_);
    int n       = para[0];
    int p       = para[1];
    int maxstep = para[2];

    // Outcome
    SEXP rBeta, p_num, p_den, list, list_names;
    PROTECT(rBeta       = allocVector(REALSXP,  p));
    PROTECT(p_num       = allocVector(REALSXP,  n));
    PROTECT(p_den       = allocVector(REALSXP,  n));
    PROTECT(list_names  = allocVector(STRSXP,   3));
    PROTECT(list        = allocVector(VECSXP,   3));

    EstCoxphR(REAL(rBeta), REAL(X_), INTEGER(STATUS), REAL(p_num), REAL(p_den), n, p, maxstep, REAL(EPS0)[0], REAL(TOL)[0]);

    SET_STRING_ELT(list_names,  0,  mkChar("coef"));
    SET_STRING_ELT(list_names,  1,  mkChar("p_num"));
    SET_STRING_ELT(list_names,  2,  mkChar("p_den"));
    SET_VECTOR_ELT(list,        0,  rBeta);
    SET_VECTOR_ELT(list,        1, p_num);
    SET_VECTOR_ELT(list,        2, p_den);
    setAttrib(list, R_NamesSymbol,  list_names);

    UNPROTECT(5);
    return list;
}

//-------------- Cox subgroup test ---------------------------------------------
double** vec2matrix(double *vec, int n, int p) {
    // vec is a array with n*p elements, by row saved as a matrix
    double **M;
    M = (double**)malloc(n*sizeof(double*));
    for ( int i = 0; i < n; i++ ) M[i] = vec + i*p;
    return M;
}

int** ivec2matrix(int *vec, int n, int p) {
    // vec is a array with n*p elements, by row saved as a matrix
    int **M;
    M = (int**)malloc(n*sizeof(int*));
    for ( int i = 0; i < n; i++ ) M[i] = vec + i*p;
    return M;
}

void get_omega(double *U, double *omega, int n, int pu) {
    // U is a n by pu matrix, saved by row
    // omega is a n by n matrix
    int i, j, k;
    double *ui, *uj, tmp, *oi, *normU;

    normU = (double*)malloc(sizeof(double)*n);

    // calculate norm U
    ui = U;
    for ( i = 0; i < n; i++) {
        tmp = 0.0;
        for ( j = 0; j < pu; j++) {
            tmp += ui[j] * ui[j];
        }
        normU[i] = 1.0 / sqrt( tmp );
        ui += pu;
    }

    // calculate omwga_ij
    oi = omega;
    ui = U;
    for ( i = 0; i < n; i++) {
        uj = ui + pu;
        for ( j = i+1; j < n; j++) {
            tmp = 0.0;
            for ( k = 0; k < pu; k++) {
                tmp += ui[k] * uj[k];
            }
            tmp *= normU[i] * normU[j];
            tmp = tmp / sqrt( 1.0 - tmp*tmp );
            oi[j] = atan( tmp ) * MPI1 + 0.25;
            uj += pu;
        }
        ui += pu;
        oi += n;
    }
    
    free(normU);
}

void get_zTz(double *Z, double *ztz, int n, int pz) {
    // Z is a n by pz matrix, saved by row
    // ztz is a n by n matrix
    int i, j, k;
    double *zi, *zj, *ztzi, tmp;

    for(i = 0; i < n; i++){
        zi = Z + pz*i;
        tmp = 0.0;
        for(k = 0; k < pz; k++){
            tmp += zi[k]*zi[k];
        }
        ztzi = ztz + (n+1)*i;
        ztzi[0] = tmp/2.0;
    }

    zi = Z + pz;
    ztzi = ztz + n;
    for ( i = 1; i < n; i++) {
        zj = Z;
        for ( j = 0; j < i; j++) {
            ztzi[j] = 0.0;
            for ( k = 0; k < pz; k++) {
                ztzi[j] += zi[k] * zj[k];
            }
            zj += pz;
        }
        zi += pz;
        ztzi += n;
    }
}

void get_expxa(double *X, double *alpha, double *expxa, int n, int px) {
    // X is a n by px matrix, saved by row
    // alpha is a px-vector
    int i, j;
    double *xi;

    xi = X;
    for ( i = 0; i < n; i++) {
        expxa[i] = 0.0;
        for ( j = 0; j < px; j++) {
            expxa[i] += xi[j] * alpha[j];
        }
        expxa[i] = exp( expxa[i] );
        xi += px;
    }
}

double coxsubgrouptest(int *yrmax, int *yrmin, double *X, double *Z, double *U, int *status, double *alpha, int *param, double *omega) {
    int n  = param[0];
    int px = param[1];
    int pz = param[2];
    int pu = param[3];
    int i,j;
    double ts;

    double *c, *xi, *zeta, *theta, **OmegaM;
    c      = (double*)malloc(sizeof(double)*n);
    xi     = (double*)malloc(sizeof(double)*n);
    zeta   = (double*)malloc(sizeof(double)*n);
    theta  = (double*)malloc(sizeof(double)*n);
    OmegaM = (double**)malloc(sizeof(double*)*n);

    get_omega(U, omega, n, pu);
    get_zTz(Z, omega, n, pz);
    get_expxa(X, alpha, theta, n, px);
    for ( i = 0; i < n; i++) {
        OmegaM[i] = omega + i*n;
        // printf("omega = %lf\t",OmegaM[i][i]);
        // OmegaM[i][i] = 0.0; // not necessary, remove it will not affect result
    }

    // calculate c, xi and zeta
    c[n-1] = theta[n-1] / n;
    for ( i = n-2; i >= 0; i--) {
        c[i] = c[i+1] + theta[i] / n;
    }
    for ( i = n-1; i >= 0; i--) {
        c[i] = 1.0 / c[yrmin[i]];
    }
    xi[0] = 0.0;
    zeta[0] = 0.0;
    for ( i = 0; i < n-1; i++) {
        xi[i+1] = xi[i] + status[i] * c[i];
        zeta[i+1] = zeta[i] + status[i] * c[i] * c[i];
    }
    for ( i = 0; i < n; i++) {
        if (i == yrmax[i]) continue;
        xi[i]   = xi[yrmax[i]];
        zeta[i] = zeta[yrmax[i]];
    }
    
    // calculate test statistic
    ts = 0.0;
    double dum2 = 1.0 / (n);
    double dum3 = dum2 / (n);
    double tmp;
    double yij;
    for ( i = 0; i < n; i++) {
        for ( j = i+1; j < n; j++) {
            yij = (yrmin[i] == yrmin[j])*status[j]*c[j];
            tmp = (xi[i]+status[i]*c[i])*(xi[j]+status[j]*c[j])*theta[i]*theta[j]*dum3;
            // tmp = ( (xi[i] - yij)*(xi[j]-status[i]*c[i]) - zeta[i] + yij*c[j] )*theta[i]*theta[j]*dum3;
            if (status[i]) {
                tmp += status[j] - (xi[j]+status[j]*c[j])*theta[j]*dum2;
                // tmp += status[j] - (xi[j]-c[i])*theta[j]*dum2;
            }
            if (status[j]) {
                tmp -= (xi[i]+status[i]*c[i])*theta[i]*dum2;
                // tmp -= (xi[i] - yij)*theta[i]*dum2;
            }
            ts += OmegaM[i][j] * OmegaM[j][i] * tmp;
        }
        ts += OmegaM[i][i]*theta[i]*theta[i]*((xi[i]+status[i]*c[i])*(xi[i]+status[i]*c[i]))*dum3/(2.0)-status[i]*OmegaM[i][i]*theta[i]*(xi[i]+status[i]*c[i])*dum2+OmegaM[i][i]*status[i]/(2.0);
    }
    ts *= 2.0/(n*(n));

    // for(i =0; i < n; i++){
    //     ts += OmegaM[i][i]*theta[i]*theta[i]*((xi[i]+status[i]*c[i])*(xi[i]+status[i]*c[i])-zeta[i]-status[i]*c[i]*c[i])/(2.0*dum3)-status[i]*OmegaM[i][i]*theta[i]*xi[i]*dum2
    // }
    // ts *= 2.0/(n*(n-1.0));

    free(theta);
    free(OmegaM);
    free(xi);
    free(zeta);
    free(c);

    return ts;
}

SEXP _COX_Subgroup_Test(SEXP YRMAX_, SEXP YRMIN_, SEXP X_, SEXP Z_,SEXP U_, SEXP STATUS_, SEXP ALPHA_, SEXP PARA_INT_) {
    // Y is ascending
    int n = INTEGER(PARA_INT_)[0];

    SEXP ts, omega, list, list_names;
    PROTECT(ts          = allocVector(REALSXP,  1));
    PROTECT(omega       = allocVector(REALSXP,  n*n));
    PROTECT(list_names  = allocVector(STRSXP,   2));
    PROTECT(list        = allocVector(VECSXP,   2));

    REAL(ts)[0] = coxsubgrouptest(
        INTEGER(YRMAX_), INTEGER(YRMIN_), REAL(X_), REAL(Z_), REAL(U_), 
        INTEGER(STATUS_), REAL(ALPHA_), INTEGER(PARA_INT_), REAL(omega)
    );

    SET_STRING_ELT(list_names,  0,  mkChar("test_statistic"));
    SET_STRING_ELT(list_names,  1,  mkChar("omega"));
    SET_VECTOR_ELT(list,        0,  ts);
    SET_VECTOR_ELT(list,        1,  omega);
    setAttrib(list, R_NamesSymbol,  list_names);

    UNPROTECT(4);
    return list;
}

double coxsubgrouptest_given_omega(int *yrmax, int *yrmin, double *X, int *status, double *alpha, int *param, double *omega) {
    int n  = param[0];
    int px = param[1];
    int i, j;
    double ts; 

    double *theta, **OmegaM, *c, *xi, *zeta;
    // upper triangle saves omega, lower triangle saves ztz
    // omega  = (double*)malloc(sizeof(double)*n*n);
    theta  = (double*)malloc(sizeof(double)*n);
    OmegaM = (double**)malloc(sizeof(double*)*n);
    xi     = (double*)malloc(sizeof(double)*n);
    zeta   = (double*)malloc(sizeof(double)*n);
    c      = (double*)malloc(sizeof(double)*n);
    get_expxa(X, alpha, theta, n, px);
    for ( i = 0; i < n; i++) {
        OmegaM[i] = omega + i*n;
        OmegaM[i][i] = 0.0; // not necessary, remove it will not affect result
    }

    // calculate c, xi and zeta
    c[n-1] = theta[n-1] / n;
    for ( i = n-2; i >= 0; i--) {
        c[i] = c[i+1] + theta[i] / n;
    }
    for ( i = n-1; i >= 0; i--) {
        c[i] = 1.0 / c[yrmin[i]];
    }
    xi[0] = 0.0;
    zeta[0] = 0.0;
    for ( i = 0; i < n-1; i++) {
        xi[i+1] = xi[i] + status[i] * c[i];
        zeta[i+1] = zeta[i] + status[i] * c[i] * c[i];
    }
    for ( i = 0; i < n; i++) {
        if (i == yrmax[i]) continue;
        xi[i]   = xi[yrmax[i]];
        zeta[i] = zeta[yrmax[i]];
    }
    
    // calculate test statistic
    ts = 0.0;
    double dum2 = 1.0 / (n-2.0);
    double dum3 = dum2 / (n-3.0);
    double tmp;
    double yij;
    for ( i = 0; i < n; i++) {
        for ( j = i+1; j < n; j++) {
            yij = (yrmin[i] == yrmin[j])*status[j]*c[j];
            tmp = ( (xi[i] - yij)*(xi[j]-status[i]*c[i]) - zeta[i] + yij*c[j] )*theta[i]*theta[j]*dum3;
            if (status[i]) {
                tmp += status[j] - (xi[j]-c[i])*theta[j]*dum2;
            }
            if (status[j]) {
                tmp -= (xi[i] - yij)*theta[i]*dum2;
            }
            ts += OmegaM[i][j] * OmegaM[j][i] * tmp;
        }
    }
    ts *= 2.0/(n*(n-1.0));

    // free(omega);
    free(theta);
    free(OmegaM);
    free(xi);
    free(zeta);
    free(c);

    return ts;
}

SEXP _COX_Subgroup_Test_given_omega(SEXP YRMAX_, SEXP YRMIN_, SEXP X_, SEXP Omega_, SEXP STATUS_, SEXP ALPHA_, SEXP PARA_INT_) {
    // Y is ascending
    SEXP ts, list, list_names;
    PROTECT(ts          = allocVector(REALSXP,  1));
    PROTECT(list_names  = allocVector(STRSXP,   1));
    PROTECT(list        = allocVector(VECSXP,   1));

    REAL(ts)[0] = coxsubgrouptest_given_omega(
        INTEGER(YRMAX_), INTEGER(YRMIN_), REAL(X_), INTEGER(STATUS_), 
        REAL(ALPHA_), INTEGER(PARA_INT_), REAL(Omega_)
    );

    SET_STRING_ELT(list_names,  0,  mkChar("test_statistic"));
    SET_VECTOR_ELT(list,        0,  ts);
    setAttrib(list, R_NamesSymbol,  list_names);

    UNPROTECT(3);
    return list;
}

void coxsubgrouptest_permutation(
    int *yrmax, int *yrmin, double *X, double *Z, double *U, int *status, 
    double *alpha, int *param, // double *omega, 
    int **idM, double *tsboot, double *pval, double *ts
    ) {
    int n  = param[0];
    int px = param[1];
    int pz = param[2];
    int pu = param[3];
    int b  = param[4];
    int i, j, bb;

    double *omega;
    omega = (double*)malloc(sizeof(double)*n*n);

    double *theta, **OmegaM, *c, *xi, *zeta;
    // upper triangle saves omega, lower triangle saves ztz
    // omega  = (double*)malloc(sizeof(double)*n*n);
    theta  = (double*)malloc(sizeof(double)*n);
    OmegaM = (double**)malloc(sizeof(double*)*n);
    xi     = (double*)malloc(sizeof(double)*n);
    zeta   = (double*)malloc(sizeof(double)*n);
    c      = (double*)malloc(sizeof(double)*n);
    get_omega(U, omega, n, pu);
    get_zTz(Z, omega, n, pz);
    get_expxa(X, alpha, theta, n, px);
    for ( i = 0; i < n; i++) {
        OmegaM[i] = omega + i*n;
        OmegaM[i][i] = 0.0; // not necessary, remove it will not affect result
    }

    // calculate c, xi and zeta
    c[n-1] = theta[n-1] / n;
    for ( i = n-2; i >= 0; i--) {
        c[i] = c[i+1] + theta[i] / n;
    }
    for ( i = n-1; i >= 0; i--) {
        c[i] = 1.0 / c[yrmin[i]];
    }
    // printArrayDouble(c, n, n);
    xi[0] = 0.0;
    zeta[0] = 0.0;
    for ( i = 0; i < n-1; i++) {
        xi[i+1] = xi[i] + status[i] * c[i];
        zeta[i+1] = zeta[i] + status[i] * c[i] * c[i];
    }
    for ( i = 0; i < n; i++) {
        if (i == yrmax[i]) continue;
        xi[i]   = xi[yrmax[i]];
        zeta[i] = zeta[yrmax[i]];
    }
    // printArrayDouble(xi, n, n);
    // printArrayDouble(zeta, n, n);
    
    // calculate test statistic
    ts[0] = 0.0;
    double dum2 = 1.0 / (n-2.0);
    double dum3 = dum2 / (n-3.0);
    double tmp;
    double yij;
    for ( bb = 0; bb < b; bb++ ) tsboot[bb]  = 0.0;
    for ( i = 0; i < n; i++) {
        for ( j = i+1; j < n; j++) {
            yij = (yrmin[i] == yrmin[j])*status[j]*c[j];
            tmp = ( (xi[i] - yij)*(xi[j]-status[i]*c[i]) - zeta[i] + yij*c[j] )*theta[i]*theta[j]*dum3;
            if (status[i]) {
                tmp += status[j] - (xi[j]-c[i])*theta[j]*dum2;
            }
            if (status[j]) {
                tmp -= (xi[i] - yij)*theta[i]*dum2;
            }
            ts[0] += OmegaM[i][j] * OmegaM[j][i] * tmp;
            for ( bb = 0; bb < b; bb++ ){
                if (idM[i][bb] < idM[j][bb]) {
                    tsboot[bb] += OmegaM[idM[i][bb]][idM[j][bb]] * OmegaM[j][i] * tmp;
                } else {
                    tsboot[bb] += OmegaM[idM[j][bb]][idM[i][bb]] * OmegaM[j][i] * tmp;
                }
            }
        }
    }
    ts[0] *= 2.0/(n*(n-1.0));
    for ( bb = 0; bb < b; bb++ ) tsboot[bb] *= 2.0/(n*(n-1.0));

    // calculate pval
    pval[0] = 0.0;
    for ( bb = 0; bb < b; bb++ ) {
        if (tsboot[bb] > ts[0]) pval[0] += 1.0;
    }
    pval[0] /= b;

    free(omega);
    free(theta);
    free(OmegaM);
    free(xi);
    free(zeta);
    free(c);
}

SEXP _COX_Subgroup_Test_Permutation(
    SEXP YRMAX_, SEXP YRMIN_, SEXP X_, SEXP Z_,SEXP U_, SEXP STATUS_, SEXP ALPHA_, 
    SEXP PARA_INT_, SEXP SAMPLE_) {
    // Y is ascending
    int n = INTEGER(PARA_INT_)[0];
    int b = INTEGER(PARA_INT_)[4];

    SEXP ts, tsboot, pval, list, list_names;
    PROTECT(ts          = allocVector(REALSXP,  1));
    PROTECT(tsboot      = allocVector(REALSXP,  b));
    PROTECT(pval        = allocVector(REALSXP,  1));
    // PROTECT(omega       = allocVector(REALSXP,  n*n));
    PROTECT(list_names  = allocVector(STRSXP,   3));
    PROTECT(list        = allocVector(VECSXP,   3));

    int **sample = ivec2matrix(INTEGER(SAMPLE_), n, b);

    coxsubgrouptest_permutation(
        INTEGER(YRMAX_), INTEGER(YRMIN_), REAL(X_), REAL(Z_), REAL(U_), 
        INTEGER(STATUS_), REAL(ALPHA_), INTEGER(PARA_INT_), // REAL(omega), 
        sample, REAL(tsboot), REAL(pval), REAL(ts)
    );

    SET_STRING_ELT(list_names,  0,  mkChar("test_statistic"));
    SET_STRING_ELT(list_names,  1,  mkChar("tsboot"));
    SET_STRING_ELT(list_names,  2,  mkChar("Pvalue"));
    // SET_STRING_ELT(list_names,  3,  mkChar("omega"));
    SET_VECTOR_ELT(list,        0,  ts);
    SET_VECTOR_ELT(list,        1,  tsboot);
    SET_VECTOR_ELT(list,        2,  pval);
    // SET_VECTOR_ELT(list,        3,  omega);
    setAttrib(list, R_NamesSymbol,  list_names);

    free(sample);
    UNPROTECT(5);
    return list;
}

// for benchmarks --------------------------------------------------------------
void calculate_v4score(double *exb, int *status, int *yrmin, int *yrmax, double *eta, double *vec, int n) {
    int i, j;
    double tmp;

    tmp = 0.0;
    j = n-1;
    for ( i = n-1; i >= 0; i--) {
        while (j >= yrmin[i]) {
            tmp += exb[j];
            j--;
        }
        eta[i] = tmp;
    }

    tmp = 0.0;
    j = 0;
    for ( i = 0; i < n; i++) {
        while (j <= yrmax[i]) {
            tmp += status[j] / eta[j];
            j++;
        }
        vec[i] = status[i] - tmp*exb[i];
    }
}

SEXP _calculate_v4score(SEXP Exb_, SEXP Status_, SEXP YRMIN_, SEXP YRMAX_, SEXP N) {
    int n = INTEGER(N)[0];

    double *eta;
    eta  = (double*)malloc(sizeof(double)*n);

    SEXP Vec;
    PROTECT(Vec = allocVector(REALSXP,  n));

    calculate_v4score(
        REAL(Exb_), INTEGER(Status_), INTEGER(YRMIN_), INTEGER(YRMAX_), eta, REAL(Vec), n
    );

    free(eta);
    UNPROTECT(1);
    return Vec;
}

void cholesky_decomposition_(double **M, int n)
{
    // M is a n by n symmetry matrix, upper triangle matrix saved is enough
    int i, j, k;
    double tmp;

    for ( i = 0; i < n; i++ ) {
        tmp = M[i][i];
        for ( j = 0; j < i; j++ ) tmp -= M[j][i] * M[j][i];
        tmp = sqrt(fabs(tmp));
        for ( j = 0; j < i; j++ ) {
            for ( k = i; k < n; k++ ) M[i][k] -= M[j][k] * M[j][i];
        }
        if ( tmp > QEPS ) {
            for ( k = i; k < n; k++ ) M[i][k] /= tmp;
        }
    }
}

void calculate_Dsqrt(double *exb, int *status, double *eta, int *yrmax, double *Dvec, int n) {
    int i, j, k;

    double **D = vec2matrix(Dvec, n, n);
    for ( i = 0; i < n*n; i++) Dvec[i] = 0.0;
    
    double dtmp1 = 0.0;
    double dtmp2 = 0.0;
    double dtmp3;
    j = 0;
    for ( i = 0; i < n; i++ ) {
        while (j <= yrmax[i]) {
            if (status[i]) {
                // dtmp3 = status[j] / eta[j];
                dtmp3 = 1.0 / eta[j];
                dtmp1 += dtmp3;
                dtmp2 += dtmp3*dtmp3;
            }
            
            j++;
        }

        // update D
        dtmp3 = -exb[i]*dtmp2;
        D[i][i] = exb[i]*(dtmp1 + dtmp3);
        for ( k = i+1; k < n; k++ ) {
            D[i][k] = exb[k]*dtmp3;
        }
    }
    cholesky_decomposition_(D, n);

    free(D);
}

SEXP _calculate_v4score_and_Dsqrt(SEXP Exb_, SEXP Status_, SEXP YRMIN_, SEXP YRMAX_, SEXP N) {
    int n = INTEGER(N)[0];

    double *eta;
    eta  = (double*)malloc(sizeof(double)*n);

    SEXP Vec, Dvec, list, list_names;
    PROTECT(Vec = allocVector(REALSXP,  n));
    PROTECT(Dvec = allocVector(REALSXP,  n*n));
    PROTECT(list_names  = allocVector(STRSXP,   2));
    PROTECT(list        = allocVector(VECSXP,   2));
    

    calculate_v4score(
        REAL(Exb_), INTEGER(Status_), INTEGER(YRMIN_), INTEGER(YRMAX_), eta, REAL(Vec), n
    );
    calculate_Dsqrt(REAL(Exb_), INTEGER(Status_), eta, INTEGER(YRMAX_), REAL(Dvec), n);

    SET_STRING_ELT(list_names,  0,  mkChar("v4score"));
    SET_STRING_ELT(list_names,  1,  mkChar("Dsqrt"));
    SET_VECTOR_ELT(list,        0,  Vec);
    SET_VECTOR_ELT(list,        1,  Dvec);
    setAttrib(list, R_NamesSymbol,  list_names);

    free(eta);
    UNPROTECT(4);
    return list;
}

int MatrixInvSymmetric(double **a,int n){
	int i,j,k,m;
    double w,g,*b;
    b = (double*)malloc(n*sizeof(double));

    for (k=0; k<=n-1; k++){
        w=a[0][0];
        if (fabs(w)+1.0==1.0){
            free(b); return(-1);
        }
        m=n-k-1;
        for (i=1; i<=n-1; i++){
            g=a[i][0]; b[i]=g/w;
            if (i<=m) b[i]=-b[i];
            for (j=1; j<=i; j++)
                a[i-1][j-1]=a[i][j]+g*b[j];
        }
        a[n-1][n-1]=1.0/w;
        for (i=1; i<=n-1; i++)
            a[n-1][i-1]=b[i];
    }
    for (i=0; i<=n-2; i++)
        for (j=i+1; j<=n-1; j++)
            a[i][j]=a[j][i];
    free(b);
    return(0);
}

void calculate_Hessian_inv_only(
    double **X, double *exb, int *status, double *eta, int *yrmin, int *yrmax, 
    double *ejXj, double **Hessian, int n, int p, int px, double **HessianH0
    ) {
    int i, j, i1, i2;
    double tmp;

    for ( i1 = 0; i1 < p; i1++ ) {
        Hessian[i1][i1] = 0.0;
        for ( i2 = i1+1; i2 < p; i2++ ) {
            Hessian[i1][i2] = 0.0;
        }
    }

    tmp = 0.0;
    j = n-1;
    for ( i = n-1; i >= 0; i--) {
        while (j >= yrmin[i]) {
            tmp += exb[j];
            j--;
        }
        eta[i] = tmp;
    }

    double dtmp1 = 0.0;
    j = 0;
    for ( i = 0; i < n; i++ ) {
        while (j <= yrmax[i]) {
            if (status[i]) {
                dtmp1 += 1.0 / eta[j];
            }
            j++;
        }

        tmp = exb[i] * dtmp1;
        for ( i1 = 0; i1 < p; i1++ ) {
            Hessian[i1][i1] += tmp * X[i][i1] * X[i][i1];
            for ( i2 = i1+1; i2 < p; i2++ ) {
                Hessian[i1][i2] += tmp * X[i][i1] * X[i][i2];
            }
        }
    }

    for ( i1 = 0; i1 < p; i1++ ) ejXj[i1] = 0.0;
    j = n-1;
    for ( i = n-1; i >= 0; i-- ) {
        while (j >= yrmin[i]) {
            for ( i1 = 0; i1 < p; i1++ ) {
                ejXj[i1] += X[j][i1] * exb[j];
            }
            j--;
        }

        if (status[i] == 0) continue;

        tmp = 1.0 / ( eta[i] * eta[i] );
        for ( i1 = 0; i1 < p; i1++ ) {
            Hessian[i1][i1] -= tmp * ejXj[i1] * ejXj[i1];
            for ( i2 = i1+1; i2 < p; i2++ ) {
                Hessian[i1][i2] -= tmp * ejXj[i1] * ejXj[i2];
            }
        }
    }

    for ( i1 = 0; i1 < p; i1++ ) {
        for ( i2 = i1+1; i2 < p; i2++ ) {
            Hessian[i2][i1] = Hessian[i1][i2];
        }
    }
    MatrixInvSymmetric(Hessian, p);
    for ( i1 = 0; i1 < px; i1++ ) {
        Hessian[i1][i1] -= HessianH0[i1][i1];
        for ( i2 = i1+1; i2 < px; i2++ ) {
            Hessian[i2][i1] -= HessianH0[i2][i1];
            Hessian[i1][i2] = Hessian[i2][i1];
        }
    }
}

SEXP _calculate_Hessian_inv_only(
    SEXP X_, SEXP Exb_, SEXP Status_, SEXP YRMIN_, SEXP YRMAX_, SEXP N, SEXP P,
    SEXP PX, SEXP Hessian_inv_H0) {
    int n  = INTEGER(N)[0];
    int p  = INTEGER(P)[0];
    int px = INTEGER(PX)[0];

    double *ejXj, *eta;
    ejXj  = (double*)malloc(sizeof(double)*p);
    eta  = (double*)malloc(sizeof(double)*n);

    SEXP Hessian_inv;
    PROTECT(Hessian_inv = allocVector(REALSXP, p*p));
    double **Hessian = vec2matrix(REAL(Hessian_inv), p, p);
    double **X = vec2matrix(REAL(X_), n, p);
    double **HessianH0 = vec2matrix(REAL(Hessian_inv_H0), px, px);

    calculate_Hessian_inv_only(
        X, REAL(Exb_), INTEGER(Status_), eta, INTEGER(YRMIN_), INTEGER(YRMAX_), 
        ejXj, Hessian, n, p, px, HessianH0
    );

    free(eta);
    free(ejXj);
    free(Hessian);
    free(HessianH0);
    free(X);
    UNPROTECT(1);
    return Hessian_inv;
}

void calculate_Hessian_inv(
    double **X, double *exb, int *status, double *eta, int *yrmin, int *yrmax, 
    double *ejXj, double **Hessian, int n, int p
    ) {
    int i, j, i1, i2;
    double tmp;

    for ( i1 = 0; i1 < p; i1++ ) {
        Hessian[i1][i1] = 0.0;
        for ( i2 = i1+1; i2 < p; i2++ ) {
            Hessian[i1][i2] = 0.0;
        }
    }

    double dtmp1 = 0.0;
    j = 0;
    for ( i = 0; i < n; i++ ) {
        while (j <= yrmax[i]) {
            if (status[i]) {
                dtmp1 += 1.0 / eta[j];
            }
            j++;
        }

        tmp = exb[i] * dtmp1;
        for ( i1 = 0; i1 < p; i1++ ) {
            Hessian[i1][i1] += tmp * X[i][i1] * X[i][i1];
            for ( i2 = i1+1; i2 < p; i2++ ) {
                Hessian[i1][i2] += tmp * X[i][i1] * X[i][i2];
            }
        }
    }

    for ( i1 = 0; i1 < p; i1++ ) ejXj[i1] = 0.0;
    j = n-1;
    for ( i = n-1; i >= 0; i-- ) {
        while (j >= yrmin[i]) {
            for ( i1 = 0; i1 < p; i1++ ) {
                ejXj[i1] += X[j][i1] * exb[j];
            }
            j--;
        }

        if (status[i] == 0) continue;

        tmp = 1.0 / ( eta[i] * eta[i] );
        for ( i1 = 0; i1 < p; i1++ ) {
            Hessian[i1][i1] -= tmp * ejXj[i1] * ejXj[i1];
            for ( i2 = i1+1; i2 < p; i2++ ) {
                Hessian[i1][i2] -= tmp * ejXj[i1] * ejXj[i2];
            }
        }
    }

    for ( i1 = 0; i1 < p; i1++ ) {
        for ( i2 = i1+1; i2 < p; i2++ ) {
            Hessian[i2][i1] = Hessian[i1][i2];
        }
    }
    MatrixInvSymmetric(Hessian, p);
}

void calculate_score(
    double **X, double *exb, int *status, int *yrmin, int *yrmax, double *eta, 
    double *score, int n, int p) {
    int i, j, k;
    double tmp;
    double vec;

    for ( i = 0; i < p; i++) score[i] = 0.0;

    tmp = 0.0;
    j = n-1;
    for ( i = n-1; i >= 0; i--) {
        while (j >= yrmin[i]) {
            tmp += exb[j];
            j--;
        }
        eta[i] = tmp;
    }

    tmp = 0.0;
    j = 0;
    for ( i = 0; i < n; i++) {
        while (j <= yrmax[i]) {
            tmp += status[j] / eta[j];
            j++;
        }
        vec = status[i] - tmp*exb[i];
        for ( k = 0; k < p; k++) score[k] += vec * X[i][k];
    }
}

void calculate_test(double *score, double **Hessian, double *test, int p) {
    int i, j;

    test[0] = 0.0;
    for ( i = 0; i < p; i++ ) {
        test[0] += Hessian[i][i]*score[i]*score[i];
        for ( j = i+1; j < p; j++ ) {
            test[0] += 2.0*Hessian[i][j]*score[i]*score[j];
        }
    }
}

SEXP _calculate_score_and_Hessianinv(
    SEXP X_, SEXP Exb_, SEXP Status_, SEXP YRMIN_, SEXP YRMAX_, SEXP N, SEXP P) {
    int n  = INTEGER(N)[0];
    int p  = INTEGER(P)[0];

    double *ejXj, *eta;
    ejXj  = (double*)malloc(sizeof(double)*p);
    eta  = (double*)malloc(sizeof(double)*n);

    SEXP ScoreR, HessianR, TestR, list, list_names;
    PROTECT(ScoreR      = allocVector(REALSXP,  p));
    PROTECT(HessianR    = allocVector(REALSXP, p*p));
    PROTECT(TestR       = allocVector(REALSXP,  1));
    PROTECT(list_names  = allocVector(STRSXP,   3));
    PROTECT(list        = allocVector(VECSXP,   3));
    double **Hessian = vec2matrix(REAL(HessianR), p, p);
    double **X = vec2matrix(REAL(X_), n, p);

    calculate_score(
        X, REAL(Exb_), INTEGER(Status_), INTEGER(YRMIN_), INTEGER(YRMAX_), eta, REAL(ScoreR), n, p
    );

    calculate_Hessian_inv(
        X, REAL(Exb_), INTEGER(Status_), eta, INTEGER(YRMIN_), INTEGER(YRMAX_), 
        ejXj, Hessian, n, p
    );

    calculate_test(REAL(ScoreR), Hessian, REAL(TestR), p);

    SET_STRING_ELT(list_names,  0,  mkChar("Score"));
    SET_STRING_ELT(list_names,  1,  mkChar("Hessian"));
    SET_STRING_ELT(list_names,  2,  mkChar("Test"));
    SET_VECTOR_ELT(list,        0,  ScoreR);
    SET_VECTOR_ELT(list,        1,  HessianR);
    SET_VECTOR_ELT(list,        2,  TestR);
    setAttrib(list, R_NamesSymbol,  list_names);

    free(eta);
    free(ejXj);
    free(Hessian);
    free(X);
    UNPROTECT(5);
    return list;
}

void calculate_D(double *exb, int *status, double *eta, int *yrmax, double **D, int n) {
    int i, j, k;

    for ( i = 0; i < n; i++) {
        D[i][i] = 0.0;
        for ( j = 0; j < n; j++ ) {
            D[i][j] = 0.0;
            D[j][i] = 0.0;
        }
    }
    
    double dtmp1 = 0.0;
    double dtmp2 = 0.0;
    double dtmp3;
    j = 0;
    for ( i = 0; i < n; i++ ) {
        while (j <= yrmax[i]) {
            if (status[i]) {
                // dtmp3 = status[j] / eta[j];
                dtmp3 = 1.0 / eta[j];
                dtmp1 += dtmp3;
                dtmp2 += dtmp3*dtmp3;
            }
            
            j++;
        }

        // update D
        dtmp3 = -exb[i]*dtmp2;
        D[i][i] = exb[i]*(dtmp1 + dtmp3);
        for ( k = i+1; k < n; k++ ) {
            D[i][k] = exb[k]*dtmp3;
            D[k][i] = D[i][k];
        }
    }
}

void calculate_test_from_start(
    double **X, double **D, double *v4score, double *score, double **Hessian, 
    double *test, int n, int p) {
    int i, j, i1, i2;

    for ( i = 0; i < p; i++) {
        score[i] = 0.0;
        Hessian[i][i] = 0.0;
        for ( j = i+1; j < p; j++) {
            Hessian[i][j] = 0.0;
            Hessian[j][i] = 0.0;
        }
    }

    for ( i = 0; i < n; i++) {
        for ( j = 0; j < p; j++) score[j] += v4score[i] * X[i][j];

        for ( i1 = 0; i1 < p; i1++ ) {
            Hessian[i1][i1] += D[i][i] * X[i][i1] * X[i][i1];
            for ( i2 = i1+1; i2 < p; i2++ ) {
                Hessian[i1][i2] += D[i][i] * X[i][i1] * X[i][i2];
            }
        }

        for ( j = i+1; j < n; j++) {
            for ( i1 = 0; i1 < p; i1++ ) {
                Hessian[i1][i1] += 2.0 * D[i][j] * X[i][i1] * X[j][i1];
                for ( i2 = i1+1; i2 < p; i2++ ) {
                    Hessian[i1][i2] += D[i][j] * (X[i][i1] * X[j][i2] + X[j][i1] * X[i][i2]);
                }
            }
        }
    }

    for ( i1 = 0; i1 < p; i1++ ) {
        for ( i2 = i1+1; i2 < p; i2++ ) {
            Hessian[i2][i1] = Hessian[i1][i2];
        }
    }
    MatrixInvSymmetric(Hessian, p);

    // calculate_test
    test[0] = 0.0;
    for ( i = 0; i < p; i++ ) {
        test[0] += Hessian[i][i]*score[i]*score[i];
        for ( j = i+1; j < p; j++ ) {
            test[0] += 2.0*Hessian[i][j]*score[i]*score[j];
        }
    }

}

SEXP _calculate_v4score_and_D2(SEXP X_, SEXP Exb_, SEXP Status_, SEXP YRMIN_, SEXP YRMAX_, SEXP N, SEXP P) {
    int n = INTEGER(N)[0];
    int p  = INTEGER(P)[0];

    double *eta;
    eta  = (double*)malloc(sizeof(double)*n);

    SEXP Vec, Dvec, list, list_names;
    SEXP ScoreR, HessianR, TestR;
    PROTECT(Vec = allocVector(REALSXP,  n));
    PROTECT(Dvec = allocVector(REALSXP,  n*n));
    PROTECT(ScoreR      = allocVector(REALSXP,  p));
    PROTECT(HessianR    = allocVector(REALSXP, p*p));
    PROTECT(TestR       = allocVector(REALSXP,  1));
    PROTECT(list_names  = allocVector(STRSXP,   5));
    PROTECT(list        = allocVector(VECSXP,   5));
    double **D = vec2matrix(REAL(Dvec), n, n);
    double **X = vec2matrix(REAL(X_), n, p);
    double **Hessian = vec2matrix(REAL(HessianR), p, p);
    
    calculate_v4score(
        REAL(Exb_), INTEGER(Status_), INTEGER(YRMIN_), INTEGER(YRMAX_), eta, REAL(Vec), n
    );
    calculate_D(REAL(Exb_), INTEGER(Status_), eta, INTEGER(YRMAX_), D, n);
    calculate_test_from_start(X, D, REAL(Vec), REAL(ScoreR), Hessian, REAL(TestR), n, p);

    SET_STRING_ELT(list_names,  0,  mkChar("v4score"));
    SET_STRING_ELT(list_names,  1,  mkChar("D"));
    SET_STRING_ELT(list_names,  2,  mkChar("ScoreR"));
    SET_STRING_ELT(list_names,  3,  mkChar("HessianinvR"));
    SET_STRING_ELT(list_names,  4,  mkChar("TestR"));
    SET_VECTOR_ELT(list,        0,  Vec);
    SET_VECTOR_ELT(list,        1,  Dvec);
    SET_VECTOR_ELT(list,        2,  ScoreR);
    SET_VECTOR_ELT(list,        3,  HessianR);
    SET_VECTOR_ELT(list,        4,  TestR);
    setAttrib(list, R_NamesSymbol,  list_names);

    free(D);
    free(eta);
    free(X);
    free(Hessian);
    UNPROTECT(7);
    return list;
}

SEXP _calculate_v4score_and_D(SEXP Exb_, SEXP Status_, SEXP YRMIN_, SEXP YRMAX_, SEXP N) {
    int n = INTEGER(N)[0];

    double *eta;
    eta  = (double*)malloc(sizeof(double)*n);

    SEXP Vec, Dvec, list, list_names;
    PROTECT(Vec = allocVector(REALSXP,  n));
    PROTECT(Dvec = allocVector(REALSXP,  n*n));
    PROTECT(list_names  = allocVector(STRSXP,   2));
    PROTECT(list        = allocVector(VECSXP,   2));
    double **D = vec2matrix(REAL(Dvec), n, n);
    
    calculate_v4score(
        REAL(Exb_), INTEGER(Status_), INTEGER(YRMIN_), INTEGER(YRMAX_), eta, REAL(Vec), n
    );
    calculate_D(REAL(Exb_), INTEGER(Status_), eta, INTEGER(YRMAX_), D, n);

    SET_STRING_ELT(list_names,  0,  mkChar("v4score"));
    SET_STRING_ELT(list_names,  1,  mkChar("D"));
    SET_VECTOR_ELT(list,        0,  Vec);
    SET_VECTOR_ELT(list,        1,  Dvec);
    setAttrib(list, R_NamesSymbol,  list_names);

    free(D);
    free(eta);
    UNPROTECT(4);
    return list;
}

void calculate_SUPtest(
    double **X, double **D, double *v4score, double *score, double **Hessian, 
    double *test, int n, int p, int px, double **Hessiancopy) {
    int i, j, i1, i2;

    for ( i = 0; i < p; i++) {
        score[i] = 0.0;
        Hessian[i][i] = 0.0;
        for ( j = i+1; j < p; j++) {
            Hessian[i][j] = 0.0;
            Hessian[j][i] = 0.0;
        }
    }

    for ( i = 0; i < n; i++) {
        for ( j = 0; j < p; j++) score[j] += v4score[i] * X[i][j];

        for ( i1 = 0; i1 < p; i1++ ) {
            Hessian[i1][i1] += D[i][i] * X[i][i1] * X[i][i1];
            for ( i2 = i1+1; i2 < p; i2++ ) {
                Hessian[i1][i2] += D[i][i] * X[i][i1] * X[i][i2];
            }
        }

        for ( j = i+1; j < n; j++) {
            for ( i1 = 0; i1 < p; i1++ ) {
                Hessian[i1][i1] += 2.0 * D[i][j] * X[i][i1] * X[j][i1];
                for ( i2 = i1+1; i2 < p; i2++ ) {
                    Hessian[i1][i2] += D[i][j] * (X[i][i1] * X[j][i2] + X[j][i1] * X[i][i2]);
                }
            }
        }
    }

    for ( i1 = 0; i1 < p; i1++ ) {
        for ( i2 = i1+1; i2 < p; i2++ ) {
            Hessian[i2][i1] = Hessian[i1][i2];
        }
    }
    for ( i1 = 0; i1 < px; i1++ ) {
        for ( i2 = 0; i2 < px; i2++ ) {
            Hessiancopy[i1][i2] = Hessian[i1][i2];
        }
    }
    MatrixInvSymmetric(Hessian, p);

    // calculate_test
    test[0] = 0.0;
    for ( i = 0; i < p; i++ ) {
        test[0] += Hessian[i][i]*score[i]*score[i];
        for ( j = i+1; j < p; j++ ) {
            test[0] += 2.0*Hessian[i][j]*score[i]*score[j];
        }
    }

}

void calculate_test_bootstrap(
    double **X, double **D, double *v4score, double *score, double **Hessian, 
    double *test, int n, int p, int **sample, int B, int px, double **Hessiancopy) {
    int i, j, i1, i2, b;
    double *xi, *xj;

    for ( b = 0; b < B; b++) {
        for ( i = px; i < p; i++) {
            score[i] = 0.0;
            Hessian[i][i] = 0.0;
            for ( j = i+1; j < p; j++) {
                Hessian[i][j] = 0.0;
                Hessian[j][i] = 0.0;
            }
            for ( j = 0; j < px; j++) {
                Hessian[i][j] = 0.0;
                Hessian[j][i] = 0.0;
            }
        }
    
        for ( i = 0; i < n; i++) {
            for ( j = px; j < p; j++) score[j] += v4score[i] * X[sample[i][b]][j];

            xi = X[sample[i][b]];
            for ( i1 = px; i1 < p; i1++ ) {
                Hessian[i1][i1] += D[i][i] * xi[i1] * xi[i1];
                for ( i2 = i1+1; i2 < p; i2++ ) {
                    Hessian[i1][i2] += D[i][i] * xi[i1] * xi[i2];
                }
                for ( i2 = 0; i2 < px; i2++ ) {
                    Hessian[i2][i1] += D[i][i] * xi[i1] * X[i][i2];
                }
            }

            for ( j = i+1; j < n; j++) {
                xj = X[sample[j][b]];
                for ( i1 = px; i1 < p; i1++ ) {
                    Hessian[i1][i1] += 2.0 * D[i][j] * xi[i1] * xj[i1];
                    for ( i2 = i1+1; i2 < p; i2++ ) {
                        Hessian[i1][i2] += D[i][j] * (xi[i1] * xj[i2] + xj[i1] * xi[i2]);
                    }
                    for ( i2 = 0; i2 < px; i2++ ) {
                        Hessian[i2][i1] += D[i][j] * (xi[i1] * X[j][i2] + xj[i1] * X[i][i2]);
                    }
                }
            }
        }
        
        for ( i1 = 0; i1 < px; i1++ ) {
            for ( i2 = i1; i2 < px; i2++ ) {
                Hessian[i1][i2] = Hessiancopy[i1][i2];
            }
        }
        for ( i1 = 0; i1 < p; i1++ ) {
            for ( i2 = i1+1; i2 < p; i2++ ) {
                Hessian[i2][i1] = Hessian[i1][i2];
            }
        }
        MatrixInvSymmetric(Hessian, p);
        
        // calculate_test
        test[b] = 0.0;
        for ( i = 0; i < p; i++ ) {
            test[b] += Hessian[i][i]*score[i]*score[i];
            for ( j = i+1; j < p; j++ ) {
                test[b] += 2.0*Hessian[i][j]*score[i]*score[j];
            }
        }
    }
}

SEXP _calculate_SUP3(
    SEXP X_, SEXP Exb_, SEXP Status_, SEXP YRMIN_, SEXP YRMAX_, SEXP N, SEXP P, 
    SEXP Sample, SEXP B, SEXP PX) {
    int n  = INTEGER(N)[0];
    int p  = INTEGER(P)[0];
    int px = INTEGER(PX)[0];
    int b  = INTEGER(B)[0];

    double *eta, *Hessiancopy;
    eta  = (double*)malloc(sizeof(double)*n);
    Hessiancopy  = (double*)malloc(sizeof(double)*p*p);

    SEXP Vec, Dvec, ScoreR, HessianR, TestR, TestB, list, list_names;
    PROTECT(Vec = allocVector(REALSXP,  n));
    PROTECT(Dvec = allocVector(REALSXP,  n*n));
    PROTECT(ScoreR      = allocVector(REALSXP,  p));
    PROTECT(HessianR    = allocVector(REALSXP, p*p));
    PROTECT(TestR       = allocVector(REALSXP,  1));
    PROTECT(TestB       = allocVector(REALSXP,  b));
    PROTECT(list_names  = allocVector(STRSXP,   6));
    PROTECT(list        = allocVector(VECSXP,   6));
    double **D = vec2matrix(REAL(Dvec), n, n);
    double **X = vec2matrix(REAL(X_), n, p);
    double **Hessian = vec2matrix(REAL(HessianR), p, p);
    int **sample = ivec2matrix(INTEGER(Sample), n, b);
    double **HessiancopyM = vec2matrix(Hessiancopy, p, p);
    
    calculate_v4score(
        REAL(Exb_), INTEGER(Status_), INTEGER(YRMIN_), INTEGER(YRMAX_), eta, REAL(Vec), n
    );
    calculate_D(REAL(Exb_), INTEGER(Status_), eta, INTEGER(YRMAX_), D, n);
    calculate_SUPtest(X, D, REAL(Vec), REAL(ScoreR), Hessian, REAL(TestR), n, p, px, HessiancopyM);
    calculate_test_bootstrap(X, D, REAL(Vec), REAL(ScoreR), Hessian, REAL(TestB), n, p, sample, b, px, HessiancopyM);

    SET_STRING_ELT(list_names,  0,  mkChar("v4score"));
    SET_STRING_ELT(list_names,  1,  mkChar("D"));
    SET_STRING_ELT(list_names,  2,  mkChar("ScoreR"));
    SET_STRING_ELT(list_names,  3,  mkChar("HessianinvR"));
    SET_STRING_ELT(list_names,  4,  mkChar("TestR"));
    SET_STRING_ELT(list_names,  5,  mkChar("TestB"));
    SET_VECTOR_ELT(list,        0,  Vec);
    SET_VECTOR_ELT(list,        1,  Dvec);
    SET_VECTOR_ELT(list,        2,  ScoreR);
    SET_VECTOR_ELT(list,        3,  HessianR);
    SET_VECTOR_ELT(list,        4,  TestR);
    SET_VECTOR_ELT(list,        5,  TestB);
    setAttrib(list, R_NamesSymbol,  list_names);

    free(D);
    free(eta);
    free(X);
    free(Hessian);
    free(Hessiancopy);
    free(HessiancopyM);
    free(sample);
    UNPROTECT(8);
    return list;
}

SEXP _calculate_SUP2(
    SEXP X_, SEXP Exb_, SEXP Status_, SEXP YRMIN_, SEXP YRMAX_, SEXP N, SEXP P, 
    SEXP Sample, SEXP B, SEXP PX) {
    int n  = INTEGER(N)[0];
    int p  = INTEGER(P)[0];
    int px = INTEGER(PX)[0];
    int b  = INTEGER(B)[0];

    double *eta, *Hessiancopy, *Vec, *Dvec, *ScoreR, *HessianR;
    eta  = (double*)malloc(sizeof(double)*n);
    Hessiancopy  = (double*)malloc(sizeof(double)*p*p);
    Vec  = (double*)malloc(sizeof(double)*n);
    Dvec  = (double*)malloc(sizeof(double)*n*n);
    ScoreR = (double*)malloc(sizeof(double)*p);
    HessianR = (double*)malloc(sizeof(double)*p*p);

    SEXP TestR, TestB, list, list_names;
    PROTECT(TestR       = allocVector(REALSXP,  1));
    PROTECT(TestB       = allocVector(REALSXP,  b));
    PROTECT(list_names  = allocVector(STRSXP,   4));
    PROTECT(list        = allocVector(VECSXP,   4));
    double **D = vec2matrix(Dvec, n, n);
    double **X = vec2matrix(REAL(X_), n, p);
    double **Hessian = vec2matrix(HessianR, p, p);
    int **sample = ivec2matrix(INTEGER(Sample), n, b);
    double **HessiancopyM = vec2matrix(Hessiancopy, p, p);
    
    calculate_v4score(
        REAL(Exb_), INTEGER(Status_), INTEGER(YRMIN_), INTEGER(YRMAX_), eta, Vec, n
    );
    calculate_D(REAL(Exb_), INTEGER(Status_), eta, INTEGER(YRMAX_), D, n);
    calculate_SUPtest(X, D, Vec, ScoreR, Hessian, REAL(TestR), n, p, px, HessiancopyM);
    calculate_test_bootstrap(X, D, Vec, ScoreR, Hessian, REAL(TestB), n, p, sample, b, px, HessiancopyM);


    SET_STRING_ELT(list_names,  0,  mkChar("TestR"));
    SET_STRING_ELT(list_names,  1,  mkChar("TestB"));
    SET_VECTOR_ELT(list,        0,  TestR);
    SET_VECTOR_ELT(list,        1,  TestB);
    setAttrib(list, R_NamesSymbol,  list_names);

    free(ScoreR);
    free(D);
    free(Dvec);
    free(Vec);
    free(eta);
    free(X);
    free(Hessian);
    free(Hessiancopy);
    free(HessiancopyM);
    free(sample);
    UNPROTECT(4);
    return list;
}

SEXP _calculate_SUP(
    SEXP X_, SEXP N, SEXP P, 
    SEXP Sample, SEXP B, SEXP PX, SEXP Vec, SEXP Dvec) {
    int n  = INTEGER(N)[0];
    int p  = INTEGER(P)[0];
    int px = INTEGER(PX)[0];
    int b  = INTEGER(B)[0];

    double *eta, *Hessiancopy, *ScoreR, *HessianR;
    eta  = (double*)malloc(sizeof(double)*n);
    Hessiancopy  = (double*)malloc(sizeof(double)*p*p);
    ScoreR = (double*)malloc(sizeof(double)*p);
    HessianR = (double*)malloc(sizeof(double)*p*p);

    SEXP TestR, TestB, list, list_names;
    PROTECT(TestR       = allocVector(REALSXP,  1));
    PROTECT(TestB       = allocVector(REALSXP,  b));
    PROTECT(list_names  = allocVector(STRSXP,   4));
    PROTECT(list        = allocVector(VECSXP,   4));
    double **D = vec2matrix(REAL(Dvec), n, n);
    double **X = vec2matrix(REAL(X_), n, p);
    double **Hessian = vec2matrix(HessianR, p, p);
    int **sample = ivec2matrix(INTEGER(Sample), n, b);
    double **HessiancopyM = vec2matrix(Hessiancopy, p, p);
    
    calculate_SUPtest(X, D, REAL(Vec), ScoreR, Hessian, REAL(TestR), n, p, px, HessiancopyM);
    calculate_test_bootstrap(X, D, REAL(Vec), ScoreR, Hessian, REAL(TestB), n, p, sample, b, px, HessiancopyM);

    SET_STRING_ELT(list_names,  0,  mkChar("TestR"));
    SET_STRING_ELT(list_names,  1,  mkChar("TestB"));
    SET_VECTOR_ELT(list,        0,  TestR);
    SET_VECTOR_ELT(list,        1,  TestB);
    setAttrib(list, R_NamesSymbol,  list_names);

    free(ScoreR);
    free(D);
    free(eta);
    free(X);
    free(Hessian);
    free(Hessiancopy);
    free(HessiancopyM);
    free(sample);
    UNPROTECT(4);
    return list;
}


// =============================================================================
// for EstCoxphNewton ----------------------------------------------------------
// =============================================================================
void calculate_exb(double *CoxCoef, double **x, int n, int px, double *exb) {
    int i, j;
    for ( i = 0; i < n; i++ ) {
        exb[i] = 0.0;
        for ( j = 0; j < px; j++ ) {
            exb[i] += x[i][j] * CoxCoef[j];
        }
        exb[i] = exp(exb[i]);
    }
}

void EstCoxphNewton(
    double *beta0, double **x, int *status, int *yrmin, int *yrmax, int n, int p, double tol
){
    int i, j;
    double *score, *exb, *ejXj, *eta, *Hessianvec;

    score = (double*)malloc(sizeof(double)*p);
    exb   = (double*)malloc(sizeof(double)*n);
    ejXj  = (double*)malloc(sizeof(double)*p);
    eta   = (double*)malloc(sizeof(double)*n);
    Hessianvec = (double*)malloc(sizeof(double)*p*p);
    double **Hessian = vec2matrix(Hessianvec, p, p);

    double tmp;
    double increment;

    for (j = 0; j < p; j++) beta0[j] = 0.0;

    increment = tol + 1.0;
    while (increment > tol) {
        calculate_exb(beta0, x, n, p, exb);
        calculate_score(x, exb, status, yrmin, yrmax, eta, score, n, p);
        calculate_Hessian_inv(x, exb, status, eta, yrmin, yrmax, ejXj, Hessian, n, p);
    
        increment = 0.0;
        for ( i = 0; i < p; i++ ) {
            tmp = 0.0;
            for ( j = 0; j < p; j++ ) tmp += Hessian[i][j] * score[j];
            if (fabs(tmp) > increment) increment = fabs(tmp);
            beta0[i] += tmp;
        }
    }

    free(Hessianvec);
    free(Hessian);
    free(eta);
    free(score);
    free(exb);
    free(ejXj);
}

SEXP _EstCoxphNewton(
    SEXP X_, SEXP Status_, SEXP YRMIN_, SEXP YRMAX_, SEXP N, SEXP P, SEXP TOL) {
    int n  = INTEGER(N)[0];
    int p  = INTEGER(P)[0];
    double tol   = REAL(TOL)[0];

    SEXP Beta, list, list_names;
    PROTECT(Beta        = allocVector(REALSXP,  p));
    PROTECT(list_names  = allocVector(STRSXP,   1));
    PROTECT(list        = allocVector(VECSXP,   1));
    double **X = vec2matrix(REAL(X_), n, p);

    EstCoxphNewton(
        REAL(Beta), X, INTEGER(Status_), INTEGER(YRMIN_), INTEGER(YRMAX_), n, p, tol
    );

    SET_STRING_ELT(list_names,  0,  mkChar("Beta"));
    SET_VECTOR_ELT(list,        0,  Beta);
    setAttrib(list, R_NamesSymbol,  list_names);

    free(X);
    UNPROTECT(3);
    return list;
}

// =============================================================================
// for SUP ---------------------------------------------------------------------
// =============================================================================
void calculate_xDx_xD_vzM_vX(
    double *Vec, double **D, double **x, double **z, double **xDx, double **xD, 
    double **vzM, double *vx, int n, int px, int pz
) {
    int i, j, k;

    // calculate xD
    for ( k = 0; k < n; k++ ) {
        for ( j = 0; j < px; j++ ) {
            xD[k][j] = 0.0;
            for ( i = 0; i < n; i++ ) {
                xD[k][j] += D[k][i] * x[i][j];
            }
        }
    }

    // calculate xDx
    for ( i = 0; i < px; i++ ) {
        for ( j = i; j < px; j++ ) {
            xDx[i][j] = 0.0;
            for ( k = 0; k < n; k++ ) {
                xDx[i][j] += x[k][i] * xD[k][j];
            }
        }
    }
    for ( i = 0; i < px; i++ ) {
        for ( j = i; j < px; j++ ) {
            xDx[j][i] = xDx[i][j];
        }
    }

    // calculate vx
    for ( j = 0; j < px; j++ ) {
        vx[j] = 0.0;
        for ( k = 0; k < n; k++ ) {
            vx[j] += Vec[k] * x[k][j];
        }
    }
    
    // calculate vzM
    for ( j = 0; j < pz; j++ ) {
        for ( k = 0; k < n; k++ ) {
            vzM[k][j] = Vec[k] * z[k][j];
        }
    }
}

void calculate_SUP_given_id(
    double **xDx, double **xD, double **D, double **z, double **vzM, 
    int *id, int n, int p, int px, int pz, double *Diz, 
    double *Score, double *Hessian, double **HessianM, double *test
) {
    int i, j, k, h;

    for ( i = px; i < p; i++ ) Score[i] = 0.0;
    for ( i = 0; i < p*p; i++ ) Hessian[i] = 0.0;

    for ( k = 0; k < n; k++ ) {
        if (id[k]) {
            for ( j = 0; j < pz; j++ ) Score[j+px] += vzM[k][j];
            for ( i = 0; i < pz; i++ ) {
                for ( j = 0; j < px; j++ ) {
                    HessianM[j][i+px] += xD[k][j] * z[k][i];
                }
            }

            for ( j = 0; j < pz; j++ ) Diz[j] = 0.0;
            for ( h = 0; h < n; h++ ) {
                if (id[h]) {
                    for ( j = 0; j < pz; j++ ) Diz[j] += D[k][h] * z[h][j];
                }
            }
            for ( i = 0; i < pz; i++ ) {
                for ( j = 0; j <= i; j++ ) {
                    HessianM[j+px][i+px] += Diz[j] * z[k][i];
                }
            }
        }
    }

    for ( i = 0; i < px; i++ ) {
        HessianM[i][i] = xDx[i][i];
        for ( j = i+1; j < px; j++ ) {
            HessianM[i][j] = xDx[i][j];
        }
    }

    for ( i = 0; i < p; i++ ) {
        for ( j = i+1; j < p; j++ ) {
            HessianM[j][i] = HessianM[i][j];
        }
    }

    // printArrayDouble(Score, p, p);
    // Rprintf("Hessian.\n");
    // printArrayDouble(Hessian, p*p, p);

    MatrixInvSymmetric(HessianM, p);

    // calculate_test
    test[0] = 0.0;
    for ( i = 0; i < p; i++ ) {
        test[0] += HessianM[i][i]*Score[i]*Score[i];
        for ( j = i+1; j < p; j++ ) {
            test[0] += 2.0*HessianM[i][j]*Score[i]*Score[j];
        }
    }
}

void calculate_SUP00(
    double **xDx, double **xD, double **D, double **z, double **vzM, double *vx, 
    int **idM, int **sample, int n, int p, int px, int pz, int K, int B, 
    double *test, double *testboot, double *pvalue
) {
    int i, j, bb;
    double testi;

    double *Score, *Hessian, *Diz;
    Score = (double*)malloc(sizeof(double)*p);
    Hessian = (double*)malloc(sizeof(double)*p*p);
    Diz = (double*)malloc(sizeof(double)*pz);
    double **HessianM = vec2matrix(Hessian, p, p);
    for ( i = 0; i < px; i++ ) Score[i] = vx[i];

    // SUP test staitsitic
    test[0] = 0.0;
    for ( i = 0; i < K; i++ ) {
        calculate_SUP_given_id(xDx, xD, D, z, vzM, idM[i], n, p, px, pz, Diz, Score, Hessian, HessianM, &testi);
        if (test[0] < testi) test[0] = testi;
    }

    // bootstrap
    int *idBoot;
    idBoot = (int*)malloc(sizeof(int)*n);
    pvalue[0] = 0.0;
    for ( bb = 0; bb < B; bb++ ) {
        testboot[bb] = 0.0;
        for ( i = 0; i < K; i++ ) {
            for ( j = 0; j < n; j++ ) idBoot[j] = idM[i][sample[j][bb]];
            calculate_SUP_given_id(xDx, xD, D, z, vzM, idBoot, n, p, px, pz, Diz, Score, Hessian, HessianM, &testi);
            if (testboot[bb] < testi) testboot[bb] = testi;
        }
        // pvalue
        if (testboot[bb] > test[0]) pvalue[0] += 1.0;
    }
    pvalue[0] /= B;
    
    free(idBoot);
    free(Score);
    free(Diz);
    free(Hessian);
    free(HessianM);
}

SEXP _calculate_SUP00(
    SEXP Status_, SEXP YRMIN_, SEXP YRMAX_, SEXP X_, SEXP Z_,  
    SEXP ID_, SEXP SAMPLE_, SEXP N, SEXP PZ, SEXP PX, SEXP B, SEXP K, SEXP TOL_) {
    int n  = INTEGER(N)[0];
    int pz = INTEGER(PZ)[0];
    int px = INTEGER(PX)[0];
    int b  = INTEGER(B)[0];
    int k  = INTEGER(K)[0];
    int p = pz + px;
    double tol = REAL(TOL_)[0];

    SEXP TestR, TestB, Pval, list, list_names;
    PROTECT(TestR       = allocVector(REALSXP,  1));
    PROTECT(TestB       = allocVector(REALSXP,  b));
    PROTECT(Pval        = allocVector(REALSXP,  1));
    PROTECT(list_names  = allocVector(STRSXP,   3));
    PROTECT(list        = allocVector(VECSXP,   3));

    double *eta, *Dvec, *Vec, *exb, *Beta;
    double **D;
    eta  = (double*)malloc(sizeof(double)*n);
    Vec = (double*)malloc(sizeof(double)*n);
    Dvec  = (double*)malloc(sizeof(double)*n*n);
    D = vec2matrix(Dvec, n, n);
    exb  = (double*)malloc(sizeof(double)*n);
    Beta  = (double*)malloc(sizeof(double)*px);
    double **x = vec2matrix(REAL(X_), n, px);
    EstCoxphNewton(Beta, x, INTEGER(Status_), INTEGER(YRMIN_), INTEGER(YRMAX_), n, px, tol);
    calculate_exb(Beta, x, n, px, exb);
    calculate_v4score(exb, INTEGER(Status_), INTEGER(YRMIN_), INTEGER(YRMAX_), eta, Vec, n);
    calculate_D(exb, INTEGER(Status_), eta, INTEGER(YRMAX_), D, n);
    free(eta);
    free(exb);

    double **z = vec2matrix(REAL(Z_), n, pz);
    double *xDxvec, *xDvec, *vzMvec, *vx;
    double **xDx, **xD, **vzM;
    vx      = (double*)malloc(sizeof(double)*px);
    xDxvec  = (double*)malloc(sizeof(double)*px*px);
    xDvec   = (double*)malloc(sizeof(double)*n*px);
    vzMvec  = (double*)malloc(sizeof(double)*n*pz);
    xDx     = vec2matrix(Dvec,   px, px);
    xD      = vec2matrix(xDvec,  n, px);
    vzM     = vec2matrix(vzMvec, n, pz);
    calculate_xDx_xD_vzM_vX(Vec, D, x, z, xDx, xD, vzM, vx, n, px, pz);
    free(Vec);
    free(x);

    int **sample = ivec2matrix(INTEGER(SAMPLE_), n, b);
    int **idM = ivec2matrix(INTEGER(ID_), k, n);
    calculate_SUP00(xDx, xD, D, z, vzM, vx, idM, sample, n, p, px, pz, k, b, REAL(TestR), REAL(TestB), REAL(Pval));

    SET_STRING_ELT(list_names,  0,  mkChar("TestR"));
    SET_STRING_ELT(list_names,  1,  mkChar("TestB"));
    SET_STRING_ELT(list_names,  2,  mkChar("Pvalue"));
    SET_VECTOR_ELT(list,        0,  TestR);
    SET_VECTOR_ELT(list,        1,  TestB);
    SET_VECTOR_ELT(list,        2,  Pval);
    setAttrib(list, R_NamesSymbol,  list_names);

    free(Beta);
    free(Dvec);
    free(D);
    free(z);
    free(xDxvec);
    free(xDx);
    free(xDvec);
    free(xD);
    free(vzMvec);
    free(vzM);
    free(vx);
    free(sample);
    free(idM);
    UNPROTECT(5);
    return list;
}

// =============================================================================
// for subgroup cox est --------------------------------------------------------
// =============================================================================
void calculate_exb_sub(double *CoxCoef, double **x, double **z, int *id, int n, int px, int pz, double *exb) {
    int i, j;
    for ( i = 0; i < n; i++ ) {
        exb[i] = 0.0;
        for ( j = 0; j < px; j++ ) {
            exb[i] += x[i][j] * CoxCoef[j];
        }
        if (id[i]) {
            for ( j = 0; j < pz; j++ ) {
                exb[i] += z[i][j] * CoxCoef[j+px];
            }
        }
        exb[i] = exp(exb[i]);
    }
}

void calculate_score_sub(
    double **X, double **Z, int *id, double *exb, int *status, int *yrmin, int *yrmax, double *eta, 
    double *score, int n, int px, int pz) {
    int i, j, k;
    double tmp;
    double vec;
    int p = px + pz;

    for ( i = 0; i < p; i++) score[i] = 0.0;

    tmp = 0.0;
    j = n-1;
    for ( i = n-1; i >= 0; i--) {
        while (j >= yrmin[i]) {
            tmp += exb[j];
            j--;
        }
        eta[i] = tmp;
    }

    tmp = 0.0;
    j = 0;
    for ( i = 0; i < n; i++) {
        while (j <= yrmax[i]) {
            tmp += status[j] / eta[j];
            j++;
        }
        vec = status[i] - tmp*exb[i];
        for ( k = 0; k < px; k++) score[k] += vec * X[i][k];
        if (id[i]) for ( k = 0; k < pz; k++) score[k+px] += vec * Z[i][k];
    }
}

void calculate_Hessian_inv_sub(
    double **x, double **z, int *id, double *exb, int *status, double *eta, int *yrmin, int *yrmax, 
    double *ejXj, double **Hessian, int n, int px, int pz
    ) {
    int i, j, i1, i2;
    double tmp;
    int p = px + pz;

    for ( i1 = 0; i1 < p; i1++ ) {
        Hessian[i1][i1] = 0.0;
        for ( i2 = i1+1; i2 < p; i2++ ) {
            Hessian[i1][i2] = 0.0;
        }
    }

    double dtmp1 = 0.0;
    j = 0;
    for ( i = 0; i < n; i++ ) {
        while (j <= yrmax[i]) {
            if (status[i]) {
                dtmp1 += 1.0 / eta[j];
            }
            j++;
        }

        tmp = exb[i] * dtmp1;
        for ( i1 = 0; i1 < px; i1++ ) {
            Hessian[i1][i1] += tmp * x[i][i1] * x[i][i1];
            for ( i2 = i1+1; i2 < px; i2++ ) {
                Hessian[i1][i2] += tmp * x[i][i1] * x[i][i2];
            }
        }
        if (id[i]) {
            for ( i1 = 0; i1 < pz; i1++ ) {
                Hessian[i1+px][i1+px] += tmp * z[i][i1] * z[i][i1];
                for ( i2 = i1+1; i2 < pz; i2++ ) {
                    Hessian[i1+px][i2+px] += tmp * z[i][i1] * z[i][i2];
                }
                for ( i2 = 0; i2 < px; i2++ ) {
                    Hessian[i2][i1+px] += tmp * x[i][i2] * z[i][i1];
                }
            }
        }
    }

    for ( i1 = 0; i1 < p; i1++ ) ejXj[i1] = 0.0;
    j = n-1;
    for ( i = n-1; i >= 0; i-- ) {
        while (j >= yrmin[i]) {
            for ( i1 = 0; i1 < px; i1++ ) {
                ejXj[i1] += x[j][i1] * exb[j];
            }
            if (id[j]) {
                for ( i1 = 0; i1 < pz; i1++ ) {
                    ejXj[i1+px] += z[j][i1] * exb[j];
                }
            }
            j--;
        }

        if (status[i] == 0) continue;

        tmp = 1.0 / ( eta[i] * eta[i] );
        for ( i1 = 0; i1 < p; i1++ ) {
            Hessian[i1][i1] -= tmp * ejXj[i1] * ejXj[i1];
            for ( i2 = i1+1; i2 < p; i2++ ) {
                Hessian[i1][i2] -= tmp * ejXj[i1] * ejXj[i2];
            }
        }
    }

    for ( i1 = 0; i1 < p; i1++ ) {
        for ( i2 = i1+1; i2 < p; i2++ ) {
            Hessian[i2][i1] = Hessian[i1][i2];
        }
    }
    MatrixInvSymmetric(Hessian, p);
}

void EstCoxphNewtonSubgroup(
    double *beta0, double **x, double **z, int *id, int *status, int *yrmin, 
    int *yrmax, int n, int px, int pz, double tol
){
    int i, j;
    double *score, *exb, *ejXj, *eta, *Hessianvec;
    int p = px + pz;

    score = (double*)malloc(sizeof(double)*p);
    exb   = (double*)malloc(sizeof(double)*n);
    ejXj  = (double*)malloc(sizeof(double)*p);
    eta   = (double*)malloc(sizeof(double)*n);
    Hessianvec = (double*)malloc(sizeof(double)*p*p);
    double **Hessian = vec2matrix(Hessianvec, p, p);

    double tmp;
    double increment;

    for (j = 0; j < p; j++) beta0[j] = 0.0;

    increment = tol + 1.0;
    while (increment > tol) {
        calculate_exb_sub(beta0, x, z, id, n, px, pz, exb);
        calculate_score_sub(x, z, id, exb, status, yrmin, yrmax, eta, score, n, px, pz);
        calculate_Hessian_inv_sub(x, z, id, exb, status, eta, yrmin, yrmax, ejXj, Hessian, n, px, pz);

        increment = 0.0;
        for ( i = 0; i < p; i++ ) {
            tmp = 0.0;
            for ( j = 0; j < p; j++ ) tmp += Hessian[i][j] * score[j];
            if (fabs(tmp) > increment) increment = fabs(tmp);
            beta0[i] += tmp;
        }
    }

    free(Hessianvec);
    free(Hessian);
    free(eta);
    free(score);
    free(exb);
    free(ejXj);
}

SEXP _EstCoxphNewtonSubgroup(
    SEXP X_, SEXP Z_, SEXP ID_, SEXP Status_, SEXP YRMIN_, SEXP YRMAX_, SEXP N, SEXP PX, SEXP PZ, SEXP TOL) {
    int n   = INTEGER(N)[0];
    int px  = INTEGER(PX)[0];
    int pz  = INTEGER(PZ)[0];
    double tol   = REAL(TOL)[0];
    int p = px + pz;

    SEXP Beta, list, list_names;
    PROTECT(Beta        = allocVector(REALSXP,  p));
    PROTECT(list_names  = allocVector(STRSXP,   1));
    PROTECT(list        = allocVector(VECSXP,   1));
    double **X = vec2matrix(REAL(X_), n, px);
    double **Z = vec2matrix(REAL(Z_), n, pz);

    EstCoxphNewtonSubgroup(
        REAL(Beta), X, Z, INTEGER(ID_), INTEGER(Status_), INTEGER(YRMIN_), INTEGER(YRMAX_), n, px, pz, tol
    );

    SET_STRING_ELT(list_names,  0,  mkChar("Beta"));
    SET_VECTOR_ELT(list,        0,  Beta);
    setAttrib(list, R_NamesSymbol,  list_names);

    free(X);
    UNPROTECT(3);
    return list;
}
 
// =============================================================================
// for LRT ---------------------------------------------------------------------
// =============================================================================
void calculate_LRT_given_id(
    double *beta0, double *score, double *exb, double *ejXj, double *eta, 
    double **x, double **z, int *id, int *status, int *yrmin, 
    int *yrmax, int n, int px, int pz, int b, double tol, double **HessianH0, double **HessianH1minusH0, 
    double **rbyDsqrt, double **rbyDsqrtX, double **vzM,  
    double *ScoreH0, double *test, double *testboot
) {
    // EstCoxphNewtonSubgroup
    int i, j, bb;
    int p = px + pz;

    double tmp;
    double increment;

    for (j = 0; j < p; j++) beta0[j] = 0.0;

    increment = tol + 1.0;
    while (increment > tol) {
        calculate_exb_sub(beta0, x, z, id, n, px, pz, exb);
        calculate_score_sub(x, z, id, exb, status, yrmin, yrmax, eta, score, n, px, pz);
        calculate_Hessian_inv_sub(x, z, id, exb, status, eta, yrmin, yrmax, ejXj, HessianH1minusH0, n, px, pz);

        increment = 0.0;
        for ( i = 0; i < p; i++ ) {
            tmp = 0.0;
            for ( j = 0; j < p; j++ ) tmp += HessianH1minusH0[i][j] * score[j];
            if (fabs(tmp) > increment) increment = fabs(tmp);
            beta0[i] += tmp;
        }
    }
    calculate_exb_sub(beta0, x, z, id, n, px, pz, exb);
    calculate_score_sub(x, z, id, exb, status, yrmin, yrmax, eta, score, n, px, pz);
    calculate_Hessian_inv_sub(x, z, id, exb, status, eta, yrmin, yrmax, ejXj, HessianH1minusH0, n, px, pz);
    // HessianH1minusH0
    for ( i = 0; i < px; i++ ) {
        HessianH1minusH0[i][i] -= HessianH0[i][i];
        for (j = i+1; j < px; j++) {
            HessianH1minusH0[i][j] -= HessianH0[i][j];
            HessianH1minusH0[j][i] = HessianH1minusH0[i][j];
        }
    }

    // calculate ScoreH0
    for ( j = px; j < p; j++ ) ScoreH0[j] = 0.0;
    for ( i = 0; i < n; i++ ) {
        if (id[i]) {
            for ( j = 0; j < pz; j++ ) ScoreH0[j+px] += vzM[i][j];
        }
    }

    // calculate_test
    test[0] = 0.0;
    for ( i = 0; i < p; i++ ) {
        test[0] += HessianH1minusH0[i][i]*ScoreH0[i]*ScoreH0[i];
        for ( j = i+1; j < p; j++ ) {
            test[0] += 2.0*HessianH1minusH0[i][j]*ScoreH0[i]*ScoreH0[j];
        }
    }

    // bootstrap
    for ( bb = 0; bb < b; bb++) {
        for ( i = 0; i < px; i++ ) score[i] = rbyDsqrtX[bb][i];
        for ( i = 0; i < pz; i++ ) {
            score[i+px] = 0.0;
            for ( j = 0; j < n; j++ ) {
                if (id[j]) {
                    score[i+px] += rbyDsqrt[bb][j] * z[j][i];
                }
            }
        }

        // calculate_test
        testboot[bb] = 0.0;
        for ( i = 0; i < p; i++ ) {
            testboot[bb] += HessianH1minusH0[i][i]*score[i]*score[i];
            for ( j = i+1; j < p; j++ ) {
                testboot[bb] += 2.0*HessianH1minusH0[i][j]*score[i]*score[j];
            }
        }
    }

    // printArrayDouble(testboot, b, b);
    // printArrayDouble(test, 1, 1);
    // printArrayDouble2(HessianH1minusH0, p, p);
}

void calculate_LRT00(
    double **x, double **z, int **idM, int *status, int *yrmin, int *yrmax, int n, 
    int px, int pz, int b, double tol, double **HessianH0, int K, double *vx, 
    double **rbyDsqrt, double **rbyDsqrtX, double **vzM, double *test, double *testboot, double *pval
) {
    int i, j;
    int p = px + pz;

    double *beta0, *score, *exb, *ejXj, *eta, *Hessianvec;
    beta0 = (double*)malloc(sizeof(double)*p);
    score = (double*)malloc(sizeof(double)*p);
    exb   = (double*)malloc(sizeof(double)*n);
    ejXj  = (double*)malloc(sizeof(double)*p);
    eta   = (double*)malloc(sizeof(double)*n);
    Hessianvec = (double*)malloc(sizeof(double)*p*p);
    double **HessianH1minusH0 = vec2matrix(Hessianvec, p, p);

    double *ScoreH0;
    ScoreH0 = (double*)malloc(sizeof(double)*p);
    for ( i = 0; i < px; i++ ) ScoreH0[i] = vx[i];

    double testi;
    double *testbooti;
    testbooti = (double*)malloc(sizeof(double)*b);
    
    test[0] = 0.0;
    for ( j = 0; j < b; j++ ) testboot[j] = 0.0;
    
    for ( i = 0; i < K; i++ ) {
        calculate_LRT_given_id(
            beta0, score, exb, ejXj, eta, x, z, idM[i], status, yrmin, yrmax, 
            n, px, pz, b, tol, HessianH0, HessianH1minusH0, rbyDsqrt, rbyDsqrtX, 
            vzM, ScoreH0, &testi, testbooti
        );
        if (test[0] < testi) test[0] = testi;
        for ( j = 0; j < b; j++ ) {
            if (testboot[j] < testbooti[j]) testboot[j] = testbooti[j];
        }
    }

    // pval
    pval[0] = 0.0;
    for ( j = 0; j < b; j++ ) {
        if (testboot[j] > test[0]) pval[0] += 1.0;
    }
    pval[0] /= b;
    
    free(Hessianvec);
    free(HessianH1minusH0);
    free(beta0);
    free(eta);
    free(score);
    free(exb);
    free(ejXj);
}

void calculate_LRT_Related(
    double **x, double **z, double *exb, int *status, double *eta, int *yrmin, int *yrmax, 
    int n, int px, int pz, int b, double *vec, double **rvmatrix,
    double **HessianH0, double *vx, double **rbyDsqrt, double **rbyDsqrtX, double **vzM
    ) {
    int i, j, k;

    // Hessian_inv
    double *ejXj;
    ejXj  = (double*)malloc(sizeof(double)*px);
    calculate_Hessian_inv(x, exb, status, eta, yrmin, yrmax, ejXj, HessianH0, n, px);
    free(ejXj);

    // calculate vx
    for ( j = 0; j < px; j++ ) {
        vx[j] = 0.0;
        for ( k = 0; k < n; k++ ) {
            vx[j] += vec[k] * x[k][j];
        }
    }

    // calculate vzM
    for ( j = 0; j < pz; j++ ) {
        for ( k = 0; k < n; k++ ) {
            vzM[k][j] = vec[k] * z[k][j];
        }
    }

    double *Dvec;
    Dvec = (double*)malloc(sizeof(double)*n*n);
    double **Dsqrt = vec2matrix(Dvec, n, n);
    for ( i = 0; i < n*n; i++) Dvec[i] = 0.0;
    
    // calculate_Dsqrt
    double dtmp1 = 0.0;
    double dtmp2 = 0.0;
    double dtmp3;
    j = 0;
    for ( i = 0; i < n; i++ ) {
        while (j <= yrmax[i]) {
            if (status[i]) {
                // dtmp3 = status[j] / eta[j];
                dtmp3 = 1.0 / eta[j];
                dtmp1 += dtmp3;
                dtmp2 += dtmp3*dtmp3;
            }
            
            j++;
        }

        // update D
        dtmp3 = -exb[i]*dtmp2;
        Dsqrt[i][i] = exb[i]*(dtmp1 + dtmp3);
        for ( k = i+1; k < n; k++ ) {
            Dsqrt[i][k] = exb[k]*dtmp3;
        }
    }
    cholesky_decomposition_(Dsqrt, n);

    // rbyDsqrt = rvmatrix %*% Dsqrt, B by n matrix
    for ( i = 0; i < b; i++) {
        for ( j = 0; j < n; j++ ) {
            rbyDsqrt[i][j] = 0.0;
            for ( k = 0; k < n; k++ ) {
                rbyDsqrt[i][j] += rvmatrix[i][k] * Dsqrt[k][j];
            }
        }
    }

    // rbyDsqrtX = rbyDsqrt %*% x, B by px matrix
    for ( i = 0; i < b; i++) {
        for ( j = 0; j < px; j++ ) {
            rbyDsqrtX[i][j] = 0.0;
            for ( k = 0; k < n; k++ ) {
                rbyDsqrtX[i][j] += rbyDsqrt[i][k] * x[k][j];
            }
        }
    }
    
    free(Dvec);
    free(Dsqrt);
}

SEXP _calculate_LRT00(
    SEXP Status_, SEXP YRMIN_, SEXP YRMAX_, SEXP X_, SEXP Z_,  
    SEXP ID_, SEXP Rvmatrix_, SEXP N, SEXP PZ, SEXP PX, SEXP B, SEXP K, SEXP TOL
    ) {
    int n  = INTEGER(N)[0];
    int pz = INTEGER(PZ)[0];
    int px = INTEGER(PX)[0];
    int b  = INTEGER(B)[0];
    int k  = INTEGER(K)[0];
    // int p = pz + px;
    double tol = REAL(TOL)[0];
    double **x = vec2matrix(REAL(X_), n, px);
    double **z = vec2matrix(REAL(Z_), n, pz);
    int *yrmin = INTEGER(YRMIN_);
    int *yrmax = INTEGER(YRMAX_);

    SEXP TestR, TestB, Pval, list, list_names;
    PROTECT(TestR       = allocVector(REALSXP,  1));
    PROTECT(TestB       = allocVector(REALSXP,  b));
    PROTECT(Pval        = allocVector(REALSXP,  1));
    PROTECT(list_names  = allocVector(STRSXP,   3));
    PROTECT(list        = allocVector(VECSXP,   3));

    // calculate H0
    double *BetaH0, *exb, *eta, *Vec;
    BetaH0 = (double*)malloc(sizeof(double)*px);
    exb    = (double*)malloc(sizeof(double)*n);
    eta    = (double*)malloc(sizeof(double)*n);
    Vec    = (double*)malloc(sizeof(double)*n);
    EstCoxphNewton(BetaH0, x, INTEGER(Status_), yrmin, yrmax, n, px, tol);
    calculate_exb(BetaH0, x, n, px, exb);
    free(BetaH0);
    calculate_v4score(exb, INTEGER(Status_), yrmin, yrmax, eta, Vec, n);
    double *HessianH0Vec, *vx, *rbyDsqrtVec, *rbyDsqrtXVec, *vzMvec;
    HessianH0Vec = (double*)malloc(sizeof(double)*px*px);
    vx           = (double*)malloc(sizeof(double)*px);
    rbyDsqrtVec  = (double*)malloc(sizeof(double)*b*n);
    rbyDsqrtXVec = (double*)malloc(sizeof(double)*b*px);
    vzMvec       = (double*)malloc(sizeof(double)*n*pz);
    double **HessianH0 = vec2matrix(HessianH0Vec, px, px);
    double **rbyDsqrt  = vec2matrix(rbyDsqrtVec, b, n);
    double **rbyDsqrtX = vec2matrix(rbyDsqrtXVec, b, px);
    double **vzM       = vec2matrix(vzMvec, n, pz);

    double **rvmatrix = vec2matrix(REAL(Rvmatrix_), b, n);
    calculate_LRT_Related(
        x, z, exb, INTEGER(Status_), eta, yrmin, yrmax, n, px, pz, b, Vec, 
        rvmatrix, HessianH0, vx, rbyDsqrt, rbyDsqrtX, vzM
    );
    free(exb);
    free(eta);
    free(Vec);

    int **idM = ivec2matrix(INTEGER(ID_), k, n);
    calculate_LRT00(
        x, z, idM, INTEGER(Status_), yrmin, yrmax, n, px, pz, b, tol, HessianH0, 
        k, vx, rbyDsqrt, rbyDsqrtX, vzM, REAL(TestR), REAL(TestB), REAL(Pval)
    );

    SET_STRING_ELT(list_names,  0,  mkChar("TestR"));
    SET_STRING_ELT(list_names,  1,  mkChar("TestB"));
    SET_STRING_ELT(list_names,  2,  mkChar("Pvalue"));
    SET_VECTOR_ELT(list,        0,  TestR);
    SET_VECTOR_ELT(list,        1,  TestB);
    SET_VECTOR_ELT(list,        2,  Pval);
    setAttrib(list, R_NamesSymbol,  list_names);

    free(x);
    free(z);
    free(HessianH0Vec);
    free(HessianH0);
    free(vx);
    free(rbyDsqrtVec);
    free(rbyDsqrt);
    free(rbyDsqrtXVec);
    free(rbyDsqrtX);
    free(vzMvec);
    free(vzM);
    free(rvmatrix);
    free(idM);
    UNPROTECT(5);
    return list;
}

// =============================================================================
// for LRT ---------------------------------------------------------------------
// =============================================================================
void calculate_LRT_given_id_permutation(
    double *beta0, double *score, double *exb, double *ejXj, double *eta, 
    double **x, double **z, int *id, int *status, int *yrmin, 
    int *yrmax, int n, int px, int pz, int b, double tol, double **HessianH0, double **HessianH1minusH0, 
    int **sample, double **vzM,  
    double *ScoreH0, double *test, double *testboot
) {
    // EstCoxphNewtonSubgroup
    int i, j, bb;
    int p = px + pz;

    double tmp;
    double increment;

    for (j = 0; j < p; j++) beta0[j] = 0.0;

    increment = tol + 1.0;
    while (increment > tol) {
        calculate_exb_sub(beta0, x, z, id, n, px, pz, exb);
        calculate_score_sub(x, z, id, exb, status, yrmin, yrmax, eta, score, n, px, pz);
        calculate_Hessian_inv_sub(x, z, id, exb, status, eta, yrmin, yrmax, ejXj, HessianH1minusH0, n, px, pz);

        increment = 0.0;
        for ( i = 0; i < p; i++ ) {
            tmp = 0.0;
            for ( j = 0; j < p; j++ ) tmp += HessianH1minusH0[i][j] * score[j];
            if (fabs(tmp) > increment) increment = fabs(tmp);
            beta0[i] += tmp;
        }
    }
    calculate_exb_sub(beta0, x, z, id, n, px, pz, exb);
    calculate_score_sub(x, z, id, exb, status, yrmin, yrmax, eta, score, n, px, pz);
    calculate_Hessian_inv_sub(x, z, id, exb, status, eta, yrmin, yrmax, ejXj, HessianH1minusH0, n, px, pz);
    // HessianH1minusH0
    for ( i = 0; i < px; i++ ) {
        HessianH1minusH0[i][i] -= HessianH0[i][i];
        for (j = i+1; j < px; j++) {
            HessianH1minusH0[i][j] -= HessianH0[i][j];
            HessianH1minusH0[j][i] = HessianH1minusH0[i][j];
        }
    }

    // calculate ScoreH0
    for ( j = px; j < p; j++ ) ScoreH0[j] = 0.0;
    for ( i = 0; i < n; i++ ) {
        if (id[i]) {
            for ( j = 0; j < pz; j++ ) ScoreH0[j+px] += vzM[i][j];
        }
    }

    // calculate_test
    test[0] = 0.0;
    for ( i = 0; i < p; i++ ) {
        test[0] += HessianH1minusH0[i][i]*ScoreH0[i]*ScoreH0[i];
        for ( j = i+1; j < p; j++ ) {
            test[0] += 2.0*HessianH1minusH0[i][j]*ScoreH0[i]*ScoreH0[j];
        }
    }

    // bootstrap
    for ( i = 0; i < px; i++ ) score[i] = ScoreH0[i];
    for ( bb = 0; bb < b; bb++) {
        for ( i = px; i < p; i++ ) score[i] = 0.0;
        for ( i = 0; i < n; i++ ) {
            if (id[sample[i][bb]]) {
                for ( j = 0; j < pz; j++ ) score[j+px] += vzM[i][j];
            }
        }

        // calculate_test
        testboot[bb] = 0.0;
        for ( i = 0; i < p; i++ ) {
            testboot[bb] += HessianH1minusH0[i][i]*score[i]*score[i];
            for ( j = i+1; j < p; j++ ) {
                testboot[bb] += 2.0*HessianH1minusH0[i][j]*score[i]*score[j];
            }
        }
    }

    // printArrayDouble(testboot, b, b);
    // printArrayDouble(test, 1, 1);
    // printArrayDouble2(HessianH1minusH0, p, p);
}

void calculate_LRT00_permutation(
    double **x, double **z, int **idM, int *status, int *yrmin, int *yrmax, int n, 
    int px, int pz, int b, double tol, double **HessianH0, int K, double *vx, 
    int **sample, double **vzM, double *test, double *testboot, double *pval
) {
    int i, j;
    int p = px + pz;

    double *beta0, *score, *exb, *ejXj, *eta, *Hessianvec;
    beta0 = (double*)malloc(sizeof(double)*p);
    score = (double*)malloc(sizeof(double)*p);
    exb   = (double*)malloc(sizeof(double)*n);
    ejXj  = (double*)malloc(sizeof(double)*p);
    eta   = (double*)malloc(sizeof(double)*n);
    Hessianvec = (double*)malloc(sizeof(double)*p*p);
    double **HessianH1minusH0 = vec2matrix(Hessianvec, p, p);

    double *ScoreH0;
    ScoreH0 = (double*)malloc(sizeof(double)*p);
    for ( i = 0; i < px; i++ ) ScoreH0[i] = vx[i];

    double testi;
    double *testbooti;
    testbooti = (double*)malloc(sizeof(double)*b);
    
    test[0] = 0.0;
    for ( j = 0; j < b; j++ ) testboot[j] = 0.0;
    
    for ( i = 0; i < K; i++ ) {
        calculate_LRT_given_id_permutation(
            beta0, score, exb, ejXj, eta, x, z, idM[i], status, yrmin, yrmax, 
            n, px, pz, b, tol, HessianH0, HessianH1minusH0, sample, 
            vzM, ScoreH0, &testi, testbooti
        );
        if (test[0] < testi) test[0] = testi;
        for ( j = 0; j < b; j++ ) {
            if (testboot[j] < testbooti[j]) testboot[j] = testbooti[j];
        }
    }

    // pval
    pval[0] = 0.0;
    for ( j = 0; j < b; j++ ) {
        if (testboot[j] > test[0]) pval[0] += 1.0;
    }
    pval[0] /= b;
    
    free(Hessianvec);
    free(HessianH1minusH0);
    free(beta0);
    free(eta);
    free(score);
    free(exb);
    free(ejXj);
}

void calculate_LRT_Related_permutation(
    double **x, double **z, double *exb, int *status, double *eta, int *yrmin, int *yrmax, 
    int n, int px, int pz, int b, double *vec,
    double **HessianH0, double *vx, double **vzM
    ) {
    int j, k;

    // Hessian_inv
    double *ejXj;
    ejXj  = (double*)malloc(sizeof(double)*px);
    calculate_Hessian_inv(x, exb, status, eta, yrmin, yrmax, ejXj, HessianH0, n, px);
    free(ejXj);

    // calculate vx
    for ( j = 0; j < px; j++ ) {
        vx[j] = 0.0;
        for ( k = 0; k < n; k++ ) {
            vx[j] += vec[k] * x[k][j];
        }
    }

    // calculate vzM
    for ( j = 0; j < pz; j++ ) {
        for ( k = 0; k < n; k++ ) {
            vzM[k][j] = vec[k] * z[k][j];
        }
    }
}

SEXP _calculate_LRT00_permutation(
    SEXP Status_, SEXP YRMIN_, SEXP YRMAX_, SEXP X_, SEXP Z_,  
    SEXP ID_, SEXP SAMPLE_, SEXP N, SEXP PZ, SEXP PX, SEXP B, SEXP K, SEXP TOL
    ) {
    int n  = INTEGER(N)[0];
    int pz = INTEGER(PZ)[0];
    int px = INTEGER(PX)[0];
    int b  = INTEGER(B)[0];
    int k  = INTEGER(K)[0];
    // int p = pz + px;
    double tol = REAL(TOL)[0];
    double **x = vec2matrix(REAL(X_), n, px);
    double **z = vec2matrix(REAL(Z_), n, pz);
    int *yrmin = INTEGER(YRMIN_);
    int *yrmax = INTEGER(YRMAX_);

    SEXP TestR, TestB, Pval, list, list_names;
    PROTECT(TestR       = allocVector(REALSXP,  1));
    PROTECT(TestB       = allocVector(REALSXP,  b));
    PROTECT(Pval        = allocVector(REALSXP,  1));
    PROTECT(list_names  = allocVector(STRSXP,   3));
    PROTECT(list        = allocVector(VECSXP,   3));

    // calculate H0
    double *BetaH0, *exb, *eta, *Vec;
    BetaH0 = (double*)malloc(sizeof(double)*px);
    exb    = (double*)malloc(sizeof(double)*n);
    eta    = (double*)malloc(sizeof(double)*n);
    Vec    = (double*)malloc(sizeof(double)*n);
    EstCoxphNewton(BetaH0, x, INTEGER(Status_), yrmin, yrmax, n, px, tol);
    calculate_exb(BetaH0, x, n, px, exb);
    free(BetaH0);
    calculate_v4score(exb, INTEGER(Status_), yrmin, yrmax, eta, Vec, n);
    double *HessianH0Vec, *vx, *vzMvec;
    HessianH0Vec = (double*)malloc(sizeof(double)*px*px);
    vx           = (double*)malloc(sizeof(double)*px);
    vzMvec       = (double*)malloc(sizeof(double)*n*pz);
    double **HessianH0 = vec2matrix(HessianH0Vec, px, px);
    double **vzM       = vec2matrix(vzMvec, n, pz);

    calculate_LRT_Related_permutation(
        x, z, exb, INTEGER(Status_), eta, yrmin, yrmax, n, px, pz, b, Vec, 
        HessianH0, vx, vzM
    );
    free(exb);
    free(eta);
    free(Vec);

    int **idM = ivec2matrix(INTEGER(ID_), k, n);
    int **sample = ivec2matrix(INTEGER(SAMPLE_), n, b);
    calculate_LRT00_permutation(
        x, z, idM, INTEGER(Status_), yrmin, yrmax, n, px, pz, b, tol, HessianH0, 
        k, vx, sample, vzM, REAL(TestR), REAL(TestB), REAL(Pval)
    );

    SET_STRING_ELT(list_names,  0,  mkChar("TestR"));
    SET_STRING_ELT(list_names,  1,  mkChar("TestB"));
    SET_STRING_ELT(list_names,  2,  mkChar("Pvalue"));
    SET_VECTOR_ELT(list,        0,  TestR);
    SET_VECTOR_ELT(list,        1,  TestB);
    SET_VECTOR_ELT(list,        2,  Pval);
    setAttrib(list, R_NamesSymbol,  list_names);

    free(x);
    free(z);
    free(HessianH0Vec);
    free(HessianH0);
    free(vx);
    free(vzMvec);
    free(vzM);
    free(sample);
    free(idM);
    UNPROTECT(5);
    return list;
}


// =============================================================================
// for ST ---------------------------------------------------------------------
// =============================================================================
void calculate_ST_given_id(
    double **z, double **vzM, double **rvmatrix, 
    int *id, int n, int pz, int b, 
    double *Score, double *Hessian, double **HessianM, double *test, double *testboot
) {
    int i, j, k, bb;

    for ( i = 0; i < pz; i++ ) Score[i] = 0.0;
    for ( i = 0; i < pz*pz; i++ ) Hessian[i] = 0.0;

    for ( k = 0; k < n; k++ ) {
        if (id[k]) {
            for ( j = 0; j < pz; j++ ) Score[j] += vzM[k][j];
            for ( i = 0; i < pz; i++ ) {
                for ( j = i; j < pz; j++ ) {
                    HessianM[i][j] += vzM[k][i] * vzM[k][j];
                }
            }
        }
    }

    for ( i = 0; i < pz; i++ ) {
        for ( j = i+1; j < pz; j++ ) {
            HessianM[j][i] = HessianM[i][j];
        }
    }

    MatrixInvSymmetric(HessianM, pz);

    // calculate_test
    test[0] = 0.0;
    for ( i = 0; i < pz; i++ ) {
        test[0] += HessianM[i][i]*Score[i]*Score[i];
        for ( j = i+1; j < pz; j++ ) {
            test[0] += 2.0*HessianM[i][j]*Score[i]*Score[j];
        }
    }

    // bootstrap
    for ( bb = 0; bb < b; bb++ ) {
        for ( i = 0; i < pz; i++ ) Score[i] = 0.0;
        for ( i = 0; i < n; i++ ) {
            if (id[i]) {
                for ( j = 0; j < pz; j++ ) Score[j] += rvmatrix[bb][i] * vzM[i][j];
            }
        }
        
        // calculate_test
        testboot[bb] = 0.0;
        for ( i = 0; i < pz; i++ ) {
            testboot[bb] += HessianM[i][i]*Score[i]*Score[i];
            for ( j = i+1; j < pz; j++ ) {
                testboot[bb] += 2.0*HessianM[i][j]*Score[i]*Score[j];
            }
        }
    }
}

void calculate_ST00(
    double **z, int **idM, int n, int pz, int b, int K, double **rvmatrix, 
    double **vzM, double *test, double *testboot, double *pval
) {
    int i, j;

    double *score, *Hessianvec;
    score = (double*)malloc(sizeof(double)*pz);
    Hessianvec = (double*)malloc(sizeof(double)*pz*pz);
    double **Hessian = vec2matrix(Hessianvec, pz, pz);

    double testi;
    double *testbooti;
    testbooti = (double*)malloc(sizeof(double)*b);
    
    test[0] = 0.0;
    for ( j = 0; j < b; j++ ) testboot[j] = 0.0;
    
    for ( i = 0; i < K; i++ ) {
        calculate_ST_given_id(
            z, vzM, rvmatrix, idM[i], n, pz, b, score, Hessianvec, Hessian, &testi, testbooti
        );
        if (test[0] < testi) test[0] = testi;
        for ( j = 0; j < b; j++ ) {
            if (testboot[j] < testbooti[j]) testboot[j] = testbooti[j];
        }
    }

    // pval
    pval[0] = 0.0;
    for ( j = 0; j < b; j++ ) {
        if (testboot[j] > test[0]) pval[0] += 1.0;
    }
    pval[0] /= b;
    
    free(Hessianvec);
    free(Hessian);
    free(score);
}

SEXP _calculate_ST00(
    SEXP Status_, SEXP YRMIN_, SEXP YRMAX_, SEXP X_, SEXP Z_,  
    SEXP ID_, SEXP Rvmatrix_, SEXP N, SEXP PZ, SEXP PX, SEXP B, SEXP K, SEXP TOL
    ) {
    int n  = INTEGER(N)[0];
    int pz = INTEGER(PZ)[0];
    int px = INTEGER(PX)[0];
    int b  = INTEGER(B)[0];
    int k  = INTEGER(K)[0];
    // int p = pz + px;
    double tol = REAL(TOL)[0];
    double **x = vec2matrix(REAL(X_), n, px);
    double **z = vec2matrix(REAL(Z_), n, pz);
    int *yrmin = INTEGER(YRMIN_);
    int *yrmax = INTEGER(YRMAX_);

    SEXP TestR, TestB, Pval, list, list_names;
    PROTECT(TestR       = allocVector(REALSXP,  1));
    PROTECT(TestB       = allocVector(REALSXP,  b));
    PROTECT(Pval        = allocVector(REALSXP,  1));
    PROTECT(list_names  = allocVector(STRSXP,   3));
    PROTECT(list        = allocVector(VECSXP,   3));

    // calculate H0
    double *BetaH0, *exb, *eta, *Vec, *vzMvec;
    BetaH0       = (double*)malloc(sizeof(double)*px);
    exb          = (double*)malloc(sizeof(double)*n);
    eta          = (double*)malloc(sizeof(double)*n);
    Vec          = (double*)malloc(sizeof(double)*n);
    vzMvec       = (double*)malloc(sizeof(double)*n*pz);
    double **vzM = vec2matrix(vzMvec, n, pz);
    EstCoxphNewton(BetaH0, x, INTEGER(Status_), yrmin, yrmax, n, px, tol);
    calculate_exb(BetaH0, x, n, px, exb);
    free(BetaH0);
    calculate_v4score(exb, INTEGER(Status_), yrmin, yrmax, eta, Vec, n);
    // calculate_ST_Related
    for ( int j = 0; j < pz; j++ ) {
        for ( int k = 0; k < n; k++ ) {
            vzM[k][j] = Vec[k] * z[k][j];
        }
    }
    free(exb);
    free(eta);
    free(Vec);

    double **rvmatrix = vec2matrix(REAL(Rvmatrix_), b, n);
    int **idM = ivec2matrix(INTEGER(ID_), k, n);
    calculate_ST00(
        z, idM, n, pz, b, k, rvmatrix, vzM, REAL(TestR), REAL(TestB), REAL(Pval)
    );

    SET_STRING_ELT(list_names,  0,  mkChar("TestR"));
    SET_STRING_ELT(list_names,  1,  mkChar("TestB"));
    SET_STRING_ELT(list_names,  2,  mkChar("Pvalue"));
    SET_VECTOR_ELT(list,        0,  TestR);
    SET_VECTOR_ELT(list,        1,  TestB);
    SET_VECTOR_ELT(list,        2,  Pval);
    setAttrib(list, R_NamesSymbol,  list_names);

    free(x);
    free(z);
    free(vzMvec);
    free(vzM);
    free(rvmatrix);
    free(idM);
    UNPROTECT(5);
    return list;
}


