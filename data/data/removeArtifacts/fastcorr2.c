/*==========================================================
 * Jiarui Wang | jwang04@g.harvard.edu
 * Gabriel Kreiman Lab | Boston Children's Hospital
 *
 * fastcorr2.c - Pearson's correlation coefficient
 *
 * Takes two 1xN matrices and outputs a scalar (r)
 *
 * The calling syntax is:
 *
 *		r = fastcorr2(v1, v2)
 *
 * This is a MEX-file for MATLAB.
 *
 *========================================================*/

#include "mex.h"
#include "math.h"

/* The computational routine */
void corr2(double *z,double *x, double *y, mwSize n)
{
    double r,xx[n],xy[n],yy[n],nr=0,dr_1=0,dr_2=0,dr_3=0,dr=0;
    double sum_y=0,sum_yy=0,sum_xy=0,sum_x=0,sum_xx=0;
    mwSize i;
    for (i=0; i<n; i++) {
        xx[i] = x[i]*x[i];
        yy[i] = y[i]*y[i];
    }
    for (i=0; i<n; i++) {
        sum_x += x[i];
        sum_y += y[i];
        sum_xx += xx[i];
        sum_yy += yy[i];
        sum_xy += x[i]*y[i];     
    }
    nr = (n*sum_xy)-(sum_x*sum_y);
    double sum_x2 = sum_x*sum_x;
    double sum_y2 = sum_y*sum_y;
    dr_1 = (n*sum_xx)-sum_x2;
    dr_2 = (n*sum_yy)-sum_y2;
    dr_3 = dr_1*dr_2;
    dr = sqrt(dr_3);
    z[0] = (nr/dr);
}

void corr2old(double *z,double *x, double *y, mwSize n)
{
    double x_sum;
    double y_sum;
    double x2_sum;
    double y2_sum;
    double xy_sum;
    /* Compute all sums */
    mwSize i;
    for (i = 0; i < n; i++) {
        x_sum = x_sum + x[i];
        y_sum = y_sum + y[i];
        x2_sum = x2_sum + x[i]*x[i];
        y2_sum = y2_sum + y[i]*y[i];
        xy_sum = xy_sum + x[i]*y[i];
    }
    
    
    z[0] = (xy_sum-(x_sum*y_sum)/n)/sqrt((x2_sum-(x_sum*x_sum)/n)*(y2_sum-(y_sum*y_sum)/n));
}

/*
 * The gateway function
 *
 * nlhs	Number of output (left-side) arguments, or the size of the plhs array.
 * plhs	Array of output arguments.
 * nrhs	Number of input (right-side) arguments, or the size of the prhs array.
 * prhs	Array of input arguments.
 *
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *V1;              /* 1xN first input matrix */
    double *V2;               /* 1xN second input matrix */
    size_t ncols_V1;                   /* size of first matrix */
    size_t ncols_V2;                   /* size of second matrix */
    double *outMatrix;              /* output matrix */

    /* check for proper number of arguments */
    if(nrhs != 2) {
        mexErrMsgIdAndTxt("MyToolbox:fastcorr2:nrhs","Two inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:fastcorr2:nlhs","One output required.");
    }
    /* make sure the first input argument is type double */
    if( !mxIsDouble(prhs[0]) || 
         mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("MyToolbox:fastcorr2:notDouble","Input matrix must be type double.");
    }
    
    /* make sure the second input argument is type double */
    if( !mxIsDouble(prhs[1]) || 
         mxIsComplex(prhs[1])) {
        mexErrMsgIdAndTxt("MyToolbox:fastcorr2:notDouble","Input matrix must be type double.");
    }
    
    /* check that number of rows in first input argument is 1 */
    if(mxGetN(prhs[0])!=1) {
        mexErrMsgIdAndTxt("MyToolbox:fastcorr2:notRowVector","Input must be a row vector.");
    }
    
    /* check that number of rows in second input argument is 1 */
    if(mxGetN(prhs[1])!=1) {
        mexErrMsgIdAndTxt("MyToolbox:fastcorr2:notRowVector","Input must be a row vector.");
    }
    
    /* create a pointer to the data in first input matrix  */
    V1 = mxGetPr(prhs[0]);

    /* create a pointer to the real data in the second input matrix  */
    V2 = mxGetPr(prhs[1]);

    /* get dimensions of the first input matrix */
    ncols_V1 = mxGetM(prhs[0]);

    /* get dimensions of the second input matrix */
    ncols_V2 = mxGetM(prhs[1]);

    /* Check to make sure the first and second input matrices are right */
    if (ncols_V1 != ncols_V2) {
        mexErrMsgIdAndTxt("MyToolbox:fastcorr2:notSameSize","Inputs must be same size.");
    }

    /*plhs[0] = mxCreateDoubleScalar(corr2(V1,V2,(mwSize)ncols_V1));*/
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    outMatrix = mxGetPr(plhs[0]);
    corr2(outMatrix,V1,V2,(mwSize)ncols_V1);
}
