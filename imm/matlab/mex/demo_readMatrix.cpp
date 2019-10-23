#include "mex.h"
#include<cstdio>
//read matrix from matlab
void mexFunction(int nlhs,mxArray* plhs[],int nrhs,mxArray* prhs[]){
	double *input=mxGetPr(prhs[0]);//传入的第一个参数 
	size_t m=mxGetM(prhs[0]); 
	size_t n=mxGetN(prhs[0]);
	printf("the row of matrix is %d\n",m);
	printf("the column of matrix is %d\n",n);
	
	size_t row=mxGetScalar(prhs[1]);
	size_t col=mxGetScalar(prhs[2]);
	
	printf("the data of row %d column %d is:%f\n",row,col,*(input+(col-1)*m+row-1));
	
	//store a matrix data 
	plhs[0] = mxCreateDoubleMatrix( (mwSize)m, (mwSize)n, mxREAL);
	double *c=mxGetPr(plhs[0]);
	c[0]=5;
	c[1]=6;
	c[10]=33;
    
    int mm = 2;
    int nn = 2;
    plhs[1] = mxCreateDoubleMatrix((mwSize)mm,(mwSize)nn,mxREAL);
    double* cc = mxGetPr(plhs[1]);
    cc[0] = 1;
    cc[3] = 2;
    
    
}