#include "mex.h"
#include "./imm.h"
#include "./imm.cpp"
#include <cstdio>

//#pragma comment(lib,"usingImmAlg.lib")

void mexFunction(int nlhs,mxArray* plhs[],int nrhs, mxArray* prhs[]){
    double *input = mxGetPr(prhs[0]);
    imm immDemo;
   // immDemo.setMeasureData(1,2);
}