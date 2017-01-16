#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <inttypes.h>
#include <string.h>
#include <mex.h>
#include "matrix.h"

typedef UINT32_T uint32;


double preSum(double *data, int maxIdx)
{
  double sumValue=0;
  for(int z=0; z < maxIdx; ++z)
  {
      sumValue = sumValue + data[z];
  }
  return sumValue;
}

double postSum(double *data, int minIdx, mwSize numElements)
{
    double sumValue = 0;
    for(int z=minIdx; z < numElements; ++z)
        sumValue = sumValue + data[z];
    return sumValue;
}

int findIndex(double *timeValues, double whichTime, mwSize numElements)
{
    int z = 0;
    for(z=0; z < numElements; ++z)
    {
        if(timeValues[z] >= whichTime){
            break;
        }
    }
    return z;
}

double fastSum(double *data, mwSize numElements)
{
    double sum =0;
    for(int j=0;j < numElements; j++)
        sum = sum + data[j];
    return sum;
}

double fastMean(double *data, mwSize numElements)
{
    double sum = 0;
    for(int j=0;j < numElements; j++)
        sum = sum + data[j];
    return sum/(numElements*1.0);
}

void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[] )
{
    // First 
    // Param 1, RT
    // Param 2, Spikes
    // Param 3, Lambda parameter
    // Param 4, beta parameter
    // Param 5, binned time values
    
    double *RT;
    double *binnedSpikes;
    double *timeValues;
    double lambda, beta;
    double meanRT;
    
    double tCurr;
    double *currTrial;
    
    double *finalSums, *preSpikes, *postSpikes;
    
    double preErr, postErr;
    
    int numPreSpikes = 0;
    int numPostSpikes = 0;
    
    int totalPreSpikes;
    int totalPostSpikes = 0;
    
    mwSize nTrials, nTimePoints;
   
    
    RT = mxGetPr(prhs[0]);
    binnedSpikes = mxGetPr(prhs[1]);
    lambda = mxGetScalar(prhs[2]);
    beta = mxGetScalar(prhs[3]);
    timeValues = mxGetPr(prhs[4]);
   
    nTrials = mxGetM(prhs[1]);
    nTimePoints = mxGetN(prhs[1]);
    meanRT = fastMean(RT, mxGetN(prhs[0]));
    
    // mexPrintf("%3.2f",meanRT);
    currTrial = (double*)mxMalloc(nTimePoints*sizeof(double));
    
    // mexPrintf("\n num Trials: %d, num Time Points:%d", nTrials, nTimePoints);
    
    plhs[0] = mxCreateNumericMatrix(1,2,mxDOUBLE_CLASS, mxREAL);
    plhs[1] = mxCreateNumericMatrix(1,nTrials*nTimePoints,mxDOUBLE_CLASS, mxREAL);
    plhs[2] = mxCreateNumericMatrix(1,nTrials*nTimePoints,mxDOUBLE_CLASS, mxREAL);
    
   
    finalSums = mxGetPr(plhs[0]);
    preSpikes = mxGetPr(plhs[1]);
    postSpikes = mxGetPr(plhs[2]);
    
    totalPreSpikes = 0;
    
    // For each trial, calculate the pre sum and and the post sum
    
    for(int k=0; k < nTrials;++k)
    {
        // mexPrintf("\n %3.2f", RT[k]);
        
        tCurr = lambda*meanRT/1000.0 + beta*(RT[k] - meanRT)/1000.0;
        // mexPrintf("\n %3.4f, %d \n", tCurr, nTimePoints);
        
        int preIdx = findIndex(timeValues,tCurr, nTimePoints);
        int postIdx = findIndex(timeValues, (RT[k]+100)/1000, nTimePoints);
        // mexPrintf("%d", postIdx);
        
        for(int timePts=0; timePts < nTimePoints; ++ timePts){
            currTrial[timePts]= binnedSpikes[(timePts)*nTrials + k];

        }
   
        // Keep copying until pre timepoints;
        for(int z=0; z < preIdx; ++z)
        {
            preSpikes[z+totalPreSpikes] = currTrial[z];
            // mexPrintf("%3.2f, ", preSpikes[z+totalPreSpikes]);
        }
        int numV = 0;
        int z = preIdx;    
        for(; z < postIdx;++z, ++numV)
        {
            // postSpikes[totalPreSpikes+z:]
            postSpikes[totalPostSpikes+numV] = currTrial[numV + preIdx];
        }

        // mexPrintf("...%d", preIdx);    
        totalPreSpikes = totalPreSpikes + preIdx;
        totalPostSpikes = totalPostSpikes + numV;
        // mexPrintf("...%d, %d",totalPreSpikes, totalPostSpikes);
        
    }
    
    
    preSpikes[totalPreSpikes+1] = 99;
    postSpikes[totalPostSpikes+1] = 99;
    
    finalSums[0] = (double) totalPreSpikes;
    finalSums[1] = (double) totalPostSpikes;

    // totalPreSpikes = numPreSpikes;
    
}