// =========================================================================
// -------------------------------------------------------------------------
//
//        Estmation of Mutual Information Using Adaptive Histogram
//
// -------------------------------------------------------------------------
// =========================================================================
//    Version : Last modified: 25 March 2016
//    Since   : 10 November 2015
//    Author  : Qi Xu
//              School of Computer Science, University of Birmingham
//    Contact : cxyandyxu@gmail.com
// =========================================================================
// If you use this code then please cite the following
//    Reference:
//      Qi Xu, M.L. Sanyang, A.Kaban. Large Scale Continuous EDA Using Mutual Information. 
//      IEEE Congress on Evolutionary Computation 2016 (CEC-2016)
// =========================================================================
//     This program is free software: you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
//
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
//
//     You should have received a copy of the GNU General Public License
//     along with this program.  If not, see <http://www.gnu.org/licenses/>.
// =========================================================================
// NOTICE:
//    This code requires compilers supporting C++11 features. 
//    
//    For Windows users, MS VC++ 2013 or higher is recommended.
//    For Linux users, gcc 4.7.4 or higher is recommended.
//    
//    This code should be compiled using 'mex' in Matlab.
//    Code optimisation options of compilers are highly recommended.
//
// =========================================================================
// INPUT:
//    N-by-D Population.
//    
//--------------------------------------------------------------------------
// OUTPUT:
//    D-by-D mutual information matrix, with diagonal equals 0.
//    
// =========================================================================



#include "mex.h"
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <stack>
#include <thread>
using std::stack;
using std::thread;

/*
calcMI
parameters:
in1: first vector a;
in2: second vector b;
N  : PopSize;
dim: equals 2 by default and DO NOT change it unless you know what happens.
*/

#define NMAX 2048
#define DMAX 1024
#define NDIM 2
#define DDIM 4
#define DIM2 4

int index[NMAX][DMAX]{ 0 };
double outx[NMAX][DMAX]{ 0.0 }; // dxd matrix.
double inx[DMAX][NMAX]{ 0.0 };
int B = 0, M, N;
bool fin1 = false, fin2 = false, fin3 = false, fin4 = false;

int RandPartition(double data[], int lo, int hi, int k) {
    double v = data[lo];
    int vi = index[lo][k];
    while (lo < hi) {
        while (lo < hi && data[hi] >= v)
            hi--;
        data[lo] = data[hi];
        index[lo][k] = index[hi][k];
        while (lo < hi && data[lo] <= v)
            lo++;
        data[hi] = data[lo];
        index[hi][k] = index[lo][k];
    }
    data[lo] = v;
    index[lo][k] = vi;
    return lo;
}

void QuickSort(double data[], int lo, int hi, int k) {
    stack<int> st;
    int key;
    do {
        while (lo < hi) {
            key = RandPartition(data, lo, hi, k);
            if ((key - lo) < (key - hi)) {
                st.push(key + 1);
                st.push(hi);
                hi = key - 1;
            }
            else
            {
                st.push(lo);
                st.push(key - 1);
                lo = key + 1;
            }
        }
        if (st.empty()) return;
        hi = st.top();
        st.pop();
        lo = st.top();
        st.pop();
    } while (1);
}

double calcMI(int in1, int in2, int Num, int dim) {
    //[Num][dim]
    /*cout << "in1: " << in1 << endl;
    cout << "in2: " << in2 << endl;
    cout << "M: " << M << endl;
    cout << "dim:" << dim << endl;*/
    int xindex[2] { 0 }; // Index for the candidate Index.
    int ydat[NMAX][NDIM]{ 0 };
    xindex[0] = in1;
    xindex[1] = in2;
    int i, j, poradi[NMAX] { 0 }; // Num
    for (i = 0; i < Num; i++) {
        poradi[i] = i + 1;
    }
    // Set ydat
    for (i = 0; i < Num; i++) {
        for (j = 0; j < dim; j++) {
            ydat[index[i][xindex[j]] - 1][j] = i + 1;
        }
    }

    int ddim = 1 << dim, dim2 = 2 * dim;
    double xcor = 0.0;
    int npar = 1, poc[4 * DDIM] { 1 }, kon[4 * DDIM] { Num };
    int NN[DDIM] { 0 };
    int marg[8 * DDIM][DIM2] { 0 };
    int amarg[DDIM][DIM2] { 0 };
    // Set marg
    marg[0][0] = marg[0][1] = 1;
    marg[0][2] = marg[0][3] = Num;
    // Set Imm, chi2, run
    int Imm[4][2] { { 0,0 },{ 0,1 },{ 1,0 },{ 1,1 } };
    double chi2[5] { 0.0,7.810,13.90,25.0,42.0 };
    int run = 0;
    int apor[NMAX] { 0 };
    int J[NMAX][NDIM] { 0 };
    int I[NMAX][DDIM] { 0 };
    int Nex = 0;
    int fi[NMAX][DDIM] { 0 };
    int Nxx = 0;
    int d, k, ind;
    int ave[2] { 0 };
    while (npar > 0) {
        run++;
        int apoc = poc[npar - 1], akon = kon[npar - 1];
        for (i = 0; i < akon - apoc + 1; i++) {
            apor[i] = poradi[i + apoc - 1];
        }
        Nex = akon - apoc + 1;  // length(apor)
        ave[0] = (marg[npar - 1][0] + marg[npar - 1][2]) / 2;
        ave[1] = (marg[npar - 1][1] + marg[npar - 1][3]) / 2;
        // Set J
        for (i = 0; i < Nex; i++) {
            for (j = 0; j < dim; j++) {
                if (ydat[apor[i] - 1][j] <= ave[j])
                    J[i][j] = 1;
                else
                    J[i][j] = 0;
            }
        }
        // amarg
        for (i = 0; i < ddim; i++) {
            for (j = 0; j < ddim; j++) {
                amarg[i][j] = marg[npar - 1][j];
            }
        }
        // Set I
        for (d = 0; d < ddim; d++) {
            for (i = 0; i < Nex; i++) {
                I[i][d] = 1;
            }
            for (k = 0; k < dim; k++) {
                if (Imm[d][k]) {
                    for (i = 0; i < Nex; i++) {
                        I[i][d] &= ~J[i][k];
                        amarg[d][k] = ave[k] + 1;
                    }
                }
                else {
                    for (i = 0; i < Nex; i++) {
                        I[i][d] &= J[i][k];
                        amarg[d][k + dim] = ave[k];
                    }
                }
            }
        }
        // NN=sum(I)
        for (i = 0; i < ddim; i++) {
            NN[i] = 0;
            for (j = 0; j < Nex; j++) {
                if (I[j][i] > 0) NN[i]++;
            }
        }
        // Compute tst
        double tmp = 0.0;
        double sum = 0.0;
        for (i = 0; i < ddim; i++) {
            sum += (((double)NN[i] - (double)Nex / (double)ddim)*((double)NN[i] - (double)Nex / (double)ddim));
        }
        double tst = (double)ddim*sum / (double)Nex;
        //cout << tst << " ";
        if ((tst > chi2[dim - 1]) || (run == 1)) {
            npar--;
            for (ind = 0; ind < ddim; ind++) {
                if (NN[ind] > ddim) {
                    npar++;
                    akon = apoc + NN[ind] - 1;
                    poc[npar - 1] = apoc;
                    kon[npar - 1] = akon;
                    for (i = 0; i < dim2; i++) {
                        marg[npar - 1][i] = amarg[ind][i];
                    }
                    // poradi(apoc:akon)=apor(find(I(:,ind))); 
                    for (i = 0, j = 0; i < Nex; i++) {
                        if (I[i][ind] != 0) {
                            fi[j][ind] = i + 1;
                            j++;
                        }
                    }
                    for (i = apoc - 1, j = 0; i < akon; i++, j++) {
                        poradi[i] = apor[fi[j][ind] - 1];
                    }
                    apoc = akon + 1;
                }
                else {
                    if (NN[ind] > 0) {
                        Nxx = 1;
                        for (i = 0; i < dim; i++) {
                            Nxx = Nxx * (amarg[ind][dim + i] - amarg[ind][i] + 1);
                        }
                        xcor = xcor + (double)NN[ind] * log((double)NN[ind] / (double)Nxx);
                    }
                }
            }
        }
        else {
            Nxx = 1;
            for (i = 0; i < dim; i++) {
                Nxx = Nxx * (marg[npar - 1][dim + i] - marg[npar - 1][i] + 1);
            }
            xcor = xcor + (double)Nex*log((double)Nex / (double)Nxx);
            npar--;
        }
    }
    // Return - yes it's finished!
    return (xcor / (double)Num + (double)(dim - 1)*log((double)Num));
    //cout << "in1 = " << in1 << ", in2 = " << in2 << std::endl;
}

void calcPart1() {
    int i, j;
    for (i = 0; i < B; i++) {
        for (j = i + 1; j < B; j++) {
            //if (mt1.try_lock()) {
                outx[i][j] = outx[j][i] = calcMI(i, j, M, NDIM);
             //  mt1.unlock();
            //}
        }
        outx[i][i] = 0.0;
    }
    fin1 = true;
    return;
}
void calcPart2() {
    int i, j;
    for (i = 0; i < B; i++) {
        for (j = B; j < B + i; j++) {
            outx[i][j] = outx[j][i] = calcMI(i, j, M, NDIM);
        }
        outx[i][i] = 0.0;
    }
    fin2 = true;
}
void calcPart3() {
    int i, j;
    for (i = 0; i < B; i++) {
        for (j = B + i; j < N; j++) {
            outx[i][j] = outx[j][i] = calcMI(i, j, M, NDIM);
        }
        outx[i][i] = 0.0;
    }
    fin3 = true;
}
void calcPart4() {
    int i, j;
    for (i = B; i < N; i++) {
        for (j = i + 1; j < N; j++) {
            outx[i][j] = outx[j][i] = calcMI(i, j, M, NDIM);
        }
        outx[i][i] = 0.0;
    }
    fin4 = true;
}

/*
Gateway function - parameters:
nlhs: Number of output (left-side) arguments, or the size of the plhs array.
plhs: Array of output arguments.
nrhs: Number of input (right-side) arguments, or the size of the prhs array.
prhs: Array of input arguments.
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *const prhs[]) {
    // Validate parameters (Just by following the guide)
//    if (nrhs != 1) { mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "One input required."); }
//    if (nlhs != 1) { mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs", "One output required."); }

    // Read input
    double *inm = mxGetPr(prhs[0]);
    double *in1 = NULL;
    double *in2 = NULL;
    double t;
    M = mxGetM(prhs[0]); // PopSize
    N = mxGetN(prhs[0]); // Dimension
//    if (M > NMAX) { mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "PopSize too large."); }
//    if (N > DMAX) { mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "Dimension too large."); }
    int i, j, k, mi, ti;

    B = N / 2;

    // Store input data
    for (i = 0; i < N; i++) {
        in1 = inm + i*M;
        memcpy(inx[i], in1, sizeof(double)*M);
    }

    // Initialise index matrix.
    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {
            index[i][j] = i + 1;
        }
    }

    // Sort the whole input matrix, as well as index matrix.
    for (k = 0; k < N; k++) {
        QuickSort(inx[k], 0, M - 1, k);
    }

    // Compute the MI matrix - MultiThread.

    thread t1(calcPart1);
    thread t2(calcPart2);
    thread t3(calcPart3);
    thread t4(calcPart4);
    if(t1.joinable())t1.join();
    if(t2.joinable())t2.join();
    if(t3.joinable())t3.join();
    if(t4.joinable())t4.join();
    //    for(i=0;i<N;i++){
    //        outx[i][i] = 0.0;
    //        for(j=i+1;j<N;j++){
    //            outx[i][j] = outx[j][i] = calcMI(i,j,M,NDIM);
    //        }
    //    }
    // Prepare output
    //while (fin1 == false || fin2 == false || fin3 == false || fin4 == false);
    
    
    plhs[0] = mxCreateDoubleMatrix(N, N, mxREAL);
    double *outMatrix = mxGetPr(plhs[0]);
    for (i = 0; i < N; i++) {
        memcpy(&outMatrix[i*N], outx[i], sizeof(double)*N);
    }
}
// End of programme.