#include <ilcplex/cplex.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include "mex.h"
using namespace std;

static void free_and_null (char **ptr) {
    if ( *ptr != NULL ) {
        free (*ptr);
        *ptr = NULL;
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs != 5) {
        mexErrMsgIdAndTxt("CPLEX:gateway:nrhs", "Need five inputs!");
        return;
    }
    if (nlhs != 2) {
        mexErrMsgIdAndTxt("CPLEX:gateway:nlhs", "Need two outputs!");
        return;
    }
    
    // read in data
    // cost matrix
    auto m = mxGetM(prhs[1]);
    auto n = mxGetN(prhs[1]);
    auto cost = vector<vector<double>>(m);
    auto costRe = mxGetPr(prhs[1]);
    // mu
    double* muRe = mxGetPr(prhs[2]);
    if (m != mxGetM(prhs[2])) {
        mexErrMsgIdAndTxt("CPLEX:gateway:mu", "Mu dimension not match with cost matrix!");
        return;
    }
    auto mu = vector<double>(m);
    auto u = vector<double>(m);
    // nu
    double* nuRe = mxGetPr(prhs[3]);
    if (n != mxGetM(prhs[3])) {
        mexErrMsgIdAndTxt("CPLEX:gateway:nu", "Nu dimension not match with cost matrix!");
        return;
    }
    auto nu = vector<double>(n);
    auto v = vector<double>(n);

    for (int i = 0; i < m; ++i) {
        cost[i].resize(n);
        for (int j = 0; j < n; ++j) {
            cost[i][j] = costRe[m*j+i];
        }
    }
    
    for (int i = 0; i < m; ++i) {
        mu[i] = muRe[i];
    }
    
    for (int j = 0; j < n; ++j) {
        nu[j] = nuRe[j];
    }

    /* calling CPLEX */

    int narcs;
    int nnodes;
    int solstat;
    double objval;
    double *x = NULL;
    double *pi = NULL;
    double *slack = NULL;
    double *dj = NULL;
    int status;

    CPXENVptr env = CPXopenCPLEX(&status);
    CPXsetintparam(env, CPXPARAM_ScreenOutput, CPX_ON);
    CPXNETptr net = CPXNETcreateprob(env, &status, "netex1");

    const int NNODES = m+n;
    const int NARCS = m*n;
#   define inf CPX_INFBOUND
    double *supply = new double[NNODES];
    for (int i = 0; i < m; ++i) supply[i] = mu[i];
    for (int j = 0; j < n; ++j) supply[m+j] = -nu[j];

    int *tail = new int[NARCS];
    int *head = new int[NARCS];
    double *obj = new double[NARCS];
    double *ub = new double[NARCS];
    double *lb = new double[NARCS];
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            tail[i*n+j] = i;
            head[i*n+j] = m+j;
            obj[i*n+j] = cost[i][j];
            ub[i*n+j] = inf;
            lb[i*n+j] = 0;
        }
    }

    CPXNETdelnodes(env, net, 0, CPXNETgetnumnodes(env, net)-1);
    CPXNETchgobjsen(env, net, CPX_MIN);
    CPXNETaddnodes(env, net, NNODES, supply, NULL);
    CPXNETaddarcs(env, net, NARCS, tail, head, lb, ub, obj, NULL);
    CPXNETprimopt(env, net);
    narcs  = CPXNETgetnumarcs  (env, net);
    nnodes = CPXNETgetnumnodes (env, net);
    x      = (double *) malloc (narcs  * sizeof (double));
    dj     = (double *) malloc (narcs  * sizeof (double));
    pi     = (double *) malloc (nnodes * sizeof (double));
    slack  = (double *) malloc (nnodes * sizeof (double));
    CPXNETsolution(env, net, &solstat, &objval, x, pi, slack, dj);
    
    vector<int> indexes;
    vector<double> solution;
    indexes.reserve(m+n-1);
    solution.reserve(m+n-1);

    for (int j = 0; j < narcs; j++) {
        // printf("x[j] is %f.\n", x[j]);
        if (x[j] > 1e-15) {
            indexes.push_back(j);
            solution.push_back(x[j]);
        }
    }

    /* end of calling */
    
    plhs[0] = mxCreateNumericMatrix(2, indexes.size(), mxINT32_CLASS, mxREAL);
    int32_t* val1 = (int32_t*) mxGetData(plhs[0]);
    for (unsigned int i = 0; i < indexes.size(); ++i) {
        val1[2*i] = indexes[i]/n;
        val1[2*i+1] = indexes[i]%n;
    }
    plhs[1] = mxCreateDoubleMatrix(1, indexes.size(), mxREAL);
    double* val2 = mxGetPr(plhs[1]);
    for (unsigned int i = 0; i < indexes.size(); ++i) {
        val2[i] = solution[i];
    }

    /* cleanup memory */

    free_and_null ((char **) &x);
    free_and_null ((char **) &dj);
    free_and_null ((char **) &pi);
    free_and_null ((char **) &slack);
    CPXNETfreeprob (env, &net);
    CPXcloseCPLEX(&env);
    delete supply, tail, head, obj, ub, lb;

    return;
}