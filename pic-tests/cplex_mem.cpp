/* watch how many mem cplex uses */
#include <ilcplex/cplex.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

using namespace std;

static void free_and_null (char **ptr) {
    if ( *ptr != NULL ) {
        free (*ptr);
        *ptr = NULL;
    }
}

int main() {
    int res = 32;
    auto m = res*res;
    auto n = res*res;

    vector<vector<double>> cost = vector<vector<double>>(m);
    for (int i = 0; i < m; ++i) {
        int ii = i%res, ij = i/res;
        cost[i].resize(n);
        for (int j = 0; j < n; ++j) {
            int ji = j%res, jj = j/res;
            cost[i][j] = (ii-ji)*(ii-ji) + (ij-jj)*(ij-jj);
        }
    }
    auto mu = vector<double>(m);
    double musum = 0;
    for (int i = 0; i < m; ++i) {
        mu[i] = (rand()+0.0) / RAND_MAX;
        musum += mu[i];
    }
    for (int i = 0; i < m; ++i) {
        mu[i] /= musum;
    }
    auto nu = vector<double>(n);
    double nusum = 0;
    for (int j = 0; j < n; ++j) {
        nu[j] = (rand()+0.0) / RAND_MAX;
        nusum += nu[j];
    }
    for (int j = 0; j < n; ++j) {
        nu[j] /= nusum;
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

    /* cleanup memory */

    free_and_null ((char **) &x);
    free_and_null ((char **) &dj);
    free_and_null ((char **) &pi);
    free_and_null ((char **) &slack);
    CPXNETfreeprob (env, &net);
    CPXcloseCPLEX(&env);
    delete supply, tail, head, obj, ub, lb;

    return 0;
}