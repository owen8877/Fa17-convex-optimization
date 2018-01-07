#include <ilcplex/cplex.h>
#include <stdlib.h>
#include <string.h>

static void free_and_null (char **ptr) {
    if ( *ptr != NULL ) {
        free (*ptr);
        *ptr = NULL;
    }
}

int main() {
    int     narcs;
    int     nnodes;
    int     solstat;
    double  objval;
    double  *x      = NULL;
    double  *pi     = NULL;
    double  *slack = NULL;
    double  *dj     = NULL;
    int     status;

    CPXENVptr env = CPXopenCPLEX(&status);
    CPXsetintparam(env, CPXPARAM_ScreenOutput, CPX_ON);
    CPXNETptr net = CPXNETcreateprob(env, &status, "netex1");

#  define NNODES  8
#  define NARCS  14
#  define inf     CPX_INFBOUND
    double  supply[NNODES] = {20.0, 0.0, 0.0, -15.0, 5.0, 0.0, 0.0, -10.0};
    int     tail[NARCS] = {    0,     1,     2,     3,     6,     5,     4,
                                     4,     2,     3,     3,     5,     5,     1};
    int     head[NARCS] = {    1,     2,     3,     6,     5,     7,     7,
                                     1,     1,     4,     5,     3,     4,     5};

    double  obj [NARCS] = { 3.0,  3.0,  4.0,  3.0,  5.0,  6.0,  7.0,
                                  4.0,  2.0,  6.0,  5.0,  4.0,  3.0,  6.0};
    double  ub  [NARCS] = {24.0, 25.0, 12.0, 10.0,  9.0,  inf, 20.0,
                                 10.0,  5.0, 15.0, 10.0, 11.0,  6.0,  inf};
    double  lb  [NARCS] = {18.0,  0.0, 12.0,  0.0,  0.0, -inf,  0.0,
                                  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0};
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
    CPXNETsolution (env, net, &solstat, &objval, x, pi, slack, dj);
    printf ("\nSolution status = %d\n", solstat);
    printf ("Solution value  = %f\n\n", objval);

    for (int i = 0; i < nnodes; i++) {
        printf ("Node %2d:  Slack = %10f  Pi = %10f\n", i, slack[i], pi[i]);
    }

    for (int j = 0; j < narcs; j++) {
        printf ("Arc  %2d:  Value = %10f  Reduced cost = %10f\n",
                  j, x[j], dj[j]);
    }

    free_and_null ((char **) &x);
    free_and_null ((char **) &dj);
    free_and_null ((char **) &pi);
    free_and_null ((char **) &slack);
    CPXNETfreeprob (env, &net);
    CPXcloseCPLEX(&env);
      
    return 0;
}