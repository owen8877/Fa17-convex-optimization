#include "mex.h"
#include <vector>
#include <limits>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <ctime>
#include <tuple>

#define _Inf_ numeric_limits<double>::infinity()

using namespace std;

class corValPair {
public:
    int p, q; double val;
    vector<corValPair*> hNeighbours;
    vector<corValPair*> vNeighbours;

    corValPair(int _p, int _q, double _val) : p(_p), q(_q), val(_val) {}
};

vector<corValPair*> x;
int m, n;
vector<vector<double>> cost;
vector<double> mu, u;
vector<double> nu, v;
int nr, nc;

void printArray(const vector<double> & array) {
    for (unsigned int i = 0; i < array.size(); ++i) {
        mexPrintf("%.2e ", array[i]);
    }
    mexPrintf("\n");
}

void printx() {
    mexPrintf("================\n");
    for (unsigned int i = 0; i < x.size(); ++i) {
        // mexPrintf("(%d, %d)\t%.1e\n", x[i]->p, x[i]->q, x[i]->val);
        // mexPrintf("\tv: ");
        mexPrintf("(%d, %d)\tv: ", x[i]->p, x[i]->q);
        for (unsigned int j = 0; j < x[i]->vNeighbours.size(); ++j) {
            mexPrintf("(%d, %d) ", x[i]->vNeighbours[j]->p, x[i]->vNeighbours[j]->q);
        }
        mexPrintf("\n\th: ");
        for (unsigned int j = 0; j < x[i]->hNeighbours.size(); ++j) {
            mexPrintf("(%d, %d) ", x[i]->hNeighbours[j]->p, x[i]->hNeighbours[j]->q);
        }
        mexPrintf("\n");
    }
    mexPrintf("================\n");
}

void setCostToInf(vector<vector<double>> & bakcost, bool direction, int index) {
    if (direction) {
        for (int j = 0; j < n; ++j) {
            bakcost[index][j] = _Inf_;
        }
    } else {
        for (int i = 0; i < m; ++i) {
            bakcost[i][index] = _Inf_;
        }
    }
}

int argmin(const vector<double> & array) {
    int where = 0; double value = _Inf_;
    for (unsigned int i = 0; i < array.size(); ++i) {
        if (array[i] < value) {
            where = i;
            value = array[i];
        }
    }
    return where;
}

void addNode(corValPair* node) {
    for (unsigned int i = 0; i < x.size(); ++i) {
        corValPair* neighbour = x[i];
        if (neighbour->p == node->p) {
            neighbour->hNeighbours.push_back(node);
            node->hNeighbours.push_back(neighbour);
        } else if (neighbour->q == node->q) {
            neighbour->vNeighbours.push_back(node);
            node->vNeighbours.push_back(neighbour);
        }
    }
    x.push_back(node);
}

void rmNode(corValPair* node) {
    for (unsigned int i = 0; i < node->vNeighbours.size(); ++i) {
        corValPair* neighbour = node->vNeighbours[i];
        neighbour->vNeighbours.erase(remove(neighbour->vNeighbours.begin(), neighbour->vNeighbours.end(), node), neighbour->vNeighbours.end());
    }
    for (unsigned int j = 0; j < node->hNeighbours.size(); ++j) {
        corValPair* neighbour = node->hNeighbours[j];
        neighbour->hNeighbours.erase(remove(neighbour->hNeighbours.begin(), neighbour->hNeighbours.end(), node), neighbour->hNeighbours.end());
    }
    x.erase(remove(x.begin(), x.end(), node), x.end());
    delete node;
}

void initialSolution() {
    vector<double> workingmu = mu;
    vector<double> workingnu = nu;
    bool searchOnRow = true; // search on some row
    vector<vector<double>> bakcost = cost;
    vector<double> workingcost = bakcost[0];
    int currentRow = 0; int currentCol = 0;
    for (int l = 0; l < m+n-1; ++l) {
        if (searchOnRow) {
            currentCol = argmin(workingcost);
            if (workingnu[currentCol] > workingmu[currentRow]) {
                addNode(new corValPair(currentRow, currentCol, workingmu[currentRow]));
                workingnu[currentCol] = workingnu[currentCol] - workingmu[currentRow];
                workingmu[currentRow] = 0;
                searchOnRow = false;
                setCostToInf(bakcost, true, currentRow);
                workingcost.resize(m);
                for (int i = 0; i < m; ++i) {
                    workingcost[i] = bakcost[i][currentCol];
                }
            } else {
                addNode(new corValPair(currentRow, currentCol, workingnu[currentCol]));
                workingmu[currentRow] = workingmu[currentRow] - workingnu[currentCol];
                workingnu[currentCol] = 0;
                workingcost[currentCol] = _Inf_;
                setCostToInf(bakcost, false, currentCol);
            }
        } else {
            currentRow = argmin(workingcost);
            if (workingmu[currentRow] > workingnu[currentCol]) {
                addNode(new corValPair(currentRow, currentCol, workingnu[currentCol]));
                workingmu[currentRow] = workingmu[currentRow] - workingnu[currentCol];
                workingnu[currentCol] = 0;
                searchOnRow = true;
                setCostToInf(bakcost, false, currentCol);
                workingcost = bakcost[currentRow];
            } else {
                addNode(new corValPair(currentRow, currentCol, workingmu[currentRow]));
                workingnu[currentCol] = workingnu[currentCol] - workingmu[currentRow];
                workingmu[currentRow] = 0;
                workingcost[currentRow] = _Inf_;
                setCostToInf(bakcost, true, currentRow);
            }
        }
    }
}

void updateuvcore(corValPair* node, bool searchOnRow) {
    // mexPrintf("(%d, %d)\n", node->p, node->q);
    if (searchOnRow) {
        u[node->p] = cost[node->p][node->q] - v[node->q];
        for (unsigned int j = 0; j < node->hNeighbours.size(); ++j) {
            updateuvcore(node->hNeighbours[j], false);
        }
    } else {
        v[node->q] = cost[node->p][node->q] - u[node->p];
        for (unsigned int i = 0; i < node->vNeighbours.size(); ++i) {
            updateuvcore(node->vNeighbours[i], true);
        }
    }
    // mexPrintf("Finish.\n");
}

void updateuv() {
    v[n-1] = 0;
    for (unsigned int i = 0; i < x.size(); ++i) {
        if (x[i]->q == n-1) {
            updateuvcore(x[i], true);
        }
    }
}

int lastrow = 0;
bool mostnegative(int& nr, int& nc) {
    for (int _i = lastrow; _i < m+lastrow ; ++_i) {
        int i = _i % m;
        for (int j = 0; j < n; ++j) {
            double tmp = cost[i][j] - u[i] - v[j];
            if (tmp < -1e-15) {
                nr = i;
                nc = j;
                lastrow = (_i+1)%m;
                return true;
            }
        }
    }
    return false;
}

double globalMin;
corValPair* globalCor;

bool happySearch(int depth, bool searchOnRow, corValPair* node, double currentMin, corValPair* wasmin) {
    if ((node->p == nr) && (node->q == nc) && depth > 0) {
        // we have finnaly found the loop!
        globalMin = currentMin;
        globalCor = wasmin;
        return true;
    }
    if (searchOnRow) {
        if (currentMin > node->val) {
            currentMin = node->val;
            wasmin = node;
        }
        for (unsigned int j = 0; j < node->hNeighbours.size(); ++j) {
            if (happySearch(depth+1, false, node->hNeighbours[j], currentMin, wasmin)) {
                node->val -= globalMin;
                return true;
            }
        }
    } else {
        for (unsigned int i = 0; i < node->vNeighbours.size(); ++i) {
            if (happySearch(depth+1, true, node->vNeighbours[i], currentMin, wasmin)) {
                node->val += globalMin;
                return true;
            }
        }
    }
    return false;
}

void core() {
    // clock_t begin = clock();

    initialSolution();

    // clock_t end = clock();
    // double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    // mexPrintf("%.4fs have been used on inital.\n", elapsed_secs);
    // begin = clock();

    int itr = 1;
    while (true) {
        if (itr % 2000 == 0) {
            double sum = 0.0;
            for (unsigned int i = 0; i < x.size(); ++i) {
                sum += x[i]->val * cost[x[i]->p][x[i]->q];
            }
            mexPrintf("Iteration %d\t%e\n", itr, sum);
        }
        // update the multipliers
        updateuv();
        // returns false if all relative cost are non-negative
        if (!mostnegative(nr, nc)) {
            break;
        }
        corValPair* newNode = new corValPair(nr, nc, 0.0);
        addNode(newNode);
        corValPair* dontForgetToDelete = new corValPair(nr, nc, _Inf_);
        if (!happySearch(0, false, newNode, _Inf_, dontForgetToDelete)) {
            mexPrintf("CTransimplex:core\tLoop detection failed!\n");
            printx();
            return;
        }
        delete dontForgetToDelete;
        rmNode(globalCor);
        itr++;
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs != 5) {
        mexErrMsgIdAndTxt("CTransimplex:gateway:nrhs", "Need five inputs!");
        return;
    }
    if (nlhs != 2) {
        mexErrMsgIdAndTxt("CTransimplex:gateway:nlhs", "Need two outputs!");
        return;
    }
    
    // read in data
    // cost matrix
    double* costRe;
    m = mxGetM(prhs[1]);
    n = mxGetN(prhs[1]);
    cost = vector<vector<double>>(m);
    costRe = mxGetPr(prhs[1]);
    // mu
    double* muRe = mxGetPr(prhs[2]);
    if (m != mxGetM(prhs[2])) {
        mexErrMsgIdAndTxt("CTransimplex:gateway:mu", "Mu dimension not match with cost matrix!");
        return;
    }
    mu = vector<double>(m);
    u = vector<double>(m);
    // nu
    double* nuRe = mxGetPr(prhs[3]);
    if (n != mxGetM(prhs[3])) {
        mexErrMsgIdAndTxt("CTransimplex:gateway:nu", "Nu dimension not match with cost matrix!");
        return;
    }
    nu = vector<double>(n);
    v = vector<double>(n);

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

    // some magic here
    core();
    
    plhs[0] = mxCreateNumericMatrix(2, x.size(), mxINT32_CLASS, mxREAL);
    int32_t* val1 = (int32_t*) mxGetData(plhs[0]);
    for (unsigned int i = 0; i < x.size(); ++i) {
        val1[2*i] = x[i]->p;
        val1[2*i+1] = x[i]->q;
    }
    plhs[1] = mxCreateDoubleMatrix(1, x.size(), mxREAL);
    double* val2 = mxGetPr(plhs[1]);
    for (unsigned int i = 0; i < x.size(); ++i) {
        val2[i] = x[i]->val;
        delete x[i];
    }

    return;
}