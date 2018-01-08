#include <vector>
#include <limits>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <ctime>
#include <tuple>

#define _Inf_ numeric_limits<double>::infinity()

using namespace std;

namespace shield {

class corValPair {
public:
    int p, q; double val;
    vector<corValPair*> hNeighbours;
    vector<corValPair*> vNeighbours;

    corValPair(int _p, int _q, double _val) : p(_p), q(_q), val(_val) {}
};

vector<corValPair*> x;
int m, n, res;
vector<vector<double>> cost;
vector<vector<bool>> isBasis;
vector<double> mu, u;
vector<double> nu, v;
int nr, nc;

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
    isBasis[node->p][node->q] = true;
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
    isBasis[node->p][node->q] = false;
    delete node;
}

void updateuvcore(corValPair* node, bool searchOnRow) {
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
}

void updateuv() {
    v[n-1] = 0;
    for (unsigned int i = 0; i < x.size(); ++i) {
        if (x[i]->q == n-1) {
            updateuvcore(x[i], true);
        }
    }
}

bool isNegativeCost(int i, int j) {
    return cost[i][j] - u[i] - v[j] < -1e-10;
}

int lastindex = 0;
bool mostnegativeusingN(int& nr, int& nc) {
    for (unsigned int __i = lastindex; __i < lastindex+x.size(); ++__i) {
        int _i = __i % x.size();
        corValPair* node = x[_i];
        int p = node->p, q = node->q;
        int s = p % res, t = p / res;
        if (s > 0) {
            if (isNegativeCost(p-1, q)) {
                nr = p-1;
                nc = q;
                lastindex = (__i+1) % x.size();
                return true;
            }
        }
        if (t > 0) {
            if (isNegativeCost(p-res, q)) {
                nr = p-res;
                nc = q;
                lastindex = (__i+1) % x.size();
                return true;
            }
        }
        if (s < res-1) {
            if (isNegativeCost(p+1, q)) {
                nr = p+1;
                nc = q;
                lastindex = (__i+1) % x.size();
                return true;
            }
        }
        if (t < res-1) {
            if (isNegativeCost(p+res, q)) {
                nr = p+res;
                nc = q;
                lastindex = (__i+1) % x.size();
                return true;
            }
        }
        if (s > 0 && t > 0) {
            if (isNegativeCost(p-1-res, q)) {
                nr = p-1-res;
                nc = q;
                lastindex = (__i+1) % x.size();
                return true;
            }
        }
        if (t > 0 && s < res-1) {
            if (isNegativeCost(p+1-res, q)) {
                nr = p+1-res;
                nc = q;
                lastindex = (__i+1) % x.size();
                return true;
            }
        }
        if (s < res-1 && t < res-1) {
            if (isNegativeCost(p+1+res, q)) {
                nr = p+1+res;
                nc = q;
                lastindex = (__i+1) % x.size();
                return true;
            }
        }
        if (t < res-1 && s > 0) {
            if (isNegativeCost(p-1+res, q)) {
                nr = p-1+res;
                nc = q;
                lastindex = (__i+1) % x.size();
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

void printx() {
    printf("================\n");
    for (unsigned int i = 0; i < x.size(); ++i) {
        // printf("(%d, %d)\t%.1e\n", x[i]->p, x[i]->q, x[i]->val);
        // printf("\tv: ");
        printf("(%d, %d)\tv: ", x[i]->p, x[i]->q);
        for (unsigned int j = 0; j < x[i]->vNeighbours.size(); ++j) {
            printf("(%d, %d) ", x[i]->vNeighbours[j]->p, x[i]->vNeighbours[j]->q);
        }
        printf("\n\th: ");
        for (unsigned int j = 0; j < x[i]->hNeighbours.size(); ++j) {
            printf("(%d, %d) ", x[i]->hNeighbours[j]->p, x[i]->hNeighbours[j]->q);
        }
        printf("\n");
    }
    printf("================\n");
}

void core() {
    int itr = 1;
    while (true) {
        if (itr % 1000 == 0) {
            double sum = 0.0;
            for (unsigned int i = 0; i < x.size(); ++i) {
                sum += x[i]->val * cost[x[i]->p][x[i]->q];
            }
            printf("Iteration %d\t%e\n", itr, sum);
        }
        // update the multipliers
        updateuv();
        // returns false if all relative cost are non-negative
        if (!mostnegativeusingN(nr, nc)) {
            break;
        }
        corValPair* newNode = new corValPair(nr, nc, 0.0);
        addNode(newNode);
        corValPair* dontForgetToDelete = new corValPair(nr, nc, _Inf_);
        if (!happySearch(0, false, newNode, _Inf_, dontForgetToDelete)) {
            printf("CTransimplex:core\tLoop detection failed!\n");
            // printx();
            return;
        }
        delete dontForgetToDelete;
        rmNode(globalCor);
        itr++;
        // if (itr > 1000)
        //     break;
    }
    printf("%d iterations in all.\n", itr);
}

vector<tuple<int, int, double>> wrapper_s(const vector<vector<double>> &_cost, const vector<double> &XM, const vector<double> &YM, const vector<tuple<int, int, double>>& x0) {
    m = XM.size();
    res = sqrt(m);
    n = YM.size();
    mu = XM;
    u = vector<double>(m);
    nu = YM;
    v = vector<double>(n);
    cost = _cost;
    isBasis.resize(m);
    for (int i = 0; i < m; ++i) {
        isBasis[i] = vector<bool>(n, false);
    }
    printf("(%d, %d, %zd)\n", m, n, x0.size());

    nr = 0; nc = 0;

    for (auto itr = x0.begin(); itr != x0.end(); ++itr) {
        addNode(new corValPair(get<0>(*itr), get<1>(*itr), get<2>(*itr)));
    }

    core();
    vector<tuple<int, int, double>> r;
    r.reserve(x.size());
    for (auto itr = x.begin(); itr != x.end(); ++itr) {
        r.push_back(make_tuple((*itr)->p, (*itr)->q, (*itr)->val));
        delete *itr;
    }
    x.clear();
    x.resize(0);
    cost.clear();
    cost.resize(0);
    isBasis.clear();
    isBasis.resize(0);
    return r;
}

};