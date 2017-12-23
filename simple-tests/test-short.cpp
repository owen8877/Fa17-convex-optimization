#include <vector>
#include <limits>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <ctime>
#include <tuple>

using namespace std;

// cycle per percent
const double __p__ = 0.05;
// discover items
const int __k__ = 30;
// short list length
const int __s__ = 40;

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
        printf("%.2e ", array[i]);
    }
    printf("\n");
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

void setCostToInf(vector<vector<double>> & bakcost, bool direction, int index) {
    if (direction) {
        for (int j = 0; j < n; ++j) {
            bakcost[index][j] = numeric_limits<double>::infinity();
        }
    } else {
        for (int i = 0; i < m; ++i) {
            bakcost[i][index] = numeric_limits<double>::infinity();
        }
    }
}

int argmin(const vector<double> & array) {
    int where = 0; double value = numeric_limits<double>::infinity();
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
                workingcost[currentCol] = numeric_limits<double>::infinity();
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
                workingcost[currentRow] = numeric_limits<double>::infinity();
                setCostToInf(bakcost, true, currentRow);
            }
        }
    }
}

void updateuvcore(corValPair* node, bool searchOnRow) {
    // printf("(%d, %d)\n", node->p, node->q);
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
    // printf("Finish.\n");
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

int lastrowS = 0;
vector<vector<int>> shortList;

void buildShortList() {
    shortList.resize(m);
    for (int i = 0; i < m; ++i) {
        shortList[i].resize(__s__);
        vector<bool> feasible(n, true);
        for (int l = 0; l < __s__; ++l) {
            int pos; double val = numeric_limits<double>::infinity();
            for (int j = 0; j < n; ++j) {
                if (cost[i][j] < val && feasible[j]) {
                    pos = j;
                    val = cost[i][j];
                }
            }
            shortList[i][l] = pos;
            feasible[pos] = false;
        }
    }
}

vector<tuple<int, int, double>> newlyList; // the last is alaways positive
void mostnegativeShortCore() {
    int rowN = __p__*m;
    for (int _i = lastrowS; _i < rowN+lastrowS ; ++_i) {
        int i = _i % m;
        for (unsigned int _j = 0; _j < shortList[i].size(); ++_j) {
            int j = shortList[i][_j];
            double tmp = cost[i][j] - u[i] - v[j];
            if (tmp < -1e-15) {
                newlyList.push_back(make_tuple(i, j, -tmp));
                if (newlyList.size() > __k__) {
                    //printf("This time discovered %d items before iterating %f percent.\n", __k__, __p__);
                    return;
                }
            }
        }
    }
    return;
}

bool mostnegativeShort(int& nr, int& nc) {
    newlyList.clear();
    mostnegativeShortCore();
    if (newlyList.size() > 0) {
        double val = 0.0;
        for (unsigned int l = 0; l < newlyList.size(); ++l) {
            if (val < get<2>(newlyList[l])) {
                nr = get<0>(newlyList[l]);
                nc = get<1>(newlyList[l]);
            }
        }
        return true;
    } else {
        return false;
    }
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
    clock_t begin = clock();

    initialSolution();

    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    printf("%.4fs have been used on inital.\n", elapsed_secs);
    begin = clock();

    buildShortList();

    end = clock();
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    printf("%.4fs have been used on build short list.\n", elapsed_secs);
    begin = clock();

    int itr = 1;
    bool shortSearch = true;
    while (true) {
        if (itr % 2000 == 0) {
            double sum = 0.0;
            for (unsigned int i = 0; i < x.size(); ++i) {
                sum += x[i]->val * cost[x[i]->p][x[i]->q];
            }
            printf("Iteration %d\t%e\n", itr, sum);
        }
        // update the multipliers
        updateuv();
        // returns false if all relative cost are non-negative
        if (shortSearch) {
            if (!mostnegativeShort(nr, nc)) {
                shortSearch = false;
                printf("Phase 3 ends at iteration %d.\n", itr);
                if (!mostnegative(nr, nc)) {
                    break;
                }
            }
        } else {
            if (!mostnegative(nr, nc)) {
                break;
            }
        }
        corValPair* newNode = new corValPair(nr, nc, 0.0);
        addNode(newNode);
        corValPair* dontForgetToDelete = new corValPair(nr, nc, numeric_limits<double>::infinity());
        if (!happySearch(0, false, newNode, numeric_limits<double>::infinity(), dontForgetToDelete)) {
            printf("CTransimplex:core\tLoop detection failed!\n");
            printx();
            return;
        }
        delete dontForgetToDelete;
        rmNode(globalCor);
        itr++;
    }
}

int main() {
    clock_t begin = clock();
    m = 64*64; n = 64*64;
    cost = vector<vector<double>>(m);
    for (int i = 0; i < m; ++i) {
        cost[i].resize(n);
        for (int j = 0; j < n; ++j) {
            cost[i][j] = (rand()+0.0) / RAND_MAX;
        }
    }
    mu = vector<double>(m);
    u = vector<double>(m);
    double musum = 0;
    for (int i = 0; i < m; ++i) {
        mu[i] = (rand()+0.0) / RAND_MAX;
        musum += mu[i];
    }
    for (int i = 0; i < m; ++i) {
        mu[i] /= musum;
    }
    nu = vector<double>(n);
    v = vector<double>(n);
    double nusum = 0;
    for (int j = 0; j < n; ++j) {
        nu[j] = (rand()+0.0) / RAND_MAX;
        nusum += nu[j];
    }
    for (int j = 0; j < n; ++j) {
        nu[j] /= nusum;
    }

    for (unsigned int i = 0; i < x.size(); ++i) {
        delete x[i];
    }

    core();
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    printf("%f\t%d\t%d\n", __p__, __k__, __s__);
    printf("%.4fs have been used in all.\n", elapsed_secs);
}
