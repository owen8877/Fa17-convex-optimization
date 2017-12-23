#include <vector>
#include <limits>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <ctime>

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
        for (int i = 0; i < n; ++i) {
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

bool mostnegative(int& nr, int& nc) {
    // printf("Relative cost.\n");
    double val = 0;
    for (int i = 0; i < m ; ++i) {
        for (int j = 0; j < n; ++j) {
            double tmp = cost[i][j] - u[i] - v[j];
            // printf("%.1e\t", tmp);
            if (tmp < val) {
                val = tmp;
                nr = i;
                nc = j;
            }
        }
        // printf("\n");
    }
    return (val < -1e-15);
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
    initialSolution();
    // printx();
    int itr = 1;
    while (true) {
        if (itr % 500 == 0) {
            double sum = 0.0;
            for (unsigned int i = 0; i < x.size(); ++i) {
                sum += x[i]->val * cost[x[i]->p][x[i]->q];
            }
            printf("Iteration %d\t%e\n", itr, sum);
        }
        // update the multipliers
        updateuv();
        // printArray(u);
        // printArray(v);
        // returns false if all relative cost are non-negative
        if (!mostnegative(nr, nc)) {
            break;
        }
        corValPair* newNode = new corValPair(nr, nc, 0.0);
        // printf("New location: %d %d\n", nr, nc);
        addNode(newNode);
        // printx();
        // printf("=======\n");
        corValPair* dontForgetToDelete = new corValPair(nr, nc, numeric_limits<double>::infinity());
        if (!happySearch(0, false, newNode, numeric_limits<double>::infinity(), dontForgetToDelete)) {
            printf("CTransimplex:core\tLoop detection failed!\n");
            printx();
            return;
        }
        // printf("[%d, %d] => [%d, %d]\n", globalCor->p, globalCor->q, nr, nc);
        delete dontForgetToDelete;
        rmNode(globalCor);
        // printx();
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
    printf("%.1fs have been used.\n", elapsed_secs);
}