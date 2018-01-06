#include <stdio.h>
#include <vector>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <tuple>

using namespace std;

class DataNode {
private:
    double measure;
    vector<DataNode*> children;
    tuple<double, double> center;
    int index;
public:
    DataNode(double _measure, tuple<double, double> _center, int _index) {
        measure = _measure;
        center = _center;
        index = _index;
    }
    DataNode(double _measure, vector<DataNode*> _children, int _index) {
        measure = _measure;
        children = _children;
        double xsum = 0.0, ysum = 0.0;
        for (unsigned int i = 0; i < children.size(); ++i) {
            xsum += get<0>(_children[i]->center);
            ysum += get<1>(_children[i]->center);
        }
        center = make_tuple(xsum/children.size(), ysum/children.size());
        index = _index;
    }

    double getMeasure() const {
        return measure;
    }

    const vector<DataNode*>& getChildren() const {
        return children;
    }

    void print() {
        printf("DataNode %d, ", this);
        if (children.size() == 0) {
            printf("leaf with measure %f.\n", measure);
        } else {
            printf("internode with measure %f.\n\tChildren: ", measure);
            for (unsigned int i = 0; i < children.size(); ++i) {
                printf("%d ,", children[i]);
            }
            printf("\n");
        }
    }

    tuple<double, double> getCenter() const {
        return center;
    }

    int getIndex() const {
        return index;
    }
};

class Decomposition {
private:
    int scale;
    vector<DataNode*> partition;
public:
    Decomposition(int _scale, vector<DataNode*> _partition) {
        scale = _scale;
        partition = _partition;
    }

    void assignScale(int _scale) {
        scale = _scale;
    }

    void print() {
        printf("Level %d decomposition, with %d partition(s).\n", scale, partition.size());
        for (unsigned int i = 0; i < partition.size(); ++i) {
            partition[i]->print();
        }
    }

    const vector<DataNode*>& getDataNodes() const {
        return partition;
    }
};

class DecompositionChain {
private:
    vector<Decomposition> chain;
public:
    DecompositionChain(const vector<vector<double>> &X, int res, int clusterSize) {
        // uses the 2-d square structure
        // start from scale J, the finest one; however we don't know J in advance, so we will assign that later.
        vector<DataNode*> dataNodes;
        dataNodes.reserve(res*res);
        for (int i = 0; i < res; ++i) {
            for (int j = 0; j < res; ++j) {
                // we visit the image X in column-first order, that (0, 0) -> (1, 0) -> ...
                dataNodes.push_back(new DataNode(X[j][i], make_tuple(j, i), j+i*res));
            }
        }
        chain.push_back(Decomposition(0, dataNodes));
        int scaledRes = res;

        while (scaledRes != 1) {
            vector<int> partitionInIndex;
            if (scaledRes < 2*clusterSize) {
                partitionInIndex = vector<int>(1, scaledRes);
            } else {
                partitionInIndex = vector<int>(std::max(0, scaledRes/clusterSize-1), clusterSize);
                partitionInIndex.push_back(scaledRes%clusterSize+clusterSize);
            }
            vector<DataNode*> parentDataNodes;
            parentDataNodes.reserve(partitionInIndex.size()*partitionInIndex.size());

            for (unsigned int i = 0; i < partitionInIndex.size(); ++i) {
                for (unsigned int j = 0; j < partitionInIndex.size(); ++j) {
                    double mearsure = 0.0;
                    vector<DataNode*> children;
                    for (int k = 0; k < partitionInIndex[i]; ++k) {
                        for (int l = 0; l < partitionInIndex[j]; ++l) {
                            int index = j*clusterSize+l + (i*clusterSize+k)*scaledRes;
                            mearsure += dataNodes[index]->getMeasure();
                            children.push_back(dataNodes[index]);
                        }
                    }
                    parentDataNodes.push_back(new DataNode(mearsure, children, j+i*partitionInIndex.size()));
                }
            }

            dataNodes = parentDataNodes;
            scaledRes = partitionInIndex.size();
            chain.push_back(Decomposition(0, dataNodes));
        }

        // never forget to reverse the chain (because we build the finest one first)
        std::reverse(chain.begin(), chain.end());
        for (unsigned int i = 0; i < chain.size(); ++i) {
            chain[i].assignScale(i);
        }
    }

    void print() {
        printf("%d decomposition in chains:\n", chain.size());
        for (unsigned int i = 0; i < chain.size(); ++i) {
            chain[i].print();
        }
        printf("========================\n");
    }

    const vector<Decomposition>& getDecompositions() const {
        return chain;
    }
};

class Transport {
public:
    DataNode *x, *y;
    double amount;

    Transport(DataNode *_x, DataNode *_y, double _amount):x(_x), y(_y), amount(_amount) {}
};

class TransportPlan {
private:
    DecompositionChain *X, *Y;
    vector<vector<vector<double>>> *costChain;
    vector<vector<double>> cost;
    int level;
    vector<Transport> category;

    TransportPlan(vector<vector<vector<double>>> *_costChain, DecompositionChain *_X, DecompositionChain *_Y, int _level) {
        X = _X;
        Y = _Y;
        costChain = _costChain;
        level = _level + 1;
    }
public:
    TransportPlan(vector<vector<vector<double>>> *_costChain, DecompositionChain *_X, DecompositionChain *_Y) {
        X = _X;
        Y = _Y;
        costChain = _costChain;
        level = 0;
        category.push_back(Transport(X->getDecompositions()[0].getDataNodes()[0], Y->getDecompositions()[0].getDataNodes()[0], 1.0));
    }

    TransportPlan propagate() {
        TransportPlan newPlan = TransportPlan(costChain, X, Y, level + 1);
        for (auto iterator = category.begin(); iterator != category.end(); ++iterator) {

        }
        return newPlan;
    }

    void refine() {
        // blabla
    }
};

vector<vector<vector<double>>> decomposeCost(const vector<vector<double>> &cost, const DecompositionChain &Xchain, const DecompositionChain &Ychain) {
    int depth = Xchain.getDecompositions().size();
    vector<vector<vector<double>>> costChain(depth);
    costChain[depth-1] = cost;

    // very naive estimate of the distance cost
    for (int l = depth-2; l >= 0; --l) {
        auto xDecompNodes = Xchain.getDecompositions()[l].getDataNodes(),
            yDecompNodes = Ychain.getDecompositions()[l].getDataNodes();
        costChain[l].resize(xDecompNodes.size());
        for (auto i = 0; i < xDecompNodes.size(); ++i) {
            auto xNodeChildren = xDecompNodes[i]->getChildren();
            costChain[l][i].resize(yDecompNodes.size());
            for (auto j = 0; j < yDecompNodes.size(); ++j) {
                auto yNodeChildren = yDecompNodes[j]->getChildren();
                double costSum = 0;
                for (auto xiter = xNodeChildren.begin(); xiter != xNodeChildren.end(); ++xiter) {
                    for (auto yiter = yNodeChildren.begin(); yiter != yNodeChildren.end(); ++yiter) {
                        costSum += costChain[l+1][(*xiter)->getIndex()][(*yiter)->getIndex()];
                    }
                }
                costChain[l][i][j] = costSum / xNodeChildren.size() / yNodeChildren.size();
            }
        }
    }

    return costChain;
}

void core(const vector<vector<double>> &X, const vector<vector<double>> &Y, const vector<vector<double>> &cost, int res) {
    const int clusterSize = 2;
    DecompositionChain Xchain = DecompositionChain(X, res, clusterSize);
    DecompositionChain Ychain = DecompositionChain(Y, res, clusterSize);

    vector<vector<vector<double>>> costChain = decomposeCost(cost, Xchain, Ychain);

    /*for (auto i = 0; i < costChain.size(); ++i) {
        printf("Level %d cost matrix.\n", i);
        for (auto k = 0; k < costChain[i].size(); ++k) {
            for (auto j = 0; j < costChain[i][k].size(); ++j) {
                printf("% 5.1f ", costChain[i][k][j]);
            }
            printf("\n");
        }
    }*/

    TransportPlan plan(&costChain, &Xchain, &Ychain);

    for (int i = 1; i < Xchain.getDecompositions().size(); ++i) {
        plan = plan.propagate();
        plan.refine();
    }
}

int main() {
    clock_t begin = clock();
    int res = 5;
    int m = res*res, n = res*res;
    vector<vector<double>> cost = vector<vector<double>>(m);
    for (int i = 0; i < m; ++i) {
        int ii = i%res, ij = i/res;
        cost[i].resize(n);
        for (int j = 0; j < n; ++j) {
            int ji = j%res, jj = j/res;
            cost[i][j] = (ii-ji)*(ii-ji) + (ij-jj)*(ij-jj);
        }
    }
    vector<vector<double>> X = vector<vector<double>>(res);
    double Xsum = 0;
    for (int i = 0; i < res; ++i) {
        X[i].resize(res);
        for (int j = 0; j < res; ++j) {
            X[i][j] = (rand()+0.0) / RAND_MAX;
            Xsum += X[i][j];
        }
    }
    for (int i = 0; i < res; ++i) {
        for (int j = 0; j < res; ++j) {
            X[i][j] /= Xsum;
        }
    }
    vector<vector<double>> Y = vector<vector<double>>(res);
    double Ysum = 0;
    for (int i = 0; i < res; ++i) {
        Y[i].resize(res);
        for (int j = 0; j < res; ++j) {
            Y[i][j] = (rand()+0.0) / RAND_MAX;
            Ysum += Y[i][j];
        }
    }
    for (int i = 0; i < res; ++i) {
        for (int j = 0; j < res; ++j) {
            Y[i][j] /= Ysum;
        }
    }

    core(X, Y, cost, res);

    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    printf("%.4fs have been used in all.\n", elapsed_secs);
    return 0;
}