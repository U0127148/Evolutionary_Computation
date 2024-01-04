#include <iostream>
#include <cmath>
#include <bitset>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <time.h>
#include <unordered_map>
using namespace std;
const int POP_SIZE = 500;
const int N = 100;
const int K = 450;
const double P_CROSSOVER = 1;
const double P_MUTATION = (double)1 / N;
const int TERMINATION = 500;
typedef vector<bitset<10> > genoType;
typedef pair<genoType, genoType> parents;
vector<genoType> pool;
unordered_map<int, float> recordMap;

inline double fitness(const genoType& p) {
    double sum = 0;
    for (int i = 0; i < N; i++) {
        int xi = int(p[i].to_ulong()) - 512;
        if (!recordMap.count(xi)) {
            recordMap[xi] = xi * sin(sqrt(abs(xi)));
        }
        sum += recordMap[xi];
    }
    return 418.98291 * N - sum;
}

bool cmp(const genoType& p1, const genoType& p2) { 
    return fitness(p1) < fitness(p2);
}

void initPopulation() {
    for (int i = 0; i < POP_SIZE; i++) { // population size
        genoType p;
        for (int j = 0; j < N; j++) { // choose each xi
            int xj = rand() % 1024; // generate 0 ~ 1023
            bitset<10> xj_b(xj);
            p.emplace_back(xj_b);
        }
        pool.emplace_back(p);
    }
}

void parentSelection(parents& parent) {
    sort(pool.begin(), pool.end(), cmp);
    int idx1, idx2;
    idx1 = idx2 = INT_MAX;
    for (int i = 0, tmp1, tmp2; i < K; i++) {
        tmp1 = rand() % pool.size();
        tmp2 = rand() % pool.size();
        idx1 = min(idx1, tmp1);
        idx2 = min(idx2, tmp2);
    }
    parent.first = pool[idx1];
    parent.second = pool[(idx1 != idx2) ? idx2 : idx2 + 1];
}

void crossover(parents& parent, vector<genoType>& offspring) {
    double p = (double) rand() / (RAND_MAX + 1.0);
    if (p <= P_CROSSOVER) { // 2-point crossover
        int pos1 = 0, pos2 = 0;
        pos1 = rand() % (N - 1);
        pos2 = rand() % (N - 1);
        while(pos1 == pos2) {
            pos2 = rand() % (N - 1);
        }
        if (pos2 < pos1) {
            swap(pos1, pos2);
        }
        genoType child1, child2;
        for (int i = 0; i <= pos1; i++) {
            child1.emplace_back(parent.first[i]);
            child2.emplace_back(parent.second[i]);
        }
        for (int i = pos1 + 1; i <= pos2; i++) {
            child1.emplace_back(parent.second[i]);
            child2.emplace_back(parent.first[i]);
        }
        for (int i = pos2 + 1; i <= (N - 1); i++) {
            child1.emplace_back(parent.first[i]);
            child2.emplace_back(parent.second[i]);
        }
        offspring.emplace_back(child1);
        offspring.emplace_back(child2);
    } else { // copy parents
        offspring.emplace_back(parent.first);
        offspring.emplace_back(parent.second);
    }
}

void mutation(vector<genoType>& offspring) {
    double p;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < 10; j++) {
            p = (double) rand() / (RAND_MAX + 1.0);
            if (p <= P_MUTATION) {
                offspring[offspring.size() - 1][i][j].flip();
            }
            p = (double) rand() / (RAND_MAX + 1.0);
            if (p <= P_MUTATION) {
                offspring[offspring.size() - 2][i][j].flip();
            }
        }
    }
}

void survivorSelection(vector<genoType>& newGenes) {
    pool.insert(pool.end(), newGenes.begin(), newGenes.end());
    sort(pool.begin(), pool.end(), cmp);
    pool = vector<genoType>(pool.begin(), pool.begin() + 100);
}

int main() {
    srand(0);
    initPopulation();
    cout << fitness(pool[0]) << endl;
    parents parent;
    vector<genoType> offspring;
    for (int t = 0; t < 5; t++) {
        for (int i = 0; i < TERMINATION; i++) {
            for (int j = 0; j < POP_SIZE / 2; j++) {
                parentSelection(parent);
                crossover(parent, offspring);
                mutation(offspring);
            }
            survivorSelection(offspring);
            offspring.clear();
            // cout << "[round " << i+1 << "]: " <<fitness(pool[0]) << endl;
        }
        cout << "[round " << t+1 << "]: " <<fitness(pool[0]) << endl;
    }
}
