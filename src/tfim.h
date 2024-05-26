/* ==========================================
 *  SSE for the 1D FM-TFIM (OBC)
 *  Author: Yi-Ming Ding
 *  Email: dingyiming@westlake.edu.cn
 *  Last updated: April 22, 2024
 * =========================================*/
#include <iostream>
#include <cstdlib>
#include <ctime>

static double randProb() {
    return static_cast<double>(rand() % RAND_MAX) / static_cast<double>(RAND_MAX);
}

class TFIM {
public:
    TFIM(int L, double beta, double J, double h);
    ~TFIM();
    // ===================
    //  Basics
    // ===================
    int L, nSites, nBonds;
    double beta;
    int n, M;
    double J, h;

    int *spins;
    int **bSites;
    int *opString;
    int *vertexList;
    int *vFirst;
    int *vLast;

    // ==========================
    //  Config Update
    // ==========================
    void diagUpdate();
    void makeVertexList() const;
    void clusterUpdate();
    void adjustM();
    void updateConfig();

    // ==========================
    //  Measurements
    // ==========================
    double energy;
    double *corr_s0_sj;
    void iniMeasure();      // initialize relevant quantities for a new measurement
    void measure();         // do measurements
    void statisticize();    // do calculations
    void saveData() const;
    void reportEnv() const;

private:
    // ===================
    //  Basics
    // ===================
    double selectionProb;
    double addFactor, removeFactor;

    // ==========================
    //  Internal stack
    // ==========================
    // Pre-allocating memory makes it faster than "std::stack"
    int *stack;
    int top;
    void stack_initialize();
    void stack_push(int x);
    int stack_pop();
    int flip;

    // ==========================
    //  Auxiliary modules
    // ==========================
    int nMeasure;
    [[nodiscard]] double measureCorr(int si, int sj) const;
    void makeCluster();
    [[nodiscard]] int randBond();
    [[nodiscard]] int randSite();
};

TFIM::TFIM(int L, double beta, double J, double h) {
    srand(time(nullptr));

    // ===================
    //  Assign params
    // ===================
    this->L = L;
    this->beta = beta;
    this->J = J;
    this->h = h;
    this->nSites = L;
    this->nBonds = L;
    this->n = 0;
    this->M = 20;

    this->selectionProb = (nSites * h) / (nBonds * 2 * J + nSites * h);
    this->addFactor = beta * (nSites * h + 2 * J * nBonds);
    this->removeFactor = 1.0 / addFactor;

    // ===============================
    //  Initialize data structures
    // ===============================
    spins = new int [nSites];
    for (int s = 0; s < nSites; ++s)
        randProb() < 0.5 ? spins[s] = 1 : spins[s] = -1;

    bSites = new int * [nBonds];
    for (int b = 0; b < nBonds; ++b)
        bSites[b] = new int [2];
    for (int b = 0; b < nBonds; ++b) {
        bSites[b][0] = b;
        bSites[b][1] = (b + 1) % L;
    }

    opString = new int [M];
    for (int i = 0; i < M; ++i)
        opString[i] = 0;

    vertexList = new int [4 * M];
    stack = new int [8 * M];
    vFirst = new int [nSites];
    vLast = new int [nSites];

    // ===============================
    //  Initialize for measurements
    // ===============================
    corr_s0_sj = new double [nSites];
}

TFIM::~TFIM() {
    delete[] spins;
    delete[] bSites;
    delete[] opString;
    delete[] vFirst;
    delete[] vLast;
    delete[] vertexList;
    delete[] stack;
    delete[] corr_s0_sj;
}

void timer(clock_t t0, clock_t t1) {
    int interval = (int) (t1 - t0);
    auto hour = interval / 3600;
    auto min = (interval - 3600 * hour) / 60;
    auto sec = interval - 3600 * hour - 60 * min;
    std::cout << "■ Runtime: " << hour << "h-" << min << "m-" << sec << "s" << std::endl;
}