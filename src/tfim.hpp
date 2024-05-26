/* ==========================================
 *  SSE for the 1D FM-TFIM (OBC)
 *  Author: Yi-Ming Ding
 *  Email: dingyiming@westlake.edu.cn
 *  Last updated: April 22, 2024
 * =========================================*/
#include "tfim.h"
#include <string>
#include <fstream>

void TFIM::diagUpdate() {
    // -------------------------------------------------
    //  Define op = 4 * b + a + 2,
    //  then
    //      [0, 0], op := 0, identity
    //      [-1, s], op = 4 * s + 1, off-diag, site
    //      [0, s], op = 4 * s + 2, diag, site
    //      [1, s], op = 4 * b + 3, diag, bond
    // -------------------------------------------------
    int op, newBond, newSite, theSite;

    for (int p = 0; p < M; ++p) {
        op = opString[p];
        if (op == 0) {
            // First decide whether to insert, then decide what to insert
            if (randProb() * (M - n) < addFactor) {
                if (randProb() < selectionProb) {
                    newSite = randSite();
                    opString[p] = 4 * newSite + 2;
                    n += 1;
                }

                else {
                    newBond = randBond();
                    if (spins[bSites[newBond][0]] == spins[bSites[newBond][1]]) {
                        opString[p] = 4 * newBond + 3;
                        n += 1;
                    }
                }
            }
        }

        else if ((op % 4 == 2) or (op % 4 == 3)) {
            if (removeFactor * (M - n + 1) > randProb()) {
                opString[p] = 0;
                n -= 1;
            }
        }

        else {
            theSite = op / 4;
            spins[theSite] *= -1;
        }
    }
}

void TFIM::makeVertexList() const {
    int op, b_p;
    int v_leg0;
    int s0, s1;
    int s0_vLast, s1_vLast;

    // ===========================================
    //  Initialize vertexList, vFirst, vLast
    // ===========================================
    for (int i = 0; i < 4 * M; ++i)
        vertexList[i] = -1;
    for (int s = 0; s < nSites; ++s) {
        vFirst[s] = -1;
        vLast[s] = -1;
    }

    // ===========================================
    //  Make the list
    //      Similar to HM, but consider bond-type
    //      and site-type respectively.
    // ===========================================
    for (int p = 0; p < M; ++p) {
        op = opString[p];
        if (op != 0) {
            // -.-.-.-.-.-.-.-.-.
            //  Bond operator
            // -.-.-.-.-.-.-.-.-.
            if (op % 4 == 3) {
                b_p = op / 4;
                s0 = bSites[b_p][0];
                s1 = bSites[b_p][1];
                v_leg0 = 4 * p;
                s0_vLast = vLast[s0];
                s1_vLast = vLast[s1];

                if (s0_vLast > -1) {
                    vertexList[s0_vLast] = v_leg0;
                    vertexList[v_leg0] = s0_vLast;
                } else
                    vFirst[s0] = v_leg0;
                vLast[s0] = v_leg0 + 2;

                if (s1_vLast > -1) {
                    vertexList[s1_vLast] = v_leg0 + 1;
                    vertexList[v_leg0 + 1] = s1_vLast;
                } else
                    vFirst[s1] = v_leg0 + 1;
                vLast[s1] = v_leg0 + 3;
            }

            // -.-.-.-.-.-.-.-.-.
            //  Site operator
            // -.-.-.-.-.-.-.-.-.
            else if ((op % 4 == 2) or (op % 4 == 1)) {
                s0 = op / 4;
                v_leg0 = 4 * p;
                s0_vLast = vLast[s0];

                if (s0_vLast > -1) {
                    vertexList[s0_vLast] = v_leg0;
                    vertexList[v_leg0] = s0_vLast;
                } else
                    vFirst[s0] = v_leg0;
                vLast[s0] = v_leg0 + 2;
            }
        }
    }

    // PBC correction
    int s_vFirst;
    int s_vLast;
    for (int s = 0; s < nSites; ++s) {
        s_vFirst = vFirst[s];
        s_vLast = vLast[s];

        if (s_vFirst != -1) {
            vertexList[s_vFirst] = s_vLast;
            vertexList[s_vLast] = s_vFirst;
        }
    }
}


void TFIM::clusterUpdate() {
    // =============================================
    //  Form and flip cluster
    // =============================================
    stack_initialize();

    for (int v = 0; v < 4 * M; v += 2) {
        if (vertexList[v] < 0)  // with no vertex or already in some cluster
            continue;

        // Decide whether to flip site operators
        // "-2" means doing flip and "-1" for otherwise
        (randProb() > 0.5) ? flip = -1 : flip = -2;
        stack_push(v);
        while (top != -1) {
            makeCluster();
        }
    }

    // =================================================
    //  Update spins
    // =================================================
    for (int i = 0; i < nSites; ++i) {
        if (vFirst[i] != -1) {
            if (vertexList[vFirst[i]] == -2)
                spins[i] *= -1;
        }

        else if (randProb() < 0.5)
            spins[i] *= -1;
    }
}

void TFIM::makeCluster() {
    int p, op;
    int vs = stack_pop();   // v start
    int v1 = vertexList[vs];

    if (v1 < 0)
        // This means "vs" links to no vertex now,
        //      which CAN ONLY HAPPEN if "vs" HAS BEEN processed (see codes below).
        // so then no need to operate "vs" or link it to other vertices, just return.
        return;
    else if (vertexList[v1] != flip)
        // The state "flip" also denotes that the vertex has been processed,
        //      so "v1" has not been processed in this case, and we can push "v1".
        // It is possible that "vs" is not processed yet but "v1" has been processed:
        //      when the vertex goes to some end point and checks back.
        // Therefore, an extra judgement is needed.
        stack_push(v1);

    // ======================================
    //  Process "vs" (I) flip the operator
    // ======================================
    if (flip == -2) {
        p = vs / 4;
        op = opString[p];
        if (op % 4 == 1)
            opString[p] += 1;
        else if (op % 4 == 2)
            opString[p] -= 1;
    }

    // ==============================================================================
    //  Process "vs" (II) link to other vertices if "vs" belong to a bond operator
    // ==============================================================================
    if (opString[vs / 4] % 4 == 3) {
        int v2 = vs ^ 1;
        int v3 = vs ^ 2;
        int v4 = v3 ^ 1;
        if (vertexList[v2] >= 0)
            stack_push(v2);
        if (vertexList[v3] >= 0)
            stack_push(v3);
        if (vertexList[v4] >= 0)
            stack_push(v4);
    }

    // ==============================================================================
    //  Mark that we have processed "vs"
    // ==============================================================================
    vertexList[vs] = flip;
}

void TFIM::adjustM() {
    int newM = n + n / 3;
    if (M < newM){
        int *opString_copy = new int [M];
        for (int i = 0; i < M; ++i)
            opString_copy[i] = opString[i];
        delete[] opString;
        opString = new int [newM];
        for (int i = 0; i < M; ++i)
            opString[i] = opString_copy[i];
        for (int i = M; i < newM; ++i)
            opString[i] = 0;
        M = newM;
        delete[] vertexList;
        vertexList = new int[4 * M];
        delete[] stack;
        stack = new int [8 * M];
    }
}

void TFIM::updateConfig() {
    diagUpdate();
    makeVertexList();
    clusterUpdate();
}

// ==========================
//  Internal stack
// ==========================
void TFIM::stack_initialize() {
    top = -1;
}

void TFIM::stack_push(int x) {
    top++;
    stack[top] = x;
}

int TFIM::stack_pop() {
    int top_val = stack[top];
    top--;
    return top_val;
}

// ==========================
//  Measurements
// ==========================
double TFIM::measureCorr(int si, int sj) const {
    return static_cast<double>(spins[si] * spins[sj]);
}

void TFIM::iniMeasure() {
    for (int i = 0; i < nSites; ++i)
        corr_s0_sj[i] = 0.0;

    energy = 0.0;
    nMeasure = 0;
}

void TFIM::measure() {
    for (int i = 0; i < nSites; ++i)
        corr_s0_sj[i] += measureCorr(0, i);
    energy += -static_cast<double>(n) / beta;   // + nBonds * J + nSites * h;
    nMeasure += 1;
}


void TFIM::statisticize() {
    for (int i = 0; i < nSites; ++i)
        corr_s0_sj[i] /= static_cast<double>(nMeasure);

    energy /= static_cast<double>(nMeasure);
}


void TFIM::saveData() const {
    std::ofstream f0;
    std::string path0 = "../data/energy.dat";
    std::string path1 = "../data/corr_s0_sj.dat";
    f0.precision(16);

    f0.open(path0, std::ios::app);
    f0 << energy << std::endl;
    f0.close();

    f0.open(path1, std::ios::app);
    for (int i = 0; i < nSites; ++i)
        f0 << corr_s0_sj[i] << std::endl;
    f0.close();
}


// ==========================
//  Auxiliary modules
// ==========================

int TFIM::randBond() {
    return (rand() % nBonds);
}

int TFIM::randSite() {
    return (rand() % nSites);
}

void TFIM::reportEnv() const {
    std::cout << "■ SSE for 1D FM-TFIM (PBC)\n\tH = -J * sumZZ - h * sumX:" << std::endl;
    std::cout << "\tL = " << L << ", beta = " << beta << ", J = " << J << ", h = " << h << std::endl;
}