/* ==========================================
 *  SSE for the 1D FM-TFIM (OBC)
 *  Author: Yi-Ming Ding
 *  Email: dingyiming@westlake.edu.cn
 *  Last updated: April 22, 2024
 * =========================================*/
#include "tfim.hpp"

int main(int argc, char **argv) {
    clock_t t0 = time(nullptr);

    const int L = std::stoi(argv[1]);
    const double beta = std::stod(argv[2]);
    const double J = std::stod(argv[3]);
    const double h = std::stod(argv[4]);
    const int nThm = std::stoi(argv[5]);
    const int nStat = std::stoi(argv[6]);
    const int nBins = std::stoi(argv[7]);

    auto *model = new TFIM(L, beta, J, h);
    model->reportEnv();
    // Thermalization
    for (int i = 0; i < nThm; ++i) {
        model->updateConfig();
        model->adjustM();
    }
    // Measurements
    for (int k = 0; k < nBins; ++k) {
        model->iniMeasure();
        for (int i = 0; i < nStat; ++i) {
            model->updateConfig();
            model->measure();
        }
        model->statisticize();
        model->saveData();
    }

    delete model;
    clock_t t1 = time(nullptr);
    timer(t0, t1);
    return 0;
}


