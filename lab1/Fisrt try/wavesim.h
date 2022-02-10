#ifndef LAB1_WAVESIM_H
#define LAB1_WAVESIM_H

#include "waveparams.h"
#include "../exceptions.h"
#include <cmath>
#include <utility>
#include <fstream>
#include <cerrno>
#include <cstring>
#include <chrono>

using namespace std::chrono;

class WaveSim {
private:
    WaveParams params;

    double pi = 3.1415926535897;
    double f0 = 1.0;
    double t0 = 1.5;
    double lambda = 4.0;
    double Xa = 0.0, Xb = 4.0, Ya = 0.0, Yb = 4.0;
    double hx;
    double hy;
    double tau;

    double **UNext;
    double **UCurrent;
    double **UPrev;
    double **P;

    int x, y, n;

    void init();

    double phi();

    void saveBinary();

    void saveToFile();

    void printArray(double **arr);

public:
    explicit WaveSim(WaveParams *params);

    void Run();
};


#endif
