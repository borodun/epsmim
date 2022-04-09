#ifndef LAB1_WAVESIM_H
#define LAB1_WAVESIM_H

#include "waveparams.h"
#include "exceptions.h"
#include <cmath>
#include <utility>
#include <fstream>
#include <cerrno>
#include <cstring>
#include <chrono>
#include <immintrin.h>

using namespace std::chrono;

class WaveSim {
private:
    WaveParams params;

    double pi = 3.1415926535897;
    double dpi = 2 * 3.1415926535897;
    double f0 = 1.0;
    double t0 = 1.5;
    double lambda = 1 / (4.0 * 4.0);
    double Xa = 0.0, Xb = 4.0, Ya = 0.0, Yb = 4.0;
    double hx;
    double hy;
    double hxsqr;
    double hysqr;
    double tau;
    double tausqr;
    int vectorSize = 4;

    int Sy, Sx;

    double *UNext;
    double *UCurrent;
    double *UPrev;
    double *P;

    __m256i leftShift = _mm256_set_epi64x(2, 1, 0, 0);
    __m256i rightShift = _mm256_set_epi64x(3, 3, 2, 1);
    __m256i swapFirstLast = _mm256_set_epi64x(0, 2, 1, 3);
    __m256d leftBlend = _mm256_set_pd(0, 0, 0, ~0);
    __m256d rightBlend = _mm256_set_pd(~0, 0, 0, 0);

    __m256d v_hxsqr;
    __m256d v_hysqr;
    __m256d v_tausqr;

    int x = 0, y = 0, n = 0;

    void init();

    __m256d phi();

    void saveBinary();

    void saveToFile();

    void printArrays();

    void printArray(double *arr);

public:
    explicit WaveSim(WaveParams *params);

    void Run();
};


#endif
