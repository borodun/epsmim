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
#include "omp.h"

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
    int vectorSize = 8;

    int Sy, Sx;

    double *UNext;
    double *UCurrent;
    double *UPrev;
    double *P;

    __m512i leftShift = _mm512_set_epi64(6, 5, 4, 3, 2, 1, 0, 0);
    __m512i rightShift = _mm512_set_epi64(7, 7, 6, 5, 4, 3, 2, 1);
    __m512i swapFirstLast = _mm512_set_epi64(0, 6, 5, 4, 3, 2, 1, 7);

    __m512d v_hxsqr;
    __m512d v_hysqr;
    __m512d v_tausqr;

//    int x = 0, y = 0, n = 0;

    void init();

    __m512d vphi(int x, int y, int n);

    double phi(int x, int y, int n);

    void saveBinary();

    void saveToFile();

    void printArrays();

    void printArray(double *arr);

public:
    explicit WaveSim(WaveParams *params);

    void Run();
};


#endif
