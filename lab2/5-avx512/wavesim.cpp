#include "wavesim.h"

WaveSim::WaveSim(WaveParams *params) : params(*params) {
    int Nx = params->getNx();
    int Ny = params->getNy();

    Sx = params->getSx();
    Sy = params->getSy();

    hx = (Xb - Xa) / Nx;
    hy = (Yb - Ya) / Ny;
    hxsqr = 1 / (2 * hx * hx);
    hysqr = 1 / (2 * hy * hy);
    tau = Nx <= 1000 && Ny <= 1000 ? 0.01 : 0.001;
    tausqr = tau * tau;

    UCurrent = static_cast<double *>(aligned_alloc(64, Ny * Nx * sizeof(double)));
    UPrev = static_cast<double *>(aligned_alloc(64, Ny * Nx * sizeof(double)));
    P = static_cast<double *>(aligned_alloc(64, Ny * Nx * sizeof(double)));

    UNext = UPrev;

    v_hxsqr = _mm512_set1_pd(hxsqr);
    v_hysqr = _mm512_set1_pd(hysqr);
    v_tausqr = _mm512_set1_pd(tausqr);
}

void WaveSim::init() {
    int Nx = params.getNx();
    int Ny = params.getNy();
    for (y = 0; y < Ny; ++y) {
        for (x = 0; x < Nx; ++x) {
            if (x < Nx / 2) {
                P[y * Nx + x] = 0.1 * 0.1;
            } else {
                P[y * Nx + x] = 0.2 * 0.2;
            }
            UCurrent[y * Nx + x] = 0;
            UPrev[y * Nx + x] = 0;
        }
    }
}

void WaveSim::Run() {
    init();
    //params.PrintInfo();

    int Nt = params.getNt();
    int Nx = params.getNx();
    int Ny = params.getNy();

    int leftBlend = 0b00000001;
    int rightBlend = 0b10000000;

    double UMax = 0;
    __m512d VUMax = _mm512_set1_pd(0);

    time_point<high_resolution_clock> simStart = high_resolution_clock::now();
    for (n = 0; n < Nt; ++n) {

        auto *UUp = (__m512d *) (UCurrent);
        auto *UMid = (__m512d *) (UCurrent + Nx);
        auto *ULow = (__m512d *) (UCurrent + 2 * Nx);

        auto *UPrevVec = (__m512d *) (UPrev + Nx);
        auto *UNextVec = (__m512d *) (UNext + Nx);

        auto *PMid = (__m512d *) (P + Nx);
        auto *PLow = (__m512d *) (P + 2 * Nx);

        //time_point<high_resolution_clock> stepStart = high_resolution_clock::now();
        for (y = 1; y < Ny - 1; ++y) {
            __m512d UUpLeft = UUp[0];
            __m512d UUpCentre = UUp[1];
            __m512d UUpRight = UUp[2];

            __m512d UMidLeft = UMid[0];
            __m512d UMidCentre = UMid[1];
            __m512d UMidRight = UMid[2];

            __m512d ULowLeft = ULow[0];
            __m512d ULowCentre = ULow[1];
            __m512d ULowRight = ULow[2];

            __m512d PMidLeft = PMid[0];
            __m512d PMidCentre = PMid[1];
            __m512d PMidRight = PMid[2];

            __m512d PLowLeft = PLow[0];
            __m512d PLowCentre = PLow[1];
            __m512d PLowRight = PLow[2];

            int row = y * Nx;
            int prevRow = (y - 1) * Nx;
            int nextRow = (y + 1) * Nx;

            for (x = 1; x < vectorSize; ++x) {
                int pos = row + x;
                int prevPos = prevRow + x;

                double currentPos = UCurrent[pos];
                double pPos = P[pos];
                double pm1 = P[pos - 1];
                double ppm1 = P[prevPos - 1];

                double avgx =
                        ((UCurrent[pos + 1] - currentPos) *
                         (P[prevPos] + pPos) +
                         (UCurrent[pos - 1] - currentPos) *
                         (ppm1 + pm1)) * hxsqr;
                double avgy =
                        ((UCurrent[nextRow + x] - currentPos) *
                         (pm1 + pPos) +
                         (UCurrent[prevPos] - currentPos) *
                         (ppm1 + P[prevPos])) * hysqr;
                double result = 2 * currentPos - UPrev[pos] + tausqr * (phi() + avgx + avgy);
                UNext[pos] = result;
                if (result > UMax) {
                    UMax = std::abs(result);
                }
            }

            for (x = 1; x < Nx / vectorSize; ++x) {
                __m512d U1 = _mm512_mask_blend_pd(leftBlend, _mm512_permutexvar_pd(leftShift, UMidCentre),
                                                  _mm512_permutexvar_pd(swapFirstLast, UMidRight)) - UMidCentre;
                __m512d P1 = PLowCentre + PMidCentre;
                __m512d U2 = _mm512_mask_blend_pd(rightBlend, _mm512_permutexvar_pd(rightShift, UMidCentre),
                                                  _mm512_permutexvar_pd(swapFirstLast, UMidLeft)) - UMidCentre;
                __m512d PLeftShift = _mm512_mask_blend_pd(rightBlend, _mm512_permutexvar_pd(rightShift, PMidCentre),
                                                          _mm512_permutexvar_pd(swapFirstLast, PMidLeft));
                __m512d PDiagShift = _mm512_mask_blend_pd(rightBlend, _mm512_permutexvar_pd(rightShift, PLowCentre),
                                                          _mm512_permutexvar_pd(swapFirstLast, PLowLeft));
                __m512d P2 = PDiagShift + PLeftShift;
                __m512d avgx = (U1 * P1 + U2 * P2) * v_hxsqr;

                __m512d U3 = UUpCentre - UMidCentre;
                __m512d P3 = PLeftShift + PMidCentre;
                __m512d U4 = ULowCentre - UMidCentre;
                __m512d P4 = PDiagShift + PLowCentre;
                __m512d avgy = (U3 * P3 + U4 * P4) * v_hysqr;

                __m512d res = 2 * UMidCentre - UPrevVec[x] + v_tausqr * (avgx + avgy + phi());

                UNextVec[x] = res;

                VUMax = _mm512_max_pd(VUMax, res);

                UUpLeft = UUpCentre;
                UUpCentre = UUpRight;
                UUpRight = UUp[x + 2];

                UMidLeft = UMidCentre;
                UMidCentre = UMidRight;
                UMidRight = UMid[x + 2];

                ULowLeft = ULowCentre;
                ULowCentre = ULowRight;
                ULowRight = ULow[x + 2];

                PMidLeft = PMidCentre;
                PMidCentre = PMidRight;
                PMidRight = PMid[x + 2];

                PLowLeft = PLowCentre;
                PLowCentre = PLowRight;
                PLowRight = PLow[x + 2];
            }

            for (x = 1; x < vectorSize; ++x) {
                int pos = row + x;
                int prevPos = prevRow + x;

                double currentPos = UCurrent[pos];
                double pPos = P[pos];
                double pm1 = P[pos - 1];
                double ppm1 = P[prevPos - 1];

                double avgx =
                        ((UCurrent[pos + 1] - currentPos) *
                         (P[prevPos] + pPos) +
                         (UCurrent[pos - 1] - currentPos) *
                         (ppm1 + pm1)) * hxsqr;
                double avgy =
                        ((UCurrent[nextRow + x] - currentPos) *
                         (pm1 + pPos) +
                         (UCurrent[prevPos] - currentPos) *
                         (ppm1 + P[prevPos])) * hysqr;
                double result = 2 * currentPos - UPrev[pos] + tausqr * (phi() + avgx + avgy);
                UNext[pos] = result;
                if (result > UMax) {
                    UMax = std::abs(result);
                }
            }

            ULow = UMid;
            UMid = UUp;
            UUp += Nx / vectorSize;

            UPrevVec += Nx / vectorSize;

            PLow = PMid;
            PMid += Nx / vectorSize;
        }
        /*time_point<high_resolution_clock> stepEnd = high_resolution_clock::now();
        nanoseconds stepTime = duration_cast<nanoseconds>(stepEnd - stepStart);
        nanoseconds totalTime = duration_cast<nanoseconds>(stepEnd - simStart);
        std::cout << n << ": " << stepTime.count() / 1e9 << "s (" << totalTime.count() / 1e9 << "s)" << std::endl;*/

        std::swap(UCurrent, UPrev);
        UNext = UPrev;
    }
    time_point<high_resolution_clock> simEnd = high_resolution_clock::now();
    nanoseconds simTime = duration_cast<nanoseconds>(simEnd - simStart);
    std::cout << std::endl << "Time: " << simTime.count() / 1e9 << "s" << std::endl;

    double umax = 0.0;
    auto *max = (double *) (&VUMax);
    for (int i = 0; i < vectorSize; i++) {
        umax = std::max(umax, max[i]);
    }
    umax = std::max(umax, UMax);
    std::cout << " Max: " << umax << std::endl;

    saveBinary();
    //saveToFile();
    free(UCurrent);
    free(UPrev);
    free(P);
}

__m512d WaveSim::vphi() {
    __m512d result = _mm512_set1_pd(0);
    if (x * vectorSize <= Sx && Sx < (x + 1) * vectorSize && y == Sy) {
        double part = dpi * (n * tau - t0);
        double ex = exp(-part * part * lambda);
        double sine = sin(part);
        double res = ex * sine;

        int pos = Sx - x * vectorSize;
        result[pos] = res;
    }
    return result;
}

double WaveSim::phi() {
    double result = 0;
    if (x == Sx && y == Sy) {
        double part = dpi * (n * tau - t0);
        double ex = exp(-part * part * lambda);
        double sine = sin(part);
        result = ex * sine;
    }
    return result;
}

void WaveSim::saveBinary() {
    int Nx = params.getNx();
    int Ny = params.getNy();
    std::fstream of(params.getFilename() + ".bin", std::ios::out | std::ios::binary);
    if (!of.is_open()) {
        throw fileException(params.getFilename() + ".bin cannot be open");
    }

    for (y = 0; y < Ny; ++y) {
        for (x = 0; x < Nx; ++x) {
            of.write((char *) &UCurrent[y * Nx + x], sizeof(double));
        }
    }

    of.close();
}

void WaveSim::saveToFile() {
    int Nx = params.getNx();
    int Ny = params.getNy();
    std::fstream of(params.getFilename() + ".dat", std::ios::out);
    if (!of.is_open()) {
        throw fileException(params.getFilename() + ".dat cannot be open");
    }
    for (y = 0; y < Ny; ++y) {
        for (x = 0; x < Nx; ++x) {
            of << UCurrent[y * Nx + x] << " ";
        }
        of << std::endl;
    }

    of.close();
}

void WaveSim::printArray(double *arr) {
    int Nx = params.getNx();
    int Ny = params.getNy();
    for (y = 0; y < Ny; ++y) {
        for (x = 0; x < Nx; ++x) {
            std::cout << arr[y * Nx + x] << " ";
        }
        std::cout << std::endl;
    }
}

void WaveSim::printArrays() {
    std::cout << "Next" << std::endl;
    printArray(UNext);
    std::cout << "Current" << std::endl;
    printArray(UCurrent);
    std::cout << "Prev" << std::endl;
    printArray(UPrev);
}

