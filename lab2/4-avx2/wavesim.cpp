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

    UCurrent = static_cast<double *>(aligned_alloc(32, Ny * Nx * sizeof(double)));
    UPrev = static_cast<double *>(aligned_alloc(32, Ny * Nx * sizeof(double)));
    P = static_cast<double *>(aligned_alloc(32, Ny * Nx * sizeof(double)));

    UNext = UPrev;

    v_hxsqr = _mm256_set1_pd(hxsqr);
    v_hysqr = _mm256_set1_pd(hysqr);
    v_tausqr = _mm256_set1_pd(tausqr);
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

void printVec(std::string name, __m256d vec) {
    std::cout << name << ": ";
    for (int i = 0; i < 4; i++) {
        std::cout << vec[i] << " ";
    }
    std::cout << std::endl;
}

void WaveSim::Run() {
    init();
    //params.PrintInfo();

    int Nt = params.getNt();
    int Nx = params.getNx();
    int Ny = params.getNy();
    __m256d UMax = _mm256_set1_pd(0);

    time_point<high_resolution_clock> simStart = high_resolution_clock::now();
    for (n = 0; n < Nt; ++n) {
        //std::cout << "Layer " << n << std::endl;

        auto *UUp = (__m256d *) (UCurrent);
        auto *UMid = (__m256d *) (UCurrent + Nx);
        auto *ULow = (__m256d *) (UCurrent + 2 * Nx);

        auto *UPrevVec = (__m256d *) (UPrev + Nx);
        auto *UNextVec = (__m256d *) (UNext + Nx);

        auto *PMid = (__m256d *) (P + Nx);
        auto *PLow = (__m256d *) (P + 2 * Nx);

        //time_point<high_resolution_clock> stepStart = high_resolution_clock::now();
        for (y = 1; y < Ny - 1; ++y) {
//            std::cout << "Row " << y << std::endl;
            __m256d UUpLeft = UUp[0];
            __m256d UUpCentre = UUp[1];
            __m256d UUpRight = UUp[2];

            __m256d UMidLeft = UMid[0];
            __m256d UMidCentre = UMid[1];
            __m256d UMidRight = UMid[2];

            __m256d ULowLeft = ULow[0];
            __m256d ULowCentre = ULow[1];
            __m256d ULowRight = ULow[2];

            __m256d PMidLeft = PMid[0];
            __m256d PMidCentre = PMid[1];
            __m256d PMidRight = PMid[2];

            __m256d PLowLeft = PLow[0];
            __m256d PLowCentre = PLow[1];
            __m256d PLowRight = PLow[2];

            for (x = 1; x < Nx / vectorSize; ++x) {
//                std::cout << "Pos " << x << std::endl;
                __m256d U1 = _mm256_blendv_pd(_mm256_permutevar_pd(UMidCentre, leftShift),
                                              _mm256_permutevar_pd(UMidRight, swapFirstLast),
                                              leftBlend) - UMidCentre;
                __m256d P1 = PLowCentre + PMidCentre;
                __m256d U2 = _mm256_blendv_pd(_mm256_permutevar_pd(UMidCentre, rightShift),
                                              _mm256_permutevar_pd(UMidLeft, swapFirstLast),
                                              rightBlend) - UMidCentre;
                __m256d PLeftShift = _mm256_blendv_pd(_mm256_permutevar_pd(PMidCentre, rightShift),
                                                      _mm256_permutevar_pd(PMidLeft, swapFirstLast),
                                                      rightBlend);
                __m256d PDiagShift = _mm256_blendv_pd(_mm256_permutevar_pd(PLowCentre, rightShift),
                                                      _mm256_permutevar_pd(PLowLeft, swapFirstLast),
                                                      rightBlend);
                __m256d P2 = PDiagShift + PLeftShift;
                __m256d avgx = (U1 * P1 + U2 * P2) * v_hxsqr;

                __m256d U3 = UUpCentre - UMidCentre;
                __m256d P3 = PLeftShift + PMidCentre;
                __m256d U4 = ULowCentre - UMidCentre;
                __m256d P4 = PDiagShift + PLowCentre;
                __m256d avgy = (U3 * P3 + U4 * P4) * v_hysqr;

                __m256d res = 2 * UMidCentre - UPrevVec[x] + v_tausqr * (avgx + avgy + phi());

                UNextVec[x] = res;
                if (x * vectorSize <= Sx && Sx < (x + 1) * vectorSize && y == Sy) {
                    printVec("Unextvec", UNextVec[x]);
                }

                UMax = _mm256_max_pd(UMax, res);

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
            if (y - 2 <= Sy && Sy < y + 2 && n == 1) {
                std::cout << "Unext: ";
                for (int i = 0; i < Nx; ++i) {
                    std::cout << UNext[i];
                }
                std::cout << std::endl;
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
    auto *max = (double *) (&UMax);
    for (int i = 0; i < vectorSize; i++) {
        umax = std::max(umax, max[i]);
    }
    std::cout << " Max: " << umax << std::endl;

    saveBinary();
    saveToFile();
    free(UCurrent);
    free(UPrev);
    free(P);
}

__m256d WaveSim::phi() {
    __m256d result = _mm256_set1_pd(0);
    if (x * vectorSize <= Sx && Sx < (x + 1) * vectorSize && y == Sy) {
        double part = dpi * (n * tau - t0);
        double ex = exp(-part * part * lambda);
        double sine = sin(part);
        double res = ex * sine;

        int pos = Sx - x * vectorSize;
        result[pos] = res;
        //printVec("phi", result);
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

