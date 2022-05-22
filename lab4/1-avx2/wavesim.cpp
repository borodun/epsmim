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
    for (int y = 0; y < Ny; ++y) {
        for (int x = 0; x < Nx; ++x) {
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
    int threads = params.getThreads();
    int chunk = Ny / threads;

    printf("Threads: %d, chunk: %d\n", threads, chunk);

    double UMax = 0;
    __m256d VUMax = _mm256_set1_pd(0);

    time_point<high_resolution_clock> simStart = high_resolution_clock::now();

#pragma omp parallel proc_bind(close) num_threads(threads)
    for (int n = 0; n < Nt; ++n) {
//        int thread = omp_get_thread_num();
//        printf("%d - Iter: %d\n", thread, n);
        double localMax = UMax;
        __m256d localVUMax = VUMax;

        auto *UUp = (__m256d *) (UCurrent);
        auto *UMid = (__m256d *) (UCurrent + Nx);
        auto *ULow = (__m256d *) (UCurrent + 2 * Nx);

        auto *UPrevVec = (__m256d *) (UPrev + Nx);
        auto *UNextVec = (__m256d *) (UNext + Nx);

        auto *PMid = (__m256d *) (P + Nx);
        auto *PLow = (__m256d *) (P + 2 * Nx);

//        time_point<high_resolution_clock> stepStart = high_resolution_clock::now();
#pragma omp for schedule(static, chunk)
        for (int y = 1; y < Ny - 1; ++y) {
//            printf("%d - Y: %d\n", thread, y);
            __m256d UUpLeft = UUp[y - 1];
            __m256d UUpCentre = UUp[y];
            __m256d UUpRight = UUp[y + 1];

            __m256d UMidLeft = UMid[y - 1];
            __m256d UMidCentre = UMid[y];
            __m256d UMidRight = UMid[y + 1];

            __m256d ULowLeft = ULow[y - 1];
            __m256d ULowCentre = ULow[y];
            __m256d ULowRight = ULow[y + 1];

            __m256d PMidLeft = PMid[y - 1];
            __m256d PMidCentre = PMid[y];
            __m256d PMidRight = PMid[y + 1];

            __m256d PLowLeft = PLow[y - 1];
            __m256d PLowCentre = PLow[y];
            __m256d PLowRight = PLow[y + 1];

            int row = y * Nx;
            int prevRow = (y - 1) * Nx;
            int nextRow = (y + 1) * Nx;

            for (int x = 1; x < vectorSize; ++x) {
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
                double result = 2 * currentPos - UPrev[pos] + tausqr * (phi(x, y, n) + avgx + avgy);
                UNext[pos] = result;
                if (result > localMax) {
                    localMax = std::abs(result);
                }
            }

            for (int x = 1; x < Nx / vectorSize; ++x) {
                int leftShift = 0b10010000;
                int rightShift = 0b11111001;
                int swapFirstLast = 0b00100111;
                int leftBlend = 0b0001;
                int rightBlend = 0b1000;

                __m256d U1 = _mm256_blend_pd(_mm256_permute4x64_pd(UMidCentre, leftShift),
                                             _mm256_permute4x64_pd(UMidRight, swapFirstLast),
                                             leftBlend) - UMidCentre;
                __m256d P1 = PLowCentre + PMidCentre;
                __m256d U2 = _mm256_blend_pd(_mm256_permute4x64_pd(UMidCentre, rightShift),
                                             _mm256_permute4x64_pd(UMidLeft, swapFirstLast),
                                             rightBlend) - UMidCentre;
                __m256d PLeftShift = _mm256_blend_pd(_mm256_permute4x64_pd(PMidCentre, rightShift),
                                                     _mm256_permute4x64_pd(PMidLeft, swapFirstLast),
                                                     rightBlend);
                __m256d PDiagShift = _mm256_blend_pd(_mm256_permute4x64_pd(PLowCentre, rightShift),
                                                     _mm256_permute4x64_pd(PLowLeft, swapFirstLast),
                                                     rightBlend);
                __m256d P2 = PDiagShift + PLeftShift;
                __m256d avgx = (U1 * P1 + U2 * P2) * v_hxsqr;

                __m256d U3 = UUpCentre - UMidCentre;
                __m256d P3 = PLeftShift + PMidCentre;
                __m256d U4 = ULowCentre - UMidCentre;
                __m256d P4 = PDiagShift + PLowCentre;
                __m256d avgy = (U3 * P3 + U4 * P4) * v_hysqr;

                __m256d res = 2 * UMidCentre - UPrevVec[x] + v_tausqr * (avgx + avgy + vphi(x, y, n));
                UNextVec[x] = res;
                localVUMax = _mm256_max_pd(localVUMax, res);

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

            for (int x = Nx - 1 - vectorSize; x < Nx - 1; ++x) {
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
                double result = 2 * currentPos - UPrev[pos] + tausqr * (phi(x, y, n) + avgx + avgy);
                UNext[pos] = result;
                if (result > localMax) {
                    localMax = std::abs(result);
                }
            }

            ULow = UMid;
            UMid = UUp;
            UUp += Nx / vectorSize;

            UPrevVec += Nx / vectorSize;

            PLow = PMid;
            PMid += Nx / vectorSize;

#pragma omp critical
            {
                VUMax = _mm256_max_pd(localVUMax, VUMax);
                if (localMax > UMax) {
                    UMax = localMax;
                }
            }
        }

#pragma omp single
        {
//            time_point<high_resolution_clock> stepEnd = high_resolution_clock::now();
//            nanoseconds stepTime = duration_cast<nanoseconds>(stepEnd - stepStart);
//            nanoseconds totalTime = duration_cast<nanoseconds>(stepEnd - simStart);
//            std::cout << thread << " - " << n << ": " << stepTime.count() / 1e9 << "s (" << totalTime.count() / 1e9 << "s)" << std::endl;

            std::swap(UCurrent, UPrev);
            UNext = UPrev;
        }
#pragma omp barrier
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

__m256d WaveSim::vphi(int x, int y, int n) {
    __m256d result = _mm256_set1_pd(0);
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

double WaveSim::phi(int x, int y, int n) {
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

    for (int y = 0; y < Ny; ++y) {
        for (int x = 0; x < Nx; ++x) {
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
    for (int y = 0; y < Ny; ++y) {
        for (int x = 0; x < Nx; ++x) {
            of << UCurrent[y * Nx + x] << " ";
        }
        of << std::endl;
    }

    of.close();
}

void WaveSim::printArray(double *arr) {
    int Nx = params.getNx();
    int Ny = params.getNy();
    for (int y = 0; y < Ny; ++y) {
        for (int x = 0; x < Nx; ++x) {
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

