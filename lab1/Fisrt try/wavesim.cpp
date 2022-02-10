#include "wavesim.h"

WaveSim::WaveSim(WaveParams *params) : params(*params) {
    int Nx = params->getNx();
    int Ny = params->getNy();

    hx = (Xb - Xa) / Nx;
    hy = (Yb - Ya) / Ny;
    tau = Nx <= 1000 && Ny <= 1000 ? 0.01 : 0.001;

    UNext = new double *[Ny];
    UCurrent = new double *[Ny];
    P = new double *[Ny];
    for (y = 0; y < Ny; ++y) {
        UNext[y] = new double[Nx];
        UCurrent[y] = new double[Nx];
        P[y] = new double[Nx];
    }

    UPrev = UCurrent;
}

void WaveSim::init() {
    int Nx = params.getNx();
    int Ny = params.getNy();
    for (y = 0; y < Ny; ++y) {
        for (x = 0; x < Nx; ++x) {
            if (x < Nx / 2) {
                P[y][x] = 0.1 * 0.1;
            } else {
                P[y][x] = 0.2 * 0.2;
            }
            UNext[y][x] = 0;
            UCurrent[y][x] = 0;
        }
    }
}

void WaveSim::Run() {
    init();
    params.PrintInfo();

    int Nt = params.getNt();
    int Nx = params.getNx();
    int Ny = params.getNy();
    double UMax = 0, UMaxPrev = 0;

    time_point<high_resolution_clock> simStart = high_resolution_clock::now();
    for (n = 0; n < Nt; ++n) {


        time_point<high_resolution_clock> stepStart = high_resolution_clock::now();
        for (y = 0; y < Ny; ++y) {
            for (x = 0; x < Nx; ++x) {
                double avgx = ((UCurrent[y][x + 1] - UCurrent[y][x]) * (P[y - 1][x] + P[y][x]) +
                               (UCurrent[y][x - 1] - UCurrent[y][x]) * (P[y - 1][x - 1] + P[y][x - 1])) / (2 * hx * hx);
                double avgy = ((UCurrent[y + 1][x] - UCurrent[y][x]) * (P[y][x - 1] + P[y][x]) +
                               (UCurrent[y - 1][x] - UCurrent[y][x]) * (P[y - 1][x - 1] + P[y - 1][x])) / (2 * hy * hy);
                double result = 2 * UCurrent[y][x] - UPrev[y][x] + tau * tau * round(phi() + avgx + avgy);
                UNext[y][x] = result;
                if (result > UMax) {
                    UMax = std::abs(result);
                }
            }
        }
        time_point<high_resolution_clock> stepEnd = high_resolution_clock::now();
        nanoseconds stepTime = duration_cast<nanoseconds>(stepEnd - stepStart);
        std::cout << n << ": " << stepTime.count() / 1e9 << std::endl;

        if (UMax > UMaxPrev) {
            std::cout << "\tUMax changed from " << UMaxPrev << " to " << UMax << " diff:" << UMax - UMaxPrev
                      << std::endl;
            UMaxPrev = UMax;
        }
        std::swap(UCurrent, UPrev);
        UNext = UPrev;
    }
    time_point<high_resolution_clock> simEnd = high_resolution_clock::now();
    nanoseconds simTime = duration_cast<nanoseconds>(simEnd - simStart);
    std::cout << std::endl << "Simulation time: " << simTime.count() / 1e9 << std::endl;

    saveToFile();
}

double WaveSim::phi() {
    double result = 0;
    int Sx = params.getSx();
    int Sy = params.getSy();
    if (x == Sx && y == Sy) {
        double part = 2 * pi * f0 * (n * tau - t0);
        double ex = exp(-(part * part) / (lambda * lambda));
        double sine = sin(part);
        result = ex * sine;
    }
    return result;
}

void WaveSim::saveToFile() {
    int Nx = params.getNx();
    int Ny = params.getNy();
    std::fstream of(params.getFilename(), std::ios::out | std::ios::binary);
    if (!of.is_open()) {
        throw fileException(params.getFilename() + " cannot be open");
    }

    of.write(reinterpret_cast<char *>(&UCurrent), (Ny * Nx) * sizeof(double));

    of.close();
}

