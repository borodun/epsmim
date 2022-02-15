#include "wavesim.h"

WaveSim::WaveSim(WaveParams *params) : params(*params) {
    int Nx = params->getNx();
    int Ny = params->getNy();

    Sx = params->getSx();
    Sy = params->getSy();

    double hx = (Xb - Xa) / Nx;
    double hy = (Yb - Ya) / Ny;
    hxsqr = 1 / (2 * hx * hx);
    hysqr = 1 / (2 * hy * hy);
    tau = Nx <= 1000 && Ny <= 1000 ? 0.01 : 0.001;

    UCurrent = new double[Ny * Nx];
    UPrev = new double[Ny * Nx];
    P = new double[Ny * Nx];

    UNext = UPrev;
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
    double UMax = 0, UMaxPrev = 0;
    int stepOfMax = 0;

    time_point<high_resolution_clock> simStart = high_resolution_clock::now();
    for (n = 0; n < Nt; ++n) {

        //time_point<high_resolution_clock> stepStart = high_resolution_clock::now();
        for (y = 1; y < Ny - 1; ++y) {
            for (x = 1; x < Nx - 1; ++x) {
                double avgx =
                        ((UCurrent[y * Nx + x + 1] - UCurrent[y * Nx + x]) *
                         (P[(y - 1) * Nx + x] + P[y * Nx + x]) +
                         (UCurrent[y * Nx + x - 1] - UCurrent[y * Nx + x]) *
                         (P[(y - 1) * Nx + x - 1] + P[y * Nx + x - 1])) * hxsqr;
                double avgy =
                        ((UCurrent[(y + 1) * Nx + x] - UCurrent[y * Nx + x]) *
                         (P[y * Nx + x - 1] + P[y * Nx + x]) +
                         (UCurrent[(y - 1) * Nx + x] - UCurrent[y * Nx + x]) *
                         (P[(y - 1) * Nx + x - 1] + P[(y - 1) * Nx + x])) * hysqr;
                double result = 2 * UCurrent[y * Nx + x] - UPrev[y * Nx + x] + tau * tau * (phi() + avgx + avgy);
                UNext[y * Nx + x] = result;
                if (result > UMax) {
                    UMax = std::abs(result);
                }
            }
        }
        /*time_point<high_resolution_clock> stepEnd = high_resolution_clock::now();
        nanoseconds stepTime = duration_cast<nanoseconds>(stepEnd - stepStart);
        nanoseconds totalTime = duration_cast<nanoseconds>(stepEnd - simStart);
        std::cout << n << ": " << stepTime.count() / 1e9 << "s (" << totalTime.count() / 1e9 << "s)" << std::endl;*/

        if (UMax > UMaxPrev) {
            //std::cout << "\tUMax changed to " << UMax << " diff: +" << UMax - UMaxPrev << std::endl;
            UMaxPrev = UMax;
            stepOfMax = n;
        }
        std::swap(UCurrent, UPrev);
        UNext = UPrev;
    }
    time_point<high_resolution_clock> simEnd = high_resolution_clock::now();
    nanoseconds simTime = duration_cast<nanoseconds>(simEnd - simStart);
    std::cout << std::endl << "Time: " << simTime.count() / 1e9 << "s" << std::endl;
    std::cout << "Max : " << UMax << " (since step " << stepOfMax << ")" << std::endl;

    saveBinary();
    //saveToFile();
}

double WaveSim::phi() const {
    double result = 0;
    if (x == Sx && y == Sy) {
        double part = dpi * (n * tau - t0);
        double ex = exp(-(part * part)) * lambda;
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

