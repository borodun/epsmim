#include <string>
#include "waveparams.h"

WaveParams::WaveParams(int argc, char *argv[]) {
    if (argc < 8) {
        throw badArgsException(std::string("Usage: ") + argv[0] + "Nx Ny Nt Sx Sy m filename");
    }

    Nx = std::stoi(argv[1]);
    Ny = std::stoi(argv[2]);
    Nt = std::stoi(argv[3]);
    Sx = std::stoi(argv[4]);
    Sy = std::stoi(argv[5]);
    m = std::stoi(argv[6]);
    filename = argv[7];

    if (Nx < 0 || Ny < 0) {
        throw badArgsException("Dimensions should be positive");
    }
    if (Nt < 0) {
        throw badArgsException("Step count should be positive");
    }
    if (Sx < 1 || Sy < 1 || Sx > Nx - 2 || Sy > Ny - 2) {
        throw badArgsException("Source point should be in [0; N-1]");
    }
    if (m < 1) {
        throw badArgsException("m should be in [1; Nt-1]");
    }
}

int WaveParams::getNx() const {
    return Nx;
}

int WaveParams::getNy() const {
    return Ny;
}

int WaveParams::getNt() const {
    return Nt;
}

int WaveParams::getSx() const {
    return Sx;
}

int WaveParams::getSy() const {
    return Sy;
}

const std::string &WaveParams::getFilename() const {
    return filename;
}

void WaveParams::PrintInfo() {
    std::cout << "Params: " << std::endl
              << "\tField : " << Nx << "x" << Ny << std::endl
              << "\tSteps : " << Nt << std::endl
              << "\tSource: (" << Sx << ", " << Sy << ")" << std::endl;
}

int WaveParams::getM() const {
    return m;
}
