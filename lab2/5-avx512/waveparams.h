#ifndef LAB1_WAVEPARAMS_H
#define LAB1_WAVEPARAMS_H

#include "exceptions.h"
#include <iostream>

class WaveParams {
private:
    int Nx, Ny;
    int Nt;
    int Sx, Sy;
    std::string filename;
public:
    WaveParams(int argc, char *argv[]);

    [[nodiscard]] int getNx() const;

    [[nodiscard]] int getNy() const;

    [[nodiscard]] int getNt() const;

    [[nodiscard]] int getSx() const;

    [[nodiscard]] int getSy() const;

    [[nodiscard]] const std::string &getFilename() const;

    void PrintInfo();
};


#endif
