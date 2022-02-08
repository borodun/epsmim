#include <iostream>

int main(int argc, char *argv[]) {
    if (argc < 3) {
        std::cout << "Usage: " << argv[0] << "Nx Ny" << std::endl;
        return 1;
    }

    const int Nx = std::stoi(argv[1]);
    const int Ny = std::stoi(argv[2]);
    if (Nx < 0 || Ny < 0) {
        std::cout << "Dimensions cannot be negative" << std::endl;
        return 1;
    }

    const int Xa = 0, Xb = 4, Ya = 0, Yb = 4;
    const double hx = (Xb - Xa) / Nx;
    const double hy = (Yb - Ya) / Ny;

    const double tau = Nx <= 1000 && Ny <= 1000 ? 0.01 : 0.001;

    double UNext[Nx][Ny];
    double UPrev[Nx][Ny];
    double P[Nx][Ny];
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            if (j < Nx / 2) {
                P[i][j] = 0.1 * 0.1;
            } else {
                P[i][j] = 0.2 * 0.2;
            }
        }
    }




    return 0;
}
