#include <iostream>
#include "wavesim.h"

int main(int argc, char *argv[]) {
    try {
        auto params = new WaveParams(argc, argv);
        auto sim = new WaveSim(params);
        sim->Run();
    } catch (WaveSimException &e) {
        std::cerr << e.what() << std::endl;
    }
    return 0;
}
