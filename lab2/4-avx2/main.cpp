#include <iostream>
#include "wavesim.h"

int main(int argc, char *argv[]) {
    std::cout << "Wave simulation" << std::endl;
    try {
        auto params = new WaveParams(argc, argv);
        auto sim = new WaveSim(params);
        sim->Run();
        delete params;
        delete sim;
    } catch (WaveSimException &e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    return 0;
}
