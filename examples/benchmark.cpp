#include "VOREngine.hpp"
#include <iostream>
#include <vector>
#include <array>
#include <map>

/**
 * VOR-rPPG Basic Library Benchmark
 * This example demonstrates simple library initialization and 
 * provides a mock signal processing loop for validation.
 */

int main() {
    std::cout << "--- VOR-rPPG Benchmark ---" << std::endl;

    // Initialize VOR Engine
    VOR::VOREngine engine;
    
    // Mock data: 1 ROI, 15 Hz, 12 seconds
    const double fs = 15.0;
    const int n = static_cast<int>(fs * 12.0);
    
    std::vector<double> times(n);
    std::vector<std::array<double, 3>> rgb(n);
    
    // Generate a mock 60 BPM signal (1 Hz)
    for (int i = 0; i < n; ++i) {
        times[i] = i / fs;
        double pulse = 0.5 * std::sin(2.0 * M_PI * 1.0 * times[i]);
        rgb[i] = {120.0 + pulse, 130.0 + pulse * 1.2, 110.0 + pulse * 0.8};
    }

    std::map<std::string, std::vector<std::array<double, 3>>> signals;
    signals["face"] = rgb;

    // Run VOR compute
    std::cout << "Running VOR compute on mock signal..." << std::endl;
    auto result = engine.compute(signals, times, fs);

    // Results
    if (result.valid) {
        std::cout << "Result: " << result.bpm << " BPM" << std::endl;
        std::cout << "Confidence: " << result.calibratedConf << std::endl;
        std::cout << "Trust Decision: " << (int)result.decision << std::endl;
    } else {
        std::cout << "Abstained: Insufficient confidence" << std::endl;
    }

    std::cout << "Benchmark complete." << std::endl;
    return 0;
}
