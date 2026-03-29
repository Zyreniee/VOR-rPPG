#ifndef VOR_MATH_HPP
#define VOR_MATH_HPP

#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

namespace VOR {
namespace Math {

// ============================================================================
// Statistical Utilities
// ============================================================================

double mean(const std::vector<double>& v);
double stdDev(const std::vector<double>& v);
double median(std::vector<double> v);
double mad(const std::vector<double>& v);

template <typename T>
T clip(T val, T minVal, T maxVal) {
    return std::max(minVal, std::min(maxVal, val));
}

// ============================================================================
// Signal Processing Primitives
// ============================================================================

std::vector<double> detrend(const std::vector<double>& signal);
std::vector<double> normalize(const std::vector<double>& signal);
double estimateFs(const std::vector<double>& times, double defaultFs);

struct ResampleResult {
    std::vector<double> values;
    std::vector<double> times;
};

ResampleResult resample(const std::vector<double>& signal, 
                       const std::vector<double>& times, 
                       double targetFs);

// ============================================================================
// Spectral Analysis (Welch)
// ============================================================================

struct PSDResult {
    std::vector<double> frequencies;
    std::vector<double> power;
    double df = 0.0;
    bool valid = false;
};

PSDResult welch(const std::vector<double>& signal, double fs, 
                int nperseg = 0, int noverlap = 0);

// ============================================================================
// Peak Detection
// ============================================================================

struct PeakInfo {
    int index = 0;
    double value = 0.0;
    double time = 0.0;
};

std::vector<PeakInfo> findPeaks(const std::vector<double>& signal, double fs, 
                                int minDistance = 0, double minProminence = 0.0);

std::vector<PeakInfo> filterPeaksMAD(const std::vector<PeakInfo>& peaks, 
                                     double madThreshold = 3.0);

} // namespace Math
} // namespace VOR

#endif // VOR_MATH_HPP
