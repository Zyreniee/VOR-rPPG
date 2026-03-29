#include "vor_extractors.hpp"
#include "vor_math.hpp"
#include <cmath>
#include <numeric>

namespace VOR {
namespace Extractors {

// ============================================================================
// Green Channel
// ============================================================================

std::vector<double> green(const std::vector<double>& greenValues, double fs) {
    if (greenValues.size() < 24) return {};
    auto sig = Math::detrend(greenValues);
    return Math::normalize(sig);
}

// ============================================================================
// CHROM (Chrominance-based)
// ============================================================================

std::vector<double> chrom(const std::vector<std::array<double, 3>>& rgb, double fs) {
    if (rgb.size() < 24) return {};

    int n = static_cast<int>(rgb.size());
    std::vector<double> R(n), G(n), B(n);
    for (int i = 0; i < n; ++i) {
        R[i] = rgb[i][0];
        G[i] = rgb[i][1];
        B[i] = rgb[i][2];
    }

    double mR = Math::mean(R) + 1e-8;
    double mG = Math::mean(G) + 1e-8;
    double mB = Math::mean(B) + 1e-8;

    std::vector<double> Rn(n), Gn(n), Bn(n);
    for (int i = 0; i < n; ++i) {
        Rn[i] = R[i] / mR - 1.0;
        Gn[i] = G[i] / mG - 1.0;
        Bn[i] = B[i] / mB - 1.0;
    }

    std::vector<double> Xs(n), Ys(n);
    for (int i = 0; i < n; ++i) {
        Xs[i] = 3.0 * Rn[i] - 2.0 * Gn[i];
        Ys[i] = 1.5 * Rn[i] + Gn[i] - 1.5 * Bn[i];
    }

    double alpha = Math::stdDev(Xs) / (Math::stdDev(Ys) + 1e-8);
    std::vector<double> S(n);
    for (int i = 0; i < n; ++i) {
        S[i] = Xs[i] - alpha * Ys[i];
    }

    return Math::normalize(Math::detrend(S));
}

// ============================================================================
// POS (Plane-Orthogonal-to-Skin)
// ============================================================================

std::vector<double> pos(const std::vector<std::array<double, 3>>& rgb, double fs) {
    if (rgb.size() < 24) return {};

    int n = static_cast<int>(rgb.size());
    std::vector<double> R(n), G(n), B(n);
    for (int i = 0; i < n; ++i) {
        R[i] = rgb[i][0];
        G[i] = rgb[i][1];
        B[i] = rgb[i][2];
    }

    double mR = Math::mean(R) + 1e-8;
    double mG = Math::mean(G) + 1e-8;
    double mB = Math::mean(B) + 1e-8;

    std::vector<double> r(n), g(n), b(n);
    for (int i = 0; i < n; ++i) {
        r[i] = R[i] / mR - 1.0;
        g[i] = G[i] / mG - 1.0;
        b[i] = B[i] / mB - 1.0;
    }

    std::vector<double> X(n), Y(n);
    for (int i = 0; i < n; ++i) {
        X[i] = g[i] - b[i];
        Y[i] = -2.0 * r[i] + g[i] + b[i];
    }

    double alpha = Math::stdDev(X) / (Math::stdDev(Y) + 1e-8);
    std::vector<double> S(n);
    for (int i = 0; i < n; ++i) {
        S[i] = X[i] + alpha * Y[i];
    }

    return Math::normalize(Math::detrend(S));
}

} // namespace Extractors
} // namespace VOR
