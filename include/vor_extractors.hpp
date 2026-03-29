#ifndef VOR_EXTRACTORS_HPP
#define VOR_EXTRACTORS_HPP

#include <vector>
#include <array>

namespace VOR {
namespace Extractors {

/**
 * @brief Simple Green channel extraction
 */
std::vector<double> green(const std::vector<double>& greenValues, double fs);

/**
 * @brief Standard CHROM (Chrominance-based) algorithm
 */
std::vector<double> chrom(const std::vector<std::array<double, 3>>& rgb, double fs);

/**
 * @brief Standard POS (Plane-Orthogonal-to-Skin) algorithm
 */
std::vector<double> pos(const std::vector<std::array<double, 3>>& rgb, double fs);

} // namespace Extractors
} // namespace VOR

#endif // VOR_EXTRACTORS_HPP
