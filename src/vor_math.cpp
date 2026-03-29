#include "vor_math.hpp"
#include <cmath>
#include <iostream>

namespace VOR {
namespace Math {

// ============================================================================
// Statistical Utilities
// ============================================================================

double mean(const std::vector<double> &v) {
  if (v.empty())
    return 0.0;
  return std::accumulate(v.begin(), v.end(), 0.0) /
         static_cast<double>(v.size());
}

double stdDev(const std::vector<double> &v) {
  if (v.size() < 2)
    return 0.0;
  double m = mean(v);
  double sum = 0.0;
  for (const auto &x : v) {
    sum += (x - m) * (x - m);
  }
  return std::sqrt(sum / static_cast<double>(v.size() - 1));
}

double median(std::vector<double> v) {
  if (v.empty())
    return 0.0;
  size_t n = v.size();
  std::nth_element(v.begin(), v.begin() + n / 2, v.end());
  if (n % 2 == 0) {
    double med2 = v[n / 2];
    std::nth_element(v.begin(), v.begin() + n / 2 - 1, v.end());
    return (v[n / 2 - 1] + med2) / 2.0;
  }
  return v[n / 2];
}

double mad(const std::vector<double> &v) {
  if (v.empty())
    return 0.0;
  double med = median(v);
  std::vector<double> absDeviations;
  absDeviations.reserve(v.size());
  for (const auto &x : v) {
    absDeviations.push_back(std::abs(x - med));
  }
  return median(absDeviations);
}

// ============================================================================
// Signal Processing Primitives
// ============================================================================

std::vector<double> detrend(const std::vector<double> &signal) {
  if (signal.size() < 2)
    return signal;
  int n = static_cast<int>(signal.size());
  double sumX = 0.0, sumY = 0.0, sumXY = 0.0, sumX2 = 0.0;
  for (int i = 0; i < n; ++i) {
    double x = static_cast<double>(i);
    sumX += x;
    sumY += signal[i];
    sumXY += x * signal[i];
    sumX2 += x * x;
  }
  double denom = n * sumX2 - sumX * sumX;
  if (std::abs(denom) < 1e-12)
    return signal;
  double b = (n * sumXY - sumX * sumY) / denom;
  double a = (sumY - b * sumX) / n;
  std::vector<double> result(n);
  for (int i = 0; i < n; ++i) {
    result[i] = signal[i] - (a + b * static_cast<double>(i));
  }
  return result;
}

std::vector<double> normalize(const std::vector<double> &signal) {
  if (signal.empty())
    return {};
  double m = mean(signal);
  double s = stdDev(signal);
  if (s < 1e-10) s = 1e-10;
  std::vector<double> result(signal.size());
  for (size_t i = 0; i < signal.size(); ++i) {
    result[i] = (signal[i] - m) / s;
  }
  return result;
}

double estimateFs(const std::vector<double> &times, double defaultFs) {
  if (times.size() < 6) return defaultFs;
  std::vector<double> diffs;
  diffs.reserve(times.size() - 1);
  for (size_t i = 1; i < times.size(); ++i) {
    double dt = times[i] - times[i - 1];
    if (std::isfinite(dt) && dt > 0.0) diffs.push_back(dt);
  }
  if (diffs.size() < 4) return defaultFs;
  double medDt = median(diffs);
  if (medDt <= 1e-6) return defaultFs;
  return clip(1.0 / medDt, 10.0, 90.0);
}

ResampleResult resample(const std::vector<double> &signal,
                      const std::vector<double> &times, double targetFs) {
  if (signal.empty() || times.empty() || times.size() != signal.size() || times.size() < 2)
    return {{}, {}};

  double t0 = times.front();
  double t1 = times.back();
  if (t1 <= t0) return {{}, {}};

  double dt = 1.0 / targetFs;
  std::vector<double> newTimes;
  for (double t = t0; t < t1; t += dt) newTimes.push_back(t);

  std::vector<double> newSignal(newTimes.size());
  for (size_t i = 0; i < newTimes.size(); ++i) {
    double t = newTimes[i];
    auto it = std::lower_bound(times.begin(), times.end(), t);
    size_t idx = std::distance(times.begin(), it);

    if (idx == 0) newSignal[i] = signal[0];
    else if (idx >= times.size()) newSignal[i] = signal.back();
    else {
      double t_start = times[idx - 1];
      double t_end = times[idx];
      double alpha = (t - t_start) / (t_end - t_start);
      newSignal[i] = signal[idx - 1] + alpha * (signal[idx] - signal[idx - 1]);
    }
  }
  return {newSignal, newTimes};
}

// ============================================================================
// Spectral Analysis (Simplified Welch DFT)
// ============================================================================

PSDResult welch(const std::vector<double> &signal, double fs, int nperseg,
               int noverlap) {
  PSDResult res;
  int n = static_cast<int>(signal.size());
  if (n < 16) return res;

  if (nperseg <= 0) nperseg = std::min(256, std::max(16, n / 3));
  if (noverlap <= 0) noverlap = nperseg / 2;

  int step = std::max(1, nperseg - noverlap);
  int nSegments = (n - nperseg) / step + 1;
  int nFreqs = nperseg / 2 + 1;

  res.frequencies.resize(nFreqs);
  res.power.assign(nFreqs, 0.0);
  res.df = fs / nperseg;
  for (int i = 0; i < nFreqs; ++i) res.frequencies[i] = i * res.df;

  // Hamming window
  std::vector<double> window(nperseg);
  double winSumSq = 0.0;
  for (int i = 0; i < nperseg; ++i) {
    window[i] = 0.54 - 0.46 * std::cos(2.0 * M_PI * i / (nperseg - 1));
    winSumSq += window[i] * window[i];
  }

  // Simplified DFT-based Welch for isolation (O(N^2) but portable)
  for (int seg = 0; seg < nSegments; ++seg) {
    int start = seg * step;
    for (int k = 0; k < nFreqs; ++k) {
      double real = 0.0, imag = 0.0;
      for (int t = 0; t < nperseg; ++t) {
        double angle = -2.0 * M_PI * k * t / nperseg;
        double val = signal[start + t] * window[t];
        real += val * std::cos(angle);
        imag += val * std::sin(angle);
      }
      res.power[k] += (real * real + imag * imag) / (fs * winSumSq);
    }
  }

  for (auto &p : res.power) p /= nSegments;
  res.valid = true;
  return res;
}

// ============================================================================
// Peak Detection
// ============================================================================

std::vector<PeakInfo> findPeaks(const std::vector<double> &signal, double fs,
                                int minDistance, double minProminence) {
  std::vector<PeakInfo> peaks;
  if (signal.size() < 3) return peaks;

  for (size_t i = 1; i < signal.size() - 1; ++i) {
    if (signal[i] > signal[i - 1] && signal[i] > signal[i + 1] && signal[i] > minProminence) {
      // Basic peak found
      bool distanceOk = true;
      if (minDistance > 0) {
        for (const auto& p : peaks) {
          if (static_cast<int>(i) - p.index < minDistance) {
            distanceOk = false;
            break;
          }
        }
      }
      if (distanceOk) {
        peaks.push_back({static_cast<int>(i), signal[i], i / fs});
      }
    }
  }
  return peaks;
}

std::vector<PeakInfo> filterPeaksMAD(const std::vector<PeakInfo> &peaks,
                                     double madThreshold) {
  if (peaks.size() < 3) return peaks;
  std::vector<double> intervals;
  for (size_t i = 1; i < peaks.size(); ++i) {
    intervals.push_back(peaks[i].time - peaks[i - 1].time);
  }
  double med = median(intervals);
  double m = mad(intervals);
  if (m < 1e-6) return peaks;

  std::vector<PeakInfo> filtered;
  filtered.push_back(peaks[0]);
  for (size_t i = 1; i < peaks.size(); ++i) {
    double dt = peaks[i].time - peaks[i - 1].time;
    if (std::abs(dt - med) < madThreshold * m) {
      filtered.push_back(peaks[i]);
    }
  }
  return filtered;
}

} // namespace Math
} // namespace VOR
