#ifndef VOR_ENGINE_HPP
#define VOR_ENGINE_HPP

#include <array>
#include <cmath>
#include <deque>
#include <string>
#include <vector>
#include <map>

/**
 * VOR: Variance-Optimized rPPG
 * 
 * Formal research implementation of the VOR algorithm.
 * Fuses multi-modal rPPG signals using risk-aware spectral logic.
 */

namespace VOR {

// ============================================================================
// Configuration
// ============================================================================

struct VORConfig {
  // Windowing
  double windowSeconds = 12.0;
  double minWindowSeconds = 3.5;

  // Frequency band (BPM range)
  double fMinHz = 0.7;  // ~42 BPM
  double fMaxHz = 3.5;  // ~210 BPM

  // Spectral Analysis
  int topKPeaks = 5;
  double minPeakSNR = 1.5;

  // Cluster Fusion
  double clusterRadiusBPM = 6.0;
  int minClusterSize = 2;

  // Decision Thresholds
  double acceptThreshold = 0.60;
  double degradedThreshold = 0.35;

  // Risk Penalty Parameters
  double boundaryPenaltyWidth = 0.15;
  double harmonicPenaltyRatio = 0.3;

  // Diagnostics
  bool enableConsoleLog = true;
};

// ============================================================================
// Data Structures
// ============================================================================

struct ROISignalBuffer {
  std::string roiName;
  std::vector<double> times;
  std::vector<double> green;
  std::vector<double> red;
  std::vector<double> blue;

  double fsEff = 0.0;
  std::vector<double> resGreen;
  std::vector<double> resRed;
  std::vector<double> resBlue;
  std::vector<double> resTimes;
  std::vector<std::array<double, 3>> resRGB;

  bool valid = false;
  int sampleCount = 0;
};

struct ExtractorOutput {
  std::string extractorName;
  std::string roiName;
  std::vector<double> signal;
  double snr = 0.0;
  bool valid = false;
};

struct SpectralPeak {
  double hz = 0.0;
  double bpm = 0.0;
  double power = 0.0;
  double snr = 0.0;
  int rank = 0;
  double rho_b = 0.0; // Boundary risk
  double rho_h = 0.0; // Harmonic risk
};

struct VORCandidate {
  double hz = 0.0;
  double bpm = 0.0;
  double power = 0.0;
  double snr = 0.0;
  double confidence = 0.0;
  std::string source;
  std::string roi;
  int peakRank = 0;
  double boundaryRisk = 0.0;
  double harmonicRisk = 0.0;
};

struct VORCluster {
  double centerBPM = 0.0;
  double spreadBPM = 0.0;
  int memberCount = 0;
  double totalWeight = 0.0;
  double score = 0.0;
  double avgSNR = 0.0;
  bool hasPSD = false;
  bool hasPeaks = false;
  std::vector<int> memberIdx;
};

enum class VORDecision { ACCEPT, HOLD, DEGRADED, ABSTAIN };

struct VORResult {
  double bpm = 0.0;
  double confidence = 0.0;
  double calibratedConf = 0.0;
  VORDecision decision = VORDecision::ABSTAIN;
  double fsEff = 0.0;
  int nCandidates = 0;
  int nClusters = 0;
  int winningClusterSize = 0;
  double winningClusterSpread = 0.0;
  double psdResolutionBPM = 0.0;
  double temporalT = 0.0;
  double methodM = 0.0;
  double resolutionR = 0.0;
  double trustedHR = 0.0;
  int trustAge = 0;
  std::string reasonCode;
  bool valid = false;

  // Extractor Diagnostics
  double chromBPM = 0.0, chromSNR = 0.0;
  double posBPM = 0.0, posSNR = 0.0;
  double greenBPM = 0.0, greenSNR = 0.0;
  std::vector<double> top1BPMs;
};

// ============================================================================
// VOREngine
// ============================================================================

class VOREngine {
public:
  VOREngine();

  /**
   * @brief Process multi-ROI RGB signals to estimate pulse
   * @param multiRgbSignals Map of ROI name -> raw RGB samples
   * @param timestamps Clock times for each sample
   * @param nominalFs Target sampling frequency
   */
  VORResult compute(const std::map<std::string, std::vector<std::array<double, 3>>> &multiRgbSignals,
                    const std::vector<double> &timestamps, double nominalFs);

  VORConfig &config() { return m_config; }

private:
  VORConfig m_config;
  int m_cycleCount = 0;
  std::deque<double> m_recentWinners;
  double m_trustedHR = 0.0;
  double m_trustedConf = 0.0;
  int m_trustAge = 0;
  double m_lastPsdResolution = 0.0;

  static constexpr int MAX_RECENT = 5;
  static constexpr int MAX_HOLD_CYCLES = 10;

  // Pipeline Logic
  ROISignalBuffer prepareROI(const std::string &roiName, const std::vector<std::array<double, 3>> &rgb, const std::vector<double> &times, double fs);
  std::vector<ExtractorOutput> runExtractors(const ROISignalBuffer &roi);
  std::vector<SpectralPeak> analyzeSpectrum(const std::vector<double> &signal, double fs, double fMin, double fMax, int topK);
  std::vector<VORCandidate> buildCandidates(const std::vector<ExtractorOutput> &outputs, double fs);
  std::vector<VORCluster> clusterCandidates(const std::vector<VORCandidate> &candidates);
  VORResult decideTrust(const std::vector<VORCluster> &clusters, const std::vector<VORCandidate> &candidates, double fsEff);
  void applyTemporalState(VORResult &result);
};

} // namespace VOR

#endif // VOR_ENGINE_HPP
