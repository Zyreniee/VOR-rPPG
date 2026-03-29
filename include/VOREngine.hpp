#ifndef VOR_ENGINE_HPP
#define VOR_ENGINE_HPP

#include <array>
#include <cmath>
#include <deque>
#include <string>
#include <vector>
#include <map>

// ============================================================================
// VOR: Variance-Optimized rPPG
// ============================================================================

namespace VOR {

// ============================================================================
// Configuration
// ============================================================================

struct VORConfig {
  // Window
  double windowSeconds = 12.0;
  double minWindowSeconds = 3.5;

  // Frequency band
  double fMinHz = 0.7;  // 42 BPM
  double fMaxHz = 3.5;  // 210 BPM

  // PeakAnalyzer
  int topKPeaks = 5;
  double minPeakSNR = 1.5;

  // ClusterFusion
  double clusterRadiusBPM = 6.0;
  int minClusterSize = 2;

  // TrustEngine thresholds
  double acceptThreshold = 0.60;
  double degradedThreshold = 0.35;

  // Risk penalties (Stage 2 — guarded by flags)
  double boundaryPenaltyWidth = 0.15;
  double harmonicPenaltyRatio = 0.3;

  // Feature flags
  bool enableVORProjection = false;
  bool enableBoundaryRisk = false;
  bool enableHarmonicRisk = false;
  bool enableFlickerPenalty = false;
  bool enableMotionCalib = false;

  // Logging
  bool enableCSVLog = true;
  bool enableConsoleLog = true;
};

// ============================================================================
// Data Structures
// ============================================================================

struct ROISignalBuffer {
  std::string roiName;

  // Raw (pre-slice, pre-resample)
  std::vector<double> times;
  std::vector<double> green;
  std::vector<double> red;
  std::vector<double> blue;

  // After 12s slice + resample
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
  std::string extractorName; // "green", "chrom", "pos", "vor"
  std::string roiName;       // "face"
  std::vector<double> signal;
  double snr = 0.0;
  bool valid = false;
};

struct SpectralPeak {
  double hz = 0.0;
  double bpm = 0.0;
  double power = 0.0;
  double snr = 0.0;
  int rank = 0; // 1-based

  // Risk annotations (Section 3.5 of VOR spec)
  double rho_b = 0.0; // Boundary risk: proximity to band edges [0,1]
  double rho_h = 0.0; // Harmonic risk: likelihood of being a harmonic [0,1]
};

struct VORCandidate {
  double hz = 0.0;
  double bpm = 0.0;
  double power = 0.0;
  double snr = 0.0;
  double confidence = 0.0;

  std::string source; // "chrom_psd", "pos_peaks", etc.
  std::string roi;
  int peakRank = 0;

  // Risk (Stage 2)
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
  double confidence = 0.0;       // Raw per-cycle confidence
  double calibratedConf = 0.0;   // After T×M×R adjustment
  VORDecision decision = VORDecision::ABSTAIN;

  // Diagnostics
  double fsEff = 0.0;
  int nCandidates = 0;
  int nClusters = 0;
  int winningClusterSize = 0;
  double winningClusterSpread = 0.0;
  double psdResolutionBPM = 0.0;

  // Temporal state diagnostics
  double temporalT = 0.0;        // Temporal consistency [0,1]
  double methodM = 0.0;          // Method agreement [0,1]
  double resolutionR = 0.0;      // Resolution quality [0,1]
  double trustedHR = 0.0;        // Current trusted HR
  int trustAge = 0;              // Cycles since trusted HR updated
  std::string reasonCode;        // ACCEPT_STABLE, HOLD_VOLATILE, etc.

  // Per-extractor top-1
  double chromBPM = 0.0, chromSNR = 0.0;
  double posBPM = 0.0, posSNR = 0.0;
  double greenBPM = 0.0, greenSNR = 0.0;
  std::vector<double> top1BPMs; // Store all top-1 candidates to compute Method Agreement M across ROIs

  bool valid = false;
};

// ============================================================================
// VOREngine
// ============================================================================

class VOREngine {
public:
  VOREngine();

  /// Main compute — call from RPPGEngine::computeHRHRVRR()
  VORResult compute(const std::map<std::string, std::vector<std::array<double, 3>>> &multiRgbSignals,
                    const std::vector<double> &faceTimes, double nominalFs);

  VORConfig &config() { return m_config; }
  const VORConfig &config() const { return m_config; }

private:
  VORConfig m_config;
  int m_cycleCount = 0;

  // ---- Temporal state (Patch 2) ----
  std::deque<double> m_recentWinners;  // Last N raw winning BPMs
  double m_trustedHR = 0.0;            // Last HR with high calibrated confidence
  double m_trustedConf = 0.0;          // Confidence of trusted HR
  int m_trustAge = 0;                  // Cycles since m_trustedHR was updated
  double m_lastPsdResolution = 0.0;    // BPM resolution of last PSD
  static constexpr int MAX_RECENT = 5;
  static constexpr int MAX_HOLD_CYCLES = 10;

  // ---- Pipeline stages ----

  /// Slice to last 12s, per-channel resample, build ROISignalBuffer
  ROISignalBuffer
  prepareROI(const std::string &roiName,
             const std::vector<std::array<double, 3>> &rgb,
             const std::vector<double> &times, double nominalFs);

  /// Run Green, CHROM, POS (and optionally VOR) extractors
  std::vector<ExtractorOutput> runExtractors(const ROISignalBuffer &roi);

  /// Extract top-k spectral peaks from a signal
  std::vector<SpectralPeak> analyzeSpectrum(const std::vector<double> &signal,
                                            double fs, double fMin,
                                            double fMax, int topK);

  /// Build candidates from all extractor outputs
  std::vector<VORCandidate>
  buildCandidates(const std::vector<ExtractorOutput> &outputs, double fs);

  /// Cluster candidates by BPM proximity
  std::vector<VORCluster>
  clusterCandidates(const std::vector<VORCandidate> &candidates);

  /// Decide accept/degraded/abstain from clusters (per-cycle)
  VORResult decideTrust(const std::vector<VORCluster> &clusters,
                        const std::vector<VORCandidate> &candidates,
                        double fsEff);

  /// Apply temporal state model: recalibrate confidence, decide HOLD
  void applyTemporalState(VORResult &result);

  // ---- Logging ----
  void logCSV(const VORResult &result,
              const std::vector<VORCandidate> &candidates,
              const std::vector<VORCluster> &clusters);
  void logConsole(const VORResult &result);
  void logState(const VORResult &result);

  bool m_csvHeaderWritten = false;
};

} // namespace VOR

#endif // VOR_ENGINE_HPP
