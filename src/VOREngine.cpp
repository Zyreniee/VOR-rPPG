#include "VOREngine.hpp"
#include "vor_math.hpp"
#include "vor_extractors.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>

namespace VOR {

// ============================================================================
// Constructor
// ============================================================================

VOREngine::VOREngine() {}

// ============================================================================
// Main Compute Entry Point
// ============================================================================

VORResult VOREngine::compute(const std::map<std::string, std::vector<std::array<double, 3>>> &multiRgbSignals,
                             const std::vector<double> &faceTimes,
                             double nominalFs) {
  m_cycleCount++;

  std::vector<ExtractorOutput> allOutputs;
  double fsEff = 0.0;

  for (const auto& [roiName, rgbSignal] : multiRgbSignals) {
    auto roi = prepareROI(roiName, rgbSignal, faceTimes, nominalFs);
    if (roi.valid) {
      fsEff = roi.fsEff;
      auto outputs = runExtractors(roi);
      allOutputs.insert(allOutputs.end(), outputs.begin(), outputs.end());
    } else if (m_config.enableConsoleLog) {
      std::cout << "[VOR] cycle=" << m_cycleCount
                << " SKIP ROI (" << roiName << "): insufficient data" << std::endl;
    }
  }

  if (allOutputs.empty()) {
    VORResult r;
    r.decision = VORDecision::ABSTAIN;
    r.fsEff = fsEff;
    return r;
  }

  // Build candidates from all extractors across all ROIs
  auto candidates = buildCandidates(allOutputs, fsEff);

  // Cluster candidates
  auto clusters = clusterCandidates(candidates);

  // Trust decision (per-cycle, memoryless)
  auto result = decideTrust(clusters, candidates, fsEff);

  // Save raw winning BPM before temporal state may override it
  double rawWinnerBPM = result.bpm;

  // Temporal state model — recalibrate confidence, HOLD logic
  applyTemporalState(result);

  // Update recent winners
  if (rawWinnerBPM > 0.0) {
    m_recentWinners.push_back(rawWinnerBPM);
    if (static_cast<int>(m_recentWinners.size()) > MAX_RECENT) {
      m_recentWinners.pop_front();
    }
  }

  // Fill diagnostics
  for (const auto &c : candidates) {
    if (c.peakRank == 1) {
      result.top1BPMs.push_back(c.bpm);
      if (c.source.find("chrom") != std::string::npos && result.chromBPM == 0.0) {
        result.chromBPM = c.bpm;
        result.chromSNR = c.snr;
      } else if (c.source.find("pos") != std::string::npos && result.posBPM == 0.0) {
        result.posBPM = c.bpm;
        result.posSNR = c.snr;
      } else if (c.source.find("green") != std::string::npos && result.greenBPM == 0.0) {
        result.greenBPM = c.bpm;
        result.greenSNR = c.snr;
      }
    }
  }

  return result;
}

// ============================================================================
// Stage 1: Prepare ROI
// ============================================================================

ROISignalBuffer
VOREngine::prepareROI(const std::string &roiName,
                      const std::vector<std::array<double, 3>> &rgb,
                      const std::vector<double> &times, double nominalFs) {
  ROISignalBuffer roi;
  roi.roiName = roiName;

  if (times.size() < 10 || rgb.size() != times.size()) {
    return roi;
  }

  double tEnd = times.back();
  double tStart = tEnd - m_config.windowSeconds;
  size_t sliceIdx = 0;
  for (size_t i = 0; i < times.size(); ++i) {
    if (times[i] >= tStart) {
      sliceIdx = i;
      break;
    }
  }

  roi.times.assign(times.begin() + sliceIdx, times.end());
  roi.green.resize(roi.times.size());
  roi.red.resize(roi.times.size());
  roi.blue.resize(roi.times.size());

  for (size_t i = 0; i < roi.times.size(); ++i) {
    size_t srcIdx = sliceIdx + i;
    roi.green[i] = rgb[srcIdx][1];
    roi.red[i] = rgb[srcIdx][0];
    roi.blue[i] = rgb[srcIdx][2];
  }

  double duration = roi.times.size() > 1 ? roi.times.back() - roi.times.front() : 0.0;
  if (duration < m_config.minWindowSeconds) {
    return roi;
  }

  roi.fsEff = Math::estimateFs(roi.times, nominalFs);
  roi.fsEff = Math::clip(roi.fsEff, 10.0, 60.0);

  // Per-channel resampling
  auto resR = Math::resample(roi.red, roi.times, roi.fsEff);
  auto resG = Math::resample(roi.green, roi.times, roi.fsEff);
  auto resB = Math::resample(roi.blue, roi.times, roi.fsEff);

  size_t len = std::min({resR.values.size(), resG.values.size(), resB.values.size()});
  if (len < static_cast<size_t>(roi.fsEff * m_config.minWindowSeconds)) {
    return roi;
  }

  roi.resRed.assign(resR.values.begin(), resR.values.begin() + len);
  roi.resGreen.assign(resG.values.begin(), resG.values.begin() + len);
  roi.resBlue.assign(resB.values.begin(), resB.values.begin() + len);
  roi.resTimes.assign(resR.times.begin(), resR.times.begin() + len);

  roi.resRGB.resize(len);
  for (size_t i = 0; i < len; ++i) {
    roi.resRGB[i] = {roi.resRed[i], roi.resGreen[i], roi.resBlue[i]};
  }

  roi.sampleCount = static_cast<int>(len);
  roi.valid = true;
  return roi;
}

// ============================================================================
// Stage 2: Run Extractors
// ============================================================================

std::vector<ExtractorOutput>
VOREngine::runExtractors(const ROISignalBuffer &roi) {
  std::vector<ExtractorOutput> outputs;

  // Green
  {
    ExtractorOutput out;
    out.extractorName = "green";
    out.roiName = roi.roiName;
    auto sig = Extractors::green(roi.resGreen, roi.fsEff);
    if (!sig.empty()) {
      out.signal = std::move(sig);
      out.valid = true;
    }
    outputs.push_back(std::move(out));
  }

  // CHROM
  {
    ExtractorOutput out;
    out.extractorName = "chrom";
    out.roiName = roi.roiName;
    auto sig = Extractors::chrom(roi.resRGB, roi.fsEff);
    if (!sig.empty()) {
      out.signal = std::move(sig);
      out.valid = true;
    }
    outputs.push_back(std::move(out));
  }

  // POS
  {
    ExtractorOutput out;
    out.extractorName = "pos";
    out.roiName = roi.roiName;
    auto sig = Extractors::pos(roi.resRGB, roi.fsEff);
    if (!sig.empty()) {
      out.signal = std::move(sig);
      out.valid = true;
    }
    outputs.push_back(std::move(out));
  }

  return outputs;
}

// ============================================================================
// Stage 3: Spectrum Analysis
// ============================================================================

std::vector<SpectralPeak>
VOREngine::analyzeSpectrum(const std::vector<double> &signal, double fs,
                           double fMin, double fMax, int topK) {
  std::vector<SpectralPeak> peaks;
  if (signal.size() < 16) return peaks;

  int minSegLen = static_cast<int>(fs * 4.0);
  int nperseg = std::max(minSegLen, std::min(512, static_cast<int>(signal.size())));
  
  m_lastPsdResolution = (fs / nperseg) * 60.0;
  auto psd = Math::welch(signal, fs, nperseg, nperseg / 2);

  if (!psd.valid || psd.frequencies.empty()) return peaks;

  int iMin = -1, iMax = -1;
  for (size_t i = 0; i < psd.frequencies.size(); ++i) {
    if (psd.frequencies[i] >= fMin && iMin < 0) iMin = static_cast<int>(i);
    if (psd.frequencies[i] <= fMax) iMax = static_cast<int>(i);
  }
  if (iMin < 0 || iMax <= iMin) return peaks;

  std::vector<double> bandPower;
  for (int i = iMin; i <= iMax; ++i) bandPower.push_back(psd.power[i]);
  double noiseFloor = Math::median(bandPower);
  if (noiseFloor < 1e-20) noiseFloor = 1e-20;

  struct RawPeak {
    int idx; double hz; double power; double snr; double rho_b; double rho_h;
  };
  std::vector<RawPeak> rawPeaks;

  for (int i = iMin + 1; i < iMax; ++i) {
    if (psd.power[i] > psd.power[i - 1] && psd.power[i] > psd.power[i + 1] && psd.power[i] > noiseFloor * 1.5) {
      double snr = 10.0 * std::log10(psd.power[i] / noiseFloor);
      rawPeaks.push_back({i, psd.frequencies[i], psd.power[i], snr, 0.0, 0.0});
    }
  }

  // Risk logic
  double delta_b = m_config.boundaryPenaltyWidth;
  double epsilon_h = 0.12;

  for (auto &p : rawPeaks) {
    double distToEdge = std::min(p.hz - fMin, fMax - p.hz);
    p.rho_b = std::max(0.0, 1.0 - distToEdge / delta_b);
    p.rho_h = 0.0;
    for (const auto &other : rawPeaks) {
      if (std::abs(other.hz - p.hz) < 0.01) continue;
      double ratio = p.hz / other.hz;
      if (std::abs(ratio - 2.0) < epsilon_h && other.power > p.power * 0.3) {
        p.rho_h = std::max(p.rho_h, std::min(1.0, other.power / p.power));
      }
    }
  }

  std::sort(rawPeaks.begin(), rawPeaks.end(), [](const RawPeak &a, const RawPeak &b) {
    return a.power * (1.0 - a.rho_b) * (1.0 - a.rho_h) > b.power * (1.0 - b.rho_b) * (1.0 - b.rho_h);
  });

  int k = std::min(topK, static_cast<int>(rawPeaks.size()));
  for (int i = 0; i < k; ++i) {
    SpectralPeak sp;
    sp.hz = rawPeaks[i].hz;
    sp.bpm = sp.hz * 60.0;
    sp.power = rawPeaks[i].power;
    sp.snr = rawPeaks[i].snr;
    sp.rank = i + 1;
    sp.rho_b = rawPeaks[i].rho_b;
    sp.rho_h = rawPeaks[i].rho_h;
    peaks.push_back(sp);
  }
  return peaks;
}

// ============================================================================
// Stage 4: Build Candidates
// ============================================================================

std::vector<VORCandidate>
VOREngine::buildCandidates(const std::vector<ExtractorOutput> &outputs, double fs) {
  std::vector<VORCandidate> candidates;
  for (const auto &ext : outputs) {
    if (!ext.valid) continue;
    auto peaks = analyzeSpectrum(ext.signal, fs, m_config.fMinHz, m_config.fMaxHz, m_config.topKPeaks);
    for (const auto &peak : peaks) {
      if (peak.snr < m_config.minPeakSNR) continue;
      VORCandidate c;
      c.hz = peak.hz; c.bpm = peak.bpm; c.power = peak.power; c.snr = peak.snr;
      c.source = ext.extractorName + "_psd"; c.roi = ext.roiName; c.peakRank = peak.rank;
      c.boundaryRisk = peak.rho_b; c.harmonicRisk = peak.rho_h;
      double x = (peak.snr - 3.0) / 2.0;
      c.confidence = (1.0 / (1.0 + std::exp(-x))) * (1.0 - peak.rho_b) * (1.0 - peak.rho_h) / std::sqrt(peak.rank);
      candidates.push_back(c);
    }

    auto timePeaks = Math::findPeaks(ext.signal, fs, static_cast<int>(fs * 0.35), 0.02);
    timePeaks = Math::filterPeaksMAD(timePeaks, 3.0);
    if (timePeaks.size() >= 3) {
      std::vector<double> rrs;
      for (size_t i = 1; i < timePeaks.size(); ++i) rrs.push_back(timePeaks[i].time - timePeaks[i-1].time);
      double medRR = Math::median(rrs);
      if (medRR > 0.0) {
        VORCandidate c;
        c.bpm = 60.0 / medRR; c.hz = c.bpm / 60.0; c.source = ext.extractorName + "_peaks";
        c.roi = ext.roiName; c.confidence = std::min(1.0, static_cast<double>(timePeaks.size()) / 10.0);
        candidates.push_back(c);
      }
    }
  }
  return candidates;
}

// ============================================================================
// Stage 5: Cluster Fusion
// ============================================================================

std::vector<VORCluster>
VOREngine::clusterCandidates(const std::vector<VORCandidate> &candidates) {
  std::vector<VORCluster> clusters;
  if (candidates.empty()) return clusters;

  std::vector<int> sortedIdx(candidates.size());
  std::iota(sortedIdx.begin(), sortedIdx.end(), 0);
  std::sort(sortedIdx.begin(), sortedIdx.end(), [&](int a, int b) { return candidates[a].bpm < candidates[b].bpm; });

  for (int idx : sortedIdx) {
    const auto &c = candidates[idx];
    bool merged = false;
    for (auto &cl : clusters) {
      if (std::abs(c.bpm - cl.centerBPM) <= m_config.clusterRadiusBPM) {
        cl.memberIdx.push_back(idx);
        cl.memberCount++;
        cl.totalWeight += c.confidence;
        double wSum = 0.0, bSum = 0.0;
        for (int mi : cl.memberIdx) { bSum += candidates[mi].bpm * candidates[mi].confidence; wSum += candidates[mi].confidence; }
        cl.centerBPM = (wSum > 0.0) ? bSum / wSum : c.bpm;
        merged = true; break;
      }
    }
    if (!merged) {
      VORCluster n; n.centerBPM = c.bpm; n.memberCount = 1; n.totalWeight = c.confidence; n.memberIdx.push_back(idx);
      clusters.push_back(n);
    }
  }

  for (auto &cl : clusters) {
    cl.score = cl.totalWeight * std::log(1.0 + cl.memberCount);
  }
  std::sort(clusters.begin(), clusters.end(), [](const VORCluster &a, const VORCluster &b) { return a.score > b.score; });
  return clusters;
}

// ============================================================================
// Stage 6: Decide Trust & Temporal State
// ============================================================================

VORResult VOREngine::decideTrust(const std::vector<VORCluster> &clusters, const std::vector<VORCandidate> &candidates, double fsEff) {
  VORResult r; r.fsEff = fsEff; r.nCandidates = candidates.size(); r.nClusters = clusters.size();
  if (clusters.empty()) return r;
  const auto &best = clusters[0];
  r.bpm = Math::clip(best.centerBPM, m_config.fMinHz * 60.0, m_config.fMaxHz * 60.0);
  r.confidence = Math::clip(best.totalWeight / 2.0, 0.0, 1.0);
  r.decision = (r.confidence > m_config.acceptThreshold) ? VORDecision::ACCEPT : VORDecision::ABSTAIN;
  r.valid = (r.decision != VORDecision::ABSTAIN);
  return r;
}

void VOREngine::applyTemporalState(VORResult &result) {
  result.psdResolutionBPM = m_lastPsdResolution;
  double T = 1.0;
  if (!m_recentWinners.empty() && result.bpm > 0.0) {
    double m = result.bpm;
    for (double v : m_recentWinners) m += v;
    m /= (m_recentWinners.size() + 1);
    double var = (result.bpm - m) * (result.bpm - m);
    for (double v : m_recentWinners) var += (v - m) * (v - m);
    double sd = std::sqrt(var / (m_recentWinners.size() + 1));
    T = std::max(0.0, 1.0 - sd / 15.0);
  }
  result.temporalT = T;
  result.calibratedConf = result.confidence * std::sqrt(T);
  
  m_trustAge++;
  if (result.calibratedConf > 0.5) {
    m_trustedHR = result.bpm; m_trustAge = 0;
    result.decision = VORDecision::ACCEPT;
  } else if (m_trustAge < MAX_HOLD_CYCLES && m_trustedHR > 0.0) {
    result.bpm = m_trustedHR; result.decision = VORDecision::HOLD;
  }
  result.trustedHR = m_trustedHR; result.trustAge = m_trustAge;
}

} // namespace VOR
