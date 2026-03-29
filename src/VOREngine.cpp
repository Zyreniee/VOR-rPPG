#include "VOREngine.hpp"
#include "AdvancedSignalProcessing.hpp"
#include "SignalProcessor.hpp"
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

  // Stage 3: Build candidates from all extractors across all ROIs
  auto candidates = buildCandidates(allOutputs, fsEff);

  // Stage 4: Cluster candidates
  auto clusters = clusterCandidates(candidates);

  // Stage 5: Trust decision (per-cycle, memoryless)
  auto result = decideTrust(clusters, candidates, fsEff);

  // Save raw winning BPM before temporal state may override it
  double rawWinnerBPM = result.bpm;

  // Stage 6: Temporal state model — recalibrate confidence, HOLD logic
  applyTemporalState(result);

  // Update recent winners with raw cluster centroid (not HOLD-overridden)
  if (rawWinnerBPM > 0.0) {
    m_recentWinners.push_back(rawWinnerBPM);
    if (static_cast<int>(m_recentWinners.size()) > MAX_RECENT) {
      m_recentWinners.pop_front();
    }
  }

  // Fill per-extractor diagnostics
  for (const auto &c : candidates) {
    if (c.peakRank == 1) {
      result.top1BPMs.push_back(c.bpm);
      // Keep the first one we see for the legacy dashboard columns
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

  // Logging
  if (m_config.enableCSVLog) {
    logCSV(result, candidates, clusters);
  }
  logState(result);
  if (m_config.enableConsoleLog) {
    logConsole(result);
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

  // ---- 12-second window slice ----
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
    roi.green[i] = rgb[srcIdx][1]; // Extract Green directly from RGB
    roi.red[i] = rgb[srcIdx][0];
    roi.blue[i] = rgb[srcIdx][2];
  }

  // Check minimum window
  double duration =
      roi.times.size() > 1 ? roi.times.back() - roi.times.front() : 0.0;
  if (duration < m_config.minWindowSeconds) {
    return roi;
  }

  // ---- Estimate effective fs ----
  roi.fsEff = SignalProcessing::estimateFs(roi.times, nominalFs);
  roi.fsEff = SignalProcessing::clip(roi.fsEff, 10.0, 60.0);

  // ---- Per-channel resampling ----
  // Extract individual R, G, B from rgb for proper resampling
  std::vector<double> rawR(roi.times.size()), rawG(roi.times.size()),
      rawB(roi.times.size());
  for (size_t i = 0; i < roi.times.size(); ++i) {
    size_t srcIdx = sliceIdx + i;
    rawR[i] = rgb[srcIdx][0]; // R
    rawG[i] = rgb[srcIdx][1]; // G
    rawB[i] = rgb[srcIdx][2]; // B
  }

  auto [rR, tR] = SignalProcessing::resampleSignal(rawR, roi.times, roi.fsEff);
  auto [rG, tG] = SignalProcessing::resampleSignal(rawG, roi.times, roi.fsEff);
  auto [rB, tB] = SignalProcessing::resampleSignal(rawB, roi.times, roi.fsEff);
  auto [rGreen, tGreen] =
      SignalProcessing::resampleSignal(roi.green, roi.times, roi.fsEff);

  // Use minimum common length
  size_t len = std::min({rR.size(), rG.size(), rB.size(), rGreen.size()});
  if (len < static_cast<size_t>(roi.fsEff * m_config.minWindowSeconds)) {
    return roi;
  }

  roi.resRed.assign(rR.begin(), rR.begin() + len);
  roi.resGreen.assign(rGreen.begin(), rGreen.begin() + len);
  roi.resBlue.assign(rB.begin(), rB.begin() + len);
  roi.resTimes.assign(tR.begin(), tR.begin() + len);

  // Build resampled RGB triplet
  roi.resRGB.resize(len);
  for (size_t i = 0; i < len; ++i) {
    roi.resRGB[i] = {rR[i], rG[i], rB[i]};
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

  // ---- Green channel ----
  {
    ExtractorOutput out;
    out.extractorName = "green";
    out.roiName = roi.roiName;
    auto sig = SignalProcessing::detrend(roi.resGreen);
    sig = SignalProcessing::normalizeSignal(sig);
    if (!sig.empty()) {
      out.signal = std::move(sig);
      out.valid = true;
    }
    outputs.push_back(std::move(out));
  }

  // ---- CHROM ----
  {
    ExtractorOutput out;
    out.extractorName = "chrom";
    out.roiName = roi.roiName;
    // chromRPPG now only does detrend+normalize (bandpass removed in Stage 1 fix)
    auto sig = SignalProcessing::chromRPPG(roi.resRGB, roi.fsEff);
    if (!sig.empty()) {
      out.signal = std::move(sig);
      out.valid = true;
    }
    outputs.push_back(std::move(out));
  }

  // ---- POS ----
  {
    ExtractorOutput out;
    out.extractorName = "pos";
    out.roiName = roi.roiName;
    auto sig = SignalProcessing::posRPPG(roi.resRGB, roi.fsEff);
    if (!sig.empty()) {
      out.signal = std::move(sig);
      out.valid = true;
    }
    outputs.push_back(std::move(out));
  }

  // ---- VOR adaptive projection (Stage 2 flag) ----
  // TODO: Implement when enableVORProjection is true

  return outputs;
}

// ============================================================================
// Stage 3: Spectrum Analysis (top-k peaks)
// ============================================================================

std::vector<SpectralPeak>
VOREngine::analyzeSpectrum(const std::vector<double> &signal, double fs,
                           double fMin, double fMax, int topK) {
  std::vector<SpectralPeak> peaks;

  if (signal.size() < 16)
    return peaks;

  // Compute Welch PSD — use at least 4s of data per segment for
  // adequate frequency resolution. With fs=15 and nperseg=60 (old: N/3),
  // resolution was 0.25 Hz = 15 BPM. Now with nperseg=max(fs*4, min(512,N))
  // resolution is ~0.06 Hz = ~3.7 BPM at worst.
  int minSegLen = static_cast<int>(fs * 4.0);
  int nperseg = std::max(minSegLen,
                         std::min(512, static_cast<int>(signal.size())));
  nperseg = std::min(nperseg, static_cast<int>(signal.size()));
  m_lastPsdResolution = (fs / nperseg) * 60.0; // BPM per bin
  auto psd = SignalProcessing::welchPSD(signal, fs, nperseg, nperseg / 2);

  if (!psd.valid || psd.frequencies.empty())
    return peaks;

  // Find band indices
  int iMin = -1, iMax = -1;
  for (size_t i = 0; i < psd.frequencies.size(); ++i) {
    if (psd.frequencies[i] >= fMin && iMin < 0)
      iMin = static_cast<int>(i);
    if (psd.frequencies[i] <= fMax)
      iMax = static_cast<int>(i);
  }

  if (iMin < 0 || iMax <= iMin)
    return peaks;

  // Compute noise floor (median of band power)
  std::vector<double> bandPower;
  for (int i = iMin; i <= iMax; ++i) {
    bandPower.push_back(psd.power[i]);
  }
  double noiseFloor = SignalProcessing::median(bandPower);
  if (noiseFloor < 1e-20)
    noiseFloor = 1e-20;

  // Find ALL local maxima in band
  struct RawPeak {
    int idx;
    double hz;
    double power;
    double snr;
    double rho_b;
    double rho_h;
  };
  std::vector<RawPeak> rawPeaks;

  for (int i = iMin + 1; i < iMax; ++i) {
    if (psd.power[i] > psd.power[i - 1] && psd.power[i] > psd.power[i + 1] &&
        psd.power[i] > noiseFloor * 1.5) {
      
      // Sub-bin parabolic interpolation
      double y1 = psd.power[i - 1];
      double y2 = psd.power[i];
      double y3 = psd.power[i + 1];
      double denom = (y1 - 2 * y2 + y3);
      double delta = 0.0;
      if (std::abs(denom) > 1e-12) {
        delta = 0.5 * (y1 - y3) / denom;
      }
      double df = psd.frequencies[i] - psd.frequencies[i - 1];
      double subBinHz = psd.frequencies[i] + delta * df;

      double snr = 10.0 * std::log10(y2 / noiseFloor);
      rawPeaks.push_back({i, subBinHz, y2, snr, 0.0, 0.0});
    }
  }

  // Also check the band edges (no exact mathematical neighbors for interpolation, just use bin freq)
  if (iMin < iMax && psd.power[iMin] > psd.power[iMin + 1] &&
      psd.power[iMin] > noiseFloor * 1.5) {
    double snr = 10.0 * std::log10(psd.power[iMin] / noiseFloor);
    rawPeaks.push_back({iMin, psd.frequencies[iMin], psd.power[iMin], snr, 0.0, 0.0});
  }

  // If no local maxima found, use global max in band
  if (rawPeaks.empty()) {
    int bestIdx = iMin;
    for (int i = iMin + 1; i <= iMax; ++i) {
      if (psd.power[i] > psd.power[bestIdx])
        bestIdx = i;
    }
    
    // Attempt parabolic interpolation if strictly inside bounds
    double subBinHz = psd.frequencies[bestIdx];
    if (bestIdx > 0 && bestIdx < static_cast<int>(psd.power.size()) - 1) {
      double y1 = psd.power[bestIdx - 1];
      double y2 = psd.power[bestIdx];
      double y3 = psd.power[bestIdx + 1];
      double denom = (y1 - 2 * y2 + y3);
      if (std::abs(denom) > 1e-12) {
        double delta = 0.5 * (y1 - y3) / denom;
        double df = psd.frequencies[bestIdx] - psd.frequencies[bestIdx - 1];
        subBinHz += delta * df;
      }
    }

    double snr = 10.0 * std::log10(psd.power[bestIdx] / noiseFloor);
    rawPeaks.push_back(
        {bestIdx, subBinHz, psd.power[bestIdx], snr, 0.0, 0.0});
  }

  // ---- Compute risk annotations (VOR spec §3.5) ----
  double delta_b = m_config.boundaryPenaltyWidth; // Hz from band edge
  double epsilon_h = 0.12; // relative tolerance for harmonic ratio detection

  for (auto &p : rawPeaks) {
    // Boundary risk: ρ_b = max(0, 1 - min(f-fMin, fMax-f) / Δ_b)
    double distToEdge = std::min(p.hz - fMin, fMax - p.hz);
    p.rho_b = std::max(0.0, 1.0 - distToEdge / delta_b);
    p.rho_b = std::min(1.0, p.rho_b);

    // Harmonic risk: ρ_h = max over other peaks of
    //   1(|f_j/f_i - 2| < ε) * (P(f_i)/P(f_j))  [if f_i is the fundamental]
    //   or 1(|f_j/f_i - 0.5| < ε) * (P(f_i)/P(f_j)) [if f_j is subharmonic]
    p.rho_h = 0.0;
    for (const auto &other : rawPeaks) {
      if (std::abs(other.hz - p.hz) < 0.01) continue; // skip self
      if (other.hz < 0.01) continue;
      double ratio = p.hz / other.hz;
      // Check if p is a harmonic (2×) of other
      if (std::abs(ratio - 2.0) < epsilon_h && other.power > p.power * 0.3) {
        double risk = other.power / std::max(p.power, 1e-20);
        p.rho_h = std::max(p.rho_h, std::min(1.0, risk));
      }
      // Check if p is a subharmonic (½×) of other
      if (std::abs(ratio - 0.5) < epsilon_h && other.power > p.power * 0.3) {
        double risk = other.power / std::max(p.power, 1e-20);
        p.rho_h = std::max(p.rho_h, std::min(1.0, risk));
      }
      // Check if p is a 3rd harmonic of other
      if (std::abs(ratio - 3.0) < epsilon_h && other.power > p.power * 0.3) {
        double risk = other.power / std::max(p.power, 1e-20);
        p.rho_h = std::max(p.rho_h, std::min(1.0, risk * 0.7));
      }
    }
  }

  // Sort by power descending
  std::sort(rawPeaks.begin(), rawPeaks.end(),
            [](const RawPeak &a, const RawPeak &b) {
              // Sort by risk-adjusted power: power * (1 - rho_b) * (1 - rho_h)
              double sa = a.power * (1.0 - a.rho_b * 0.5) * (1.0 - a.rho_h * 0.5);
              double sb = b.power * (1.0 - b.rho_b * 0.5) * (1.0 - b.rho_h * 0.5);
              return sa > sb;
            });

  // Take top-k
  int k = std::min(topK, static_cast<int>(rawPeaks.size()));
  for (int i = 0; i < k; ++i) {
    SpectralPeak sp;
    sp.hz = rawPeaks[i].hz;
    sp.bpm = rawPeaks[i].hz * 60.0;
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
VOREngine::buildCandidates(const std::vector<ExtractorOutput> &outputs,
                           double fs) {
  std::vector<VORCandidate> candidates;

  for (const auto &ext : outputs) {
    if (!ext.valid || ext.signal.empty())
      continue;

    // PSD-based candidates (top-k peaks)
    auto peaks =
        analyzeSpectrum(ext.signal, fs, m_config.fMinHz, m_config.fMaxHz,
                        m_config.topKPeaks);

    for (const auto &peak : peaks) {
      if (peak.snr < m_config.minPeakSNR)
        continue;

      VORCandidate c;
      c.hz = peak.hz;
      c.bpm = peak.bpm;
      c.power = peak.power;
      c.snr = peak.snr;
      c.source = ext.extractorName + "_psd";
      c.roi = ext.roiName;
      c.peakRank = peak.rank;
      c.boundaryRisk = peak.rho_b;
      c.harmonicRisk = peak.rho_h;

      // VOR spec §3.5: φ = sigmoid(SNR - 3) * (1 - ρ_b) * (1 - ρ_h) / √(rank)
      double x = (peak.snr - 3.0) / 2.0;
      c.confidence = 1.0 / (1.0 + std::exp(-x));
      c.confidence *= (1.0 - peak.rho_b);
      c.confidence *= (1.0 - peak.rho_h);
      c.confidence *= (1.0 / std::sqrt(static_cast<double>(peak.rank)));

      candidates.push_back(c);
    }

    // Time-domain peak-based candidates
    {
      auto timePeaks = SignalProcessing::findPeaks(
          ext.signal, fs, static_cast<int>(fs * 0.35), 0.02);
      timePeaks = SignalProcessing::filterPeaksMAD(timePeaks, 3.0);

      if (timePeaks.size() >= 3) {
        std::vector<double> rrIntervals;
        for (size_t i = 1; i < timePeaks.size(); ++i) {
          double dt = timePeaks[i].time - timePeaks[i - 1].time;
          if (dt > 0.0)
            rrIntervals.push_back(dt);
        }

        if (!rrIntervals.empty()) {
          double medianRR = SignalProcessing::median(rrIntervals);
          if (medianRR > 0.0) {
            double bpm = 60.0 / medianRR;
            bpm = SignalProcessing::clip(bpm, m_config.fMinHz * 60.0,
                                         m_config.fMaxHz * 60.0);

            VORCandidate c;
            c.hz = bpm / 60.0;
            c.bpm = bpm;
            c.power = 0.0;
            c.snr = 0.0; // unknown for time-domain
            c.source = ext.extractorName + "_peaks";
            c.roi = ext.roiName;
            c.peakRank = 1;
            c.confidence = std::min(
                1.0, static_cast<double>(timePeaks.size()) / 8.0); // more peaks = more trust

            candidates.push_back(c);
          }
        }
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

  if (candidates.empty())
    return clusters;

  // Sort candidates by BPM
  std::vector<int> sortedIdx(candidates.size());
  std::iota(sortedIdx.begin(), sortedIdx.end(), 0);
  std::sort(sortedIdx.begin(), sortedIdx.end(),
            [&candidates](int a, int b) {
              return candidates[a].bpm < candidates[b].bpm;
            });

  // Greedy clustering
  for (int idx : sortedIdx) {
    const auto &c = candidates[idx];
    bool merged = false;

    for (auto &cl : clusters) {
      if (std::abs(c.bpm - cl.centerBPM) <= m_config.clusterRadiusBPM) {
        // Merge into cluster
        cl.memberIdx.push_back(idx);
        cl.memberCount++;
        cl.totalWeight += c.confidence;

        // Update weighted centroid
        double weightSum = 0.0;
        double bpmSum = 0.0;
        for (int mi : cl.memberIdx) {
          bpmSum += candidates[mi].bpm * candidates[mi].confidence;
          weightSum += candidates[mi].confidence;
        }
        cl.centerBPM = (weightSum > 0.0) ? bpmSum / weightSum : c.bpm;

        // Track sources
        if (c.source.find("psd") != std::string::npos)
          cl.hasPSD = true;
        if (c.source.find("peaks") != std::string::npos)
          cl.hasPeaks = true;

        // Update average SNR
        double snrSum = 0.0;
        int snrCount = 0;
        for (int mi : cl.memberIdx) {
          if (candidates[mi].snr > 0.0) {
            snrSum += candidates[mi].snr;
            snrCount++;
          }
        }
        cl.avgSNR = (snrCount > 0) ? snrSum / snrCount : 0.0;

        merged = true;
        break;
      }
    }

    if (!merged) {
      VORCluster newCl;
      newCl.centerBPM = c.bpm;
      newCl.memberCount = 1;
      newCl.totalWeight = c.confidence;
      newCl.memberIdx.push_back(idx);
      newCl.avgSNR = c.snr;
      if (c.source.find("psd") != std::string::npos)
        newCl.hasPSD = true;
      if (c.source.find("peaks") != std::string::npos)
        newCl.hasPeaks = true;
      clusters.push_back(newCl);
    }
  }

  // Compute spread and score for each cluster
  for (auto &cl : clusters) {
    if (cl.memberCount >= 2) {
      double sumSq = 0.0;
      for (int mi : cl.memberIdx) {
        double diff = candidates[mi].bpm - cl.centerBPM;
        sumSq += diff * diff;
      }
      cl.spreadBPM = std::sqrt(sumSq / cl.memberCount);
    } else {
      cl.spreadBPM = m_config.clusterRadiusBPM; // single member = max spread
    }

    // Score: weight * diversity * tightness
    double diversityBonus = (cl.hasPSD && cl.hasPeaks) ? 1.5 : 1.0;
    double tightness = 1.0 / (1.0 + cl.spreadBPM);
    cl.score = cl.totalWeight * std::log(1.0 + cl.memberCount) * tightness *
               diversityBonus;
  }

  // Sort clusters by score descending
  std::sort(clusters.begin(), clusters.end(),
            [](const VORCluster &a, const VORCluster &b) {
              return a.score > b.score;
            });

  return clusters;
}

// ============================================================================
// Stage 6: Trust Engine
// ============================================================================

VORResult
VOREngine::decideTrust(const std::vector<VORCluster> &clusters,
                       [[maybe_unused]] const std::vector<VORCandidate> &candidates,
                       double fsEff) {
  VORResult r;
  r.fsEff = fsEff;
  r.nCandidates = static_cast<int>(candidates.size());
  r.nClusters = static_cast<int>(clusters.size());

  if (clusters.empty()) {
    r.decision = VORDecision::ABSTAIN;
    r.confidence = 0.0;
    return r;
  }

  const auto &best = clusters[0];

  // ---- Confidence components ----

  // C1: Cluster size (4+ members = max)
  double cCluster = std::min(1.0, static_cast<double>(best.memberCount) / 4.0);

  // C2: Tightness (spread < 3 BPM = good)
  double cSpread = std::max(0.0, 1.0 - best.spreadBPM / 10.0);

  // C3: SNR (above 3 dB = good)
  double xSnr = (best.avgSNR - 3.0) / 3.0;
  double cSNR = 1.0 / (1.0 + std::exp(-xSnr));

  // C4: Source diversity (both PSD and peaks = bonus)
  double cDiversity = (best.hasPSD && best.hasPeaks) ? 1.0 : 0.5;

  // C5: Separation from runner-up (if exists)
  double cSeparation = 1.0;
  if (clusters.size() >= 2) {
    double gap = std::abs(best.centerBPM - clusters[1].centerBPM);
    cSeparation =
        std::min(1.0, gap / 15.0); // 15+ BPM separation = confident
  }

  // Weighted combination
  r.confidence = 0.25 * cCluster + 0.20 * cSpread + 0.25 * cSNR +
                 0.15 * cDiversity + 0.15 * cSeparation;

  r.bpm = best.centerBPM;
  r.winningClusterSize = best.memberCount;
  r.winningClusterSpread = best.spreadBPM;

  // ---- Decision ----
  if (r.confidence >= m_config.acceptThreshold) {
    r.decision = VORDecision::ACCEPT;
    r.valid = true;
  } else if (r.confidence >= m_config.degradedThreshold) {
    r.decision = VORDecision::DEGRADED;
    r.valid = true;
  } else {
    r.decision = VORDecision::ABSTAIN;
    r.valid = false;
  }

  // Physiological bounds
  r.bpm = SignalProcessing::clip(r.bpm, m_config.fMinHz * 60.0,
                                  m_config.fMaxHz * 60.0);

  return r;
}

// ============================================================================
// Stage 7: Temporal State Model (Patches 2, 3, 5)
// ============================================================================

void VOREngine::applyTemporalState(VORResult &result) {
  // ---- Record PSD resolution ----
  result.psdResolutionBPM = m_lastPsdResolution;

  // ---- Compute temporal consistency T ----
  // T = 1 - std(recent_winners ∪ {raw_bpm}) / σ_max
  double T = 1.0;
  if (result.bpm > 0.0) {
    std::vector<double> window(m_recentWinners.begin(), m_recentWinners.end());
    window.push_back(result.bpm);

    if (window.size() >= 3) {
      double mean = 0.0;
      for (double v : window) mean += v;
      mean /= window.size();

      double variance = 0.0;
      for (double v : window) variance += (v - mean) * (v - mean);
      double stddev = std::sqrt(variance / window.size());

      constexpr double SIGMA_MAX = 15.0; // BPM
      T = std::max(0.0, 1.0 - stddev / SIGMA_MAX);
    }
  } else {
    T = 0.0; // No BPM → no temporal consistency
  }
  result.temporalT = T;

  // ---- Compute method agreement M across all ROIs and algorithms ----
  // M = 1 - range(all top-1 candidates) / Δ_max
  double M = 1.0;
  {
    if (result.top1BPMs.size() >= 2) {
      double maxVal = *std::max_element(result.top1BPMs.begin(), result.top1BPMs.end());
      double minVal = *std::min_element(result.top1BPMs.begin(), result.top1BPMs.end());
      double range = maxVal - minVal;
      constexpr double DELTA_MAX = 20.0; // BPM
      M = std::max(0.0, 1.0 - range / DELTA_MAX);
    }
  }
  result.methodM = M;

  // ---- Compute resolution quality R ----
  // R = min(1, nperseg / (fs * 4)) — already captured as psdResolutionBPM
  // Good resolution: <4 BPM/bin → R close to 1
  // Poor resolution: >10 BPM/bin → R close to 0
  double R = 1.0;
  if (m_lastPsdResolution > 0.0) {
    R = std::max(0.0, std::min(1.0, 4.0 / m_lastPsdResolution));
  }
  result.resolutionR = R;

  // ---- Recalibrate confidence ----
  // calibrated = raw × T^0.5 × M^0.5 × R
  double calibrated = result.confidence
                      * std::sqrt(T)
                      * std::sqrt(M)
                      * R;
  calibrated = std::max(0.0, std::min(1.0, calibrated));
  result.calibratedConf = calibrated;

  // ---- State transition ----
  constexpr double THETA_ACCEPT = 0.55;
  constexpr double THETA_DEGRADED = 0.30;
  constexpr double T_MIN_FOR_ACCEPT = 0.5;

  m_trustAge++;

  if (calibrated >= THETA_ACCEPT && T >= T_MIN_FOR_ACCEPT) {
    // Strong evidence, temporally consistent → ACCEPT
    result.decision = VORDecision::ACCEPT;
    result.valid = true;
    result.reasonCode = "ACCEPT_STABLE";

    // Update trusted state
    m_trustedHR = result.bpm;
    m_trustedConf = calibrated;
    m_trustAge = 0;

  } else if (m_trustedHR > 0.0 && m_trustAge < MAX_HOLD_CYCLES) {
    // Have recent trusted HR → HOLD
    result.decision = VORDecision::HOLD;
    result.bpm = m_trustedHR; // Output the trusted value
    result.valid = true;

    if (T < T_MIN_FOR_ACCEPT) {
      result.reasonCode = "HOLD_VOLATILE";
    } else {
      result.reasonCode = "HOLD_LOW_CONF";
    }

  } else if (calibrated >= THETA_DEGRADED) {
    // No trusted HR to hold, but some evidence → DEGRADED
    result.decision = VORDecision::DEGRADED;
    result.valid = true;
    result.reasonCode = "DEGRADED_NO_HISTORY";

  } else {
    // Nothing reliable
    result.decision = VORDecision::ABSTAIN;
    result.valid = false;
    result.bpm = 0.0;
    result.reasonCode = "ABSTAIN_LOW_CONF";
  }

  // Fill trusted state diagnostics
  result.trustedHR = m_trustedHR;
  result.trustAge = m_trustAge;

  // (Recent winners are pushed in compute() with the raw BPM)
}

// ============================================================================
// Logging
// ============================================================================

void VOREngine::logCSV(const VORResult &result,
                       const std::vector<VORCandidate> &candidates,
                       const std::vector<VORCluster> &clusters) {
  static std::ofstream csvFile;
  if (!m_csvHeaderWritten) {
    csvFile.open("vor_debug.csv", std::ios::trunc);
    if (csvFile.is_open()) {
      csvFile << "cycle,fs_eff,n_cands,n_clusters,"
              << "chrom_bpm,chrom_snr,pos_bpm,pos_snr,green_bpm,green_snr,"
              << "cl1_bpm,cl1_size,cl1_spread,cl1_score,"
              << "cl2_bpm,cl2_size,cl2_spread,cl2_score,"
              << "winning_bpm,confidence,decision" << std::endl;
      m_csvHeaderWritten = true;
    }
  }

  if (!csvFile.is_open())
    return;

  auto decStr = [](VORDecision d) -> std::string {
    switch (d) {
    case VORDecision::ACCEPT:
      return "ACCEPT";
    case VORDecision::HOLD:
      return "HOLD";
    case VORDecision::DEGRADED:
      return "DEGRADED";
    case VORDecision::ABSTAIN:
      return "ABSTAIN";
    }
    return "UNKNOWN";
  };

  csvFile << m_cycleCount << "," << result.fsEff << ","
          << result.nCandidates << "," << result.nClusters << ","
          << result.chromBPM << "," << result.chromSNR << ","
          << result.posBPM << "," << result.posSNR << ","
          << result.greenBPM << "," << result.greenSNR << ",";

  // Cluster 1
  if (!clusters.empty()) {
    csvFile << clusters[0].centerBPM << "," << clusters[0].memberCount << ","
            << clusters[0].spreadBPM << "," << clusters[0].score;
  } else {
    csvFile << "0,0,0,0";
  }
  csvFile << ",";

  // Cluster 2
  if (clusters.size() >= 2) {
    csvFile << clusters[1].centerBPM << "," << clusters[1].memberCount << ","
            << clusters[1].spreadBPM << "," << clusters[1].score;
  } else {
    csvFile << "0,0,0,0";
  }
  csvFile << ",";

  csvFile << result.bpm << "," << result.confidence << ","
          << decStr(result.decision) << std::endl;

  csvFile.flush();

  // ---- Per-candidate detail log (for scientific evaluation) ----
  static std::ofstream candFile;
  static bool candHeaderWritten = false;
  if (!candHeaderWritten) {
    candFile.open("vor_candidates.csv", std::ios::trunc);
    if (candFile.is_open()) {
      candFile << "cycle,method,roi,type,bpm,snr,rank,rho_b,rho_h,phi"
               << std::endl;
      candHeaderWritten = true;
    }
  }
  if (candFile.is_open()) {
    for (const auto &c : candidates) {
      std::string ctype =
          (c.source.find("peaks") != std::string::npos) ? "temporal" : "spectral";
      candFile << m_cycleCount << "," << c.source << "," << c.roi << ","
               << ctype << "," << c.bpm << "," << c.snr << "," << c.peakRank
               << "," << c.boundaryRisk << "," << c.harmonicRisk << ","
               << c.confidence << std::endl;
    }
    candFile.flush();
  }
}

void VOREngine::logConsole(const VORResult &result) {
  auto decStr = [](VORDecision d) -> std::string {
    switch (d) {
    case VORDecision::ACCEPT:
      return "ACCEPT";
    case VORDecision::HOLD:
      return "HOLD";
    case VORDecision::DEGRADED:
      return "DEGRADED";
    case VORDecision::ABSTAIN:
      return "ABSTAIN";
    }
    return "?";
  };

  std::cout << "[VOR] cy=" << m_cycleCount
            << " hr=" << static_cast<int>(result.bpm)
            << " raw_conf=" << static_cast<int>(result.confidence * 100)
            << " cal_conf=" << static_cast<int>(result.calibratedConf * 100)
            << " T=" << static_cast<int>(result.temporalT * 100)
            << " M=" << static_cast<int>(result.methodM * 100)
            << " R=" << static_cast<int>(result.resolutionR * 100)
            << " " << decStr(result.decision)
            << " (" << result.reasonCode << ")"
            << std::endl;
}

void VOREngine::logState(const VORResult &result) {
  static std::ofstream stateFile;
  static bool stateHeaderWritten = false;
  if (!stateHeaderWritten) {
    stateFile.open("vor_state.csv", std::ios::trunc);
    if (stateFile.is_open()) {
      stateFile << "cycle,raw_bpm,final_bpm,raw_conf,cal_conf,"
                << "temporal_T,method_M,resolution_R,"
                << "decision,reason,trusted_hr,trust_age,"
                << "chrom_top1,pos_top1,green_top1,"
                << "n_candidates,n_clusters,"
                << "psd_resolution_bpm,cluster_size,cluster_spread"
                << std::endl;
      stateHeaderWritten = true;
    }
  }
  if (!stateFile.is_open()) return;

  auto decStr = [](VORDecision d) -> std::string {
    switch (d) {
    case VORDecision::ACCEPT:   return "ACCEPT";
    case VORDecision::HOLD:     return "HOLD";
    case VORDecision::DEGRADED: return "DEGRADED";
    case VORDecision::ABSTAIN:  return "ABSTAIN";
    }
    return "?";
  };

  // Raw BPM = last entry in recent winners (if any)
  double rawBpm = m_recentWinners.empty() ? 0.0 : m_recentWinners.back();

  stateFile << m_cycleCount << ","
            << rawBpm << ","
            << result.bpm << ","
            << result.confidence << ","
            << result.calibratedConf << ","
            << result.temporalT << ","
            << result.methodM << ","
            << result.resolutionR << ","
            << decStr(result.decision) << ","
            << result.reasonCode << ","
            << result.trustedHR << ","
            << result.trustAge << ","
            << result.chromBPM << ","
            << result.posBPM << ","
            << result.greenBPM << ","
            << result.nCandidates << ","
            << result.nClusters << ","
            << result.psdResolutionBPM << ","
            << result.winningClusterSize << ","
            << result.winningClusterSpread
            << std::endl;
  stateFile.flush();
}

} // namespace VOR
